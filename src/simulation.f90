!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use dev_lpt_class,     only: lpt
   use ddadi_class,       only: ddadi
   use multivdscalar_class, only: multivdscalar
   use vdscalar_class,    only: vdscalar
   use hypre_uns_class,   only: hypre_uns
   use hypre_str_class,   only: hypre_str
   use lowmach_class,     only: lowmach
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use string, only: str_medium
   implicit none
   private

   !> Get an LPT solver, a lowmach solver, and corresponding time tracker, plus a couple of linear solvers
   type(hypre_uns),     public :: ps
   type(hypre_str),     public :: vs
   type(ddadi),         public :: mss,tps
   type(lowmach),       public :: fs
   type(lpt),           public :: lp
   type(multivdscalar), public :: msc  !< Mass fraction (Composition) solver
   type(vdscalar),      public :: tp   !< Temperature solver
   type(timetracker),   public :: time

   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,lptfile,tfile,mscfile,outfluxfile

   public :: simulation_init,simulation_run,simulation_final

   !> Work arrays and fluid properties
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,resRHO
   real(WP), dimension(:,:,:), allocatable :: srcUlp,srcVlp,srcWlp
   real(WP) :: visc
   real(WP), dimension(:,:,:), allocatable :: srcConti      !> Continuity source term
   real(WP), dimension(:,:,:), allocatable :: rho0          !> density field

   !> Multiscalar solver work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resMSC,srcMSC,mscTmp    !< Residuals, source, temp solution for multiscalar solver
   logical, dimension(:,:,:,:), allocatable :: bqflag                   !< Flag for bquick scheme
   real(WP), dimension(:,:,:), allocatable :: massFracSum                          !< Sum of mass fractions for check
   real(WP), dimension(:), allocatable :: SCmassOutSum             !< Sum of mass fractions for check

   !> Temperature solver work arrays
   real(WP), dimension(:,:,:), allocatable :: resTp, srcTp               !< Residuals and source for temperature solver
   real(WP), dimension(:,:,:), allocatable :: gTfield, gCpfield, gLambdafield !< Temperature; Averaged Cp and Lambda fields


   !> Add thermal properties
   real(WP) :: initTemp
   !> gas heat capacity and conductivity array for computing mean values
   real(WP), dimension(:), allocatable :: gCp, gLambda
   real(WP) :: Le ! Lewis number

   !> gas density array for computing mean density
   real(WP), dimension(:), allocatable :: gRhos
   character(len=str_medium), dimension(:), allocatable :: SCname
   integer:: speciesNum

   !> inlet gas properties
   real(WP) :: inlet_velocity
   real(WP), dimension(:), allocatable :: inletComp    !> inlet composition
   real(WP) :: inletTemp

   !> Wallclock time for monitoring
   type :: timer
      real(WP) :: time_in
      real(WP) :: time
      real(WP) :: percent
   end type timer
   type(timer) :: wt_total,wt_vel,wt_pres,wt_lpt,wt_rest

contains


   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function left_of_domain

   !> Function that localizes the left (x-) of the domain for scalar
   function left_of_domainsc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1) isIn=.true.
   end function left_of_domainsc

   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain

   !> This subroutine is to get rho for multiscalar
   !> Note: here the density already multiplied by the volume fraction
   !> It will also be used to calculate velocity field / other scalars
   subroutine get_rho()
      implicit none
      integer :: i,j,k
      do k=msc%cfg%kmino_,msc%cfg%kmaxo_
         do j=msc%cfg%jmino_,msc%cfg%jmaxo_
            do i=msc%cfg%imino_,msc%cfg%imaxo_
               rho0(i,j,k)=dot_product(msc%SC(i,j,k,:),gRhos)
               msc%rho(i,j,k)=dot_product(msc%SC(i,j,k,:),gRhos)*(1-lp%VF(i,j,k))
            end do
         end do
      end do
   end subroutine get_rho

   !< This subroutine is to get averaged lambda/Cp
   !< It's related to the diffusion term of Temp transport equation
   subroutine get_tp_diff()
      implicit none
      integer :: i,j,k
      do k=tp%cfg%kmino_,tp%cfg%kmaxo_
         do j=tp%cfg%jmino_,tp%cfg%jmaxo_
            do i=tp%cfg%imino_,tp%cfg%imaxo_
               tp%diff(i,j,k) = dot_product(msc%SC(i,j,k,:),gLambda)/dot_product(msc%SC(i,j,k,:),gCp)
            end do
         end do
      end do
   end subroutine get_tp_diff

   !< This subroutine is to apply Lewis number to the diffusion term of multi-scalar
   subroutine applyLewisNumber()
      implicit none
      integer :: i,j,k
      do k=msc%cfg%kmino_,msc%cfg%kmaxo_
         do j=msc%cfg%jmino_,msc%cfg%jmaxo_
            do i=msc%cfg%imino_,msc%cfg%imaxo_
               if (Le .gt. 0.0_WP) msc%diff(i,j,k,:)=Le*tp%diff(i,j,k)
            end do
         end do
      end do
   end subroutine applyLewisNumber

   subroutine checkMassFractionSum()
      implicit none
      integer :: i,j,k
      do k=msc%cfg%kmino_,msc%cfg%kmaxo_
         do j=msc%cfg%jmino_,msc%cfg%jmaxo_
            do i=msc%cfg%imino_,msc%cfg%imaxo_
               massFracSum(i,j,k) = sum(msc%SC(i,j,k,:))
            end do
         end do
      end do
   end subroutine checkMassFractionSum

   subroutine getFluxOutSum(time,dt)
      use multivdscalar_class, only: bcond
      implicit none
      integer :: i,j,k,isc,n
      type(bcond), pointer :: mybc
      real(WP) :: dt,time
      real(WP) :: area
      area = fs%cfg%dy(1)*fs%cfg%dz(1)
      call msc%get_bcond('outflow',mybc)
      do n=1,mybc%itr%no_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         if (j.le.msc%cfg%jmax_ .and. j.ge.msc%cfg%jmin_ .and. k.le.msc%cfg%kmax_ .and. k.ge.msc%cfg%kmin_) then
            do isc=1,msc%nscalar
               SCmassOutSum(isc) = SCmassOutSum(isc) + msc%SC(i,j,k,isc)*fs%rhoU(i,j,k)*area*dt
            end do
         end if
      end do
   end subroutine getFluxOutSum

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      use string, only: str_medium
      implicit none

      call param_read('Species Numbers',speciesNum)

      ! read species properties
      allocate(gRhos(speciesNum))
      allocate(SCname(speciesNum))
      call param_read('Species gas phase densities',gRhos)
      call param_read('Species names',SCname)

      ! read inlet flow properties
      allocate(inletComp(speciesNum))
      call param_read('Inlet composition',inletComp)
      call param_read('Inlet temperature',inletTemp)

      ! read thermal property parameters
      allocate(gCp(speciesNum))
      allocate(gLambda(speciesNum))
      call param_read('Species gas phase heat capacities',gCp)
      call param_read('Species gas phase heat conductivities',gLambda)
      call param_read('Lewis number',Le, default=-1.0_WP)

      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Initialize timers
      initialize_timers: block
         wt_total%time=0.0_WP; wt_total%percent=0.0_WP
         wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
         wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
         wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
         wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP
      end block initialize_timers

      ! Create a low Mach flow solver with bconds
      create_flow_solver: block
         use hypre_uns_class, only: gmres_amg
         use hypre_str_class, only: pcg_pfmg
         use lowmach_class,   only: dirichlet,clipped_neumann
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Define boundary conditions
         call fs%add_bcond(name='inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )

         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Configure pressure solver
         ps=hypre_uns(cfg=cfg,name='Pressure',method=gmres_amg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)

         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver

      ! Create a multiscalar solver
      create_multiscalar: block
         use multivdscalar_class, only: dirichlet,neumann,quick,bquick,upwind
         real(WP) :: diffusivity
         ! Create scalar solver
         msc=multivdscalar(cfg=cfg,scheme=bquick,nscalar=speciesNum,name='Variable density multi scalar')
         call param_read('Species names',msc%SCname)

         call msc%add_bcond(name='inflow',type=dirichlet,locator=left_of_domainsc,dir='-x')
         ! Outflow on the right
         call msc%add_bcond(name='outflow',type=neumann,locator=right_of_domain,dir='+x')

         ! Assign constant diffusivity
         call param_read('Dynamic diffusivity',diffusivity)
         msc%diff=diffusivity

         ! Configure implicit scalar solver
         mss=ddadi(cfg=cfg,name='multiScalar',nst=13)
         ! Setup the solver
         call msc%setup(implicit_solver=mss)

      end block create_multiscalar

      ! Create a temperature solver
      create_temperature_solver: block
         use vdscalar_class, only: dirichlet,neumann,quick,upwind
         ! Create scalar solver
         tp=vdscalar(cfg=cfg,scheme=upwind,name='Temperature')
         ! Define boundary conditions
         call tp%add_bcond(name='inflow',type=dirichlet,locator=left_of_domainsc   )
         ! Outflow on the right
         call tp%add_bcond(name='outflow',type=neumann,locator=right_of_domain,dir='+x')
         ! Configure implicit scalar solver
         tps=ddadi(cfg=cfg,name='Temperature',nst=13)
         ! Setup the solver
         call tp%setup(implicit_solver=tps)
      end block create_temperature_solver

      ! Allocate work arrays
      allocate_work_arrays: block
         !< FLow solver
         allocate(resRHO  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcUlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcVlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcWlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcConti(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(rho0    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)) ; rho0=0.0_WP

         !< Scalar solver
         allocate(resMSC(msc%cfg%imino_:msc%cfg%imaxo_,msc%cfg%jmino_:msc%cfg%jmaxo_,msc%cfg%kmino_:msc%cfg%kmaxo_,msc%nscalar))
         allocate(srcMSC(msc%cfg%imino_:msc%cfg%imaxo_,msc%cfg%jmino_:msc%cfg%jmaxo_,msc%cfg%kmino_:msc%cfg%kmaxo_,msc%nscalar))
         allocate(mscTmp(msc%cfg%imino_:msc%cfg%imaxo_,msc%cfg%jmino_:msc%cfg%jmaxo_,msc%cfg%kmino_:msc%cfg%kmaxo_,msc%nscalar))
         allocate(bqflag(msc%cfg%imino_:msc%cfg%imaxo_,msc%cfg%jmino_:msc%cfg%jmaxo_,msc%cfg%kmino_:msc%cfg%kmaxo_,msc%nscalar))
         allocate(massFracSum(msc%cfg%imino_:msc%cfg%imaxo_,msc%cfg%jmino_:msc%cfg%jmaxo_,msc%cfg%kmino_:msc%cfg%kmaxo_))
         allocate(SCmassOutSum(msc%nscalar))

         !< Temperature solver
         allocate(resTp   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcTp   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(gTfield (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gCpfield(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gLambdafield(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

      end block allocate_work_arrays

      srcConti=0.0_WP
      srcMSC=0.0_WP
      srcTp=0.0_WP


      ! Initialize our LPT solver
      initialize_lpt: block
         use random, only: random_uniform
         use mathtools, only: Pi
         use reactionSolver_class, only: reactionSolver
         real(WP) :: dp,Hbed,VFavg,Tp,Lpart
         integer :: i,ix,iy,iz,np
         real(WP) :: particleHeatCap
         real(WP),dimension(:),allocatable :: initalComp
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the inpit
         call param_read('Drag model',lp%drag_model,default='Tenneti')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter',dp)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize

         ! Newly added=============
         ! Allocate the properties array
         allocate(lp%densities(speciesNum))
         allocate(lp%heatCapacities(speciesNum))
         allocate(lp%nongasMask(speciesNum))
         allocate(initalComp(speciesNum))

         call param_read("Particle heat Capacity", particleHeatCap, default=2.3e3_WP)
         call param_read("Particle inital composition",initalComp)
         call param_read('Particle numbers', np, default=1)
         call param_read('Species nongas densities',lp%densities)
         call param_read('Species nongas heatCapacities',lp%heatCapacities)
         call param_read('Species nongas mask',lp%nongasMask)
         call param_read('Particle temperature',Tp,default=298.15_WP)

         lp%compNum = speciesNum

         ! Root process initializes particles uniformly
         if (lp%cfg%amRoot) then
            Lpart = dp*1.5_WP
            call lp%resize(np)
            ! Distribute particles
            do i=1,np
               ! Give position
               ix = (i-1)
               lp%p(i)%pos(1) = lp%cfg%x(lp%cfg%imin)+(real(ix,WP)+0.5_WP)*Lpart
               lp%p(i)%pos(2) = 0.0_WP
               lp%p(i)%pos(3) = 0.0_WP
               ! Give id
               lp%p(i)%id=int(i,8)
               ! Set the diameter
               lp%p(i)%d=dp
               ! Give zero velocity
               lp%p(i)%vel=0.0_WP
               ! Give zero collision force
               lp%p(i)%Acol=0.0_WP
               lp%p(i)%Tcol=0.0_WP
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Locate the particle on the mesh
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])

               ! Set the temperature
               lp%p(i)%T=Tp
               ! Set initial mass
               lp%p(i)%initMass = Pi/6*dp**3.0_WP*lp%rho
               lp%p(i)%mass = lp%p(i)%initMass
               ! Set initial particle compositions
               lp%p(i)%composition = initalComp
               ! initalize other changable properties
               lp%p(i)%avgrho &
                  = dot_product(lp%densities,initalComp)
               lp%p(i)%avgCp &
                  = dot_product(lp%heatCapacities,initalComp)

               ! Activate the particle
               lp%p(i)%flag=0
            end do
         end if

         call lp%sync()
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Set collision timescale
         call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
         ! Set coefficient of restitution
         call param_read('Coefficient of restitution',lp%e_n)
         call param_read('Wall restitution',lp%e_w,default=lp%e_n)
         call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
         ! Set gravity
         call param_read('Gravity',lp%gravity)
         if (lp%cfg%amRoot) then
            print*,"===== Particle Setup Description ====="
            print*,'Number of particles', np
         end if

      end block initialize_lpt

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i,j,nvar
         nvar = 3+speciesNum
         pmesh=partmesh(nvar=nvar,nvec=3,name='lpt')
         pmesh%varname(1)='diameter'
         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='ang_vel'
         pmesh%vecname(3)='Fcol'
         pmesh%varname(2)='Temperature'
         pmesh%varname(3)='mass'
         do j=1,speciesNum
            pmesh%varname(3+j)="frac_"//trim(SCname(j))
         end do
         call lp%update_partmesh(pmesh)
         do i=1,lp%np_
            pmesh%var(1,i)=lp%p(i)%d
            pmesh%vec(:,1,i)=lp%p(i)%vel
            pmesh%vec(:,2,i)=lp%p(i)%angVel
            pmesh%vec(:,3,i)=lp%p(i)%Acol
            pmesh%var(2,i)=lp%p(i)%T
            pmesh%var(3,i)=lp%p(i)%mass
            do j=1,speciesNum
               pmesh%var(3+j,i)=lp%p(i)%composition(j)
            end do
         end do
      end block create_pmesh

      ! Initialize our multi mixture fraction field
      initialize_multiscalar: block
         use multivdscalar_class, only: bcond
         integer :: n,i,j,k,isc
         type(bcond), pointer :: mybc
         real(WP), dimension(speciesNum) :: initComp
         call param_read('Gas phase inital composition',initComp)
         ! set initial field
         do k=msc%cfg%kmino_, msc%cfg%kmaxo_
            do j=msc%cfg%jmino_, msc%cfg%jmaxo_
               do i=msc%cfg%imino_, msc%cfg%imaxo_
                  msc%SC(i,j,k,:)=initComp
               end do
            end do
         end do
         ! Apply BCs
         call msc%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            msc%SC(i,j,k,:)=inletComp
         end do

         ! Compute density
         call get_rho()
         call msc%get_max()
         call msc%get_int()

         call msc%cfg%sync(msc%rho)
         call msc%cfg%sync(fs%visc)
         do isc=1,msc%nscalar
            call msc%cfg%sync(msc%SC  (:,:,:,isc))
            call msc%cfg%sync(msc%diff(:,:,:,isc))
         end do

         call checkMassFractionSum()

      end block initialize_multiscalar

      ! Initialize our velocity field
      initialize_velocity: block
         use lowmach_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k

         ! Set inflow velocity/momentum
         call param_read('Inlet velocity',inlet_velocity)

         ! Zero initial field
         fs%U=inlet_velocity; fs%V=0.0_WP; fs%W=0.0_WP

         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%rhoU(i,j,k)=inlet_velocity*dot_product(gRhos,inletComp)
            fs%U(i,j,k)   =inlet_velocity
         end do
         fs%rho=msc%rho

         ! Form momentum
         call fs%rho_multiply()
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         call fs%interp_vel(Ui,Vi,Wi)
         ! Here it should be changed to use mvdscalar
         resRHO=0.0_WP
         call fs%get_div(drhodt=resRHO)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
      end block initialize_velocity

      ! Initialize our temperature field
      init_temp : block
         use vdscalar_class, only: bcond
         integer :: n,i,j,k
         type(bcond), pointer :: mybc

         call param_read('inital gas phase temperature',initTemp)
         do i=tp%cfg%imin_,tp%cfg%imax_
            do j=tp%cfg%jmin_,tp%cfg%jmax_
               do k=tp%cfg%kmin_,tp%cfg%kmax_
                  tp%SC(i,j,k)=initTemp*dot_product(msc%SC(i,j,k,:),gCp)
               end do
            end do
         end do

         ! Apply BCs
         call tp%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            tp%SC(i,j,k)=(2*inletTemp-tp%SC(i+1,j,k)/dot_product(msc%SC(i+1,j,k,:),gCp))*dot_product(msc%SC(i,j,k,:),gCp)
         end do

         call tp%apply_bcond(time%t,time%dt)

         tp%rho=msc%rho

         call get_tp_diff()

         ! Compute Temperature,Cp,Lambda fields for output
         do i=tp%cfg%imino_,tp%cfg%imaxo_
            do j=tp%cfg%jmino_,tp%cfg%jmaxo_
               do k=tp%cfg%kmino_,tp%cfg%kmaxo_
                  gTfield(i,j,k)=initTemp
                  gCpfield(i,j,k)=dot_product(msc%SC(i,j,k,:),gCp)
                  gLambdafield(i,j,k)=dot_product(msc%SC(i,j,k,:),gLambda)
               end do
            end do
         end do

      end block init_temp

      ! Add Ensight output
      create_ensight: block
         integer :: i
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='pyrolysisFBed')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('particleVF',lp%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('density',rho0)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('massFractionSum',massFracSum)
         ! Add MSC Comp to output
         do i=1,speciesNum
            call ens_out%add_scalar("frac_"//trim(SCname(i)),msc%SC(:,:,:,i))
         end do
         ! Add Temperature to output
         call ens_out%add_scalar('Temperature',gTfield)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      ! Create monitor filea
      create_monitor: block
         integer :: isc
         ! Prepare some info about fields
         real(WP) :: cfl
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call fs%get_max()
         call lp%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%add_column(lp%CFL_col,'Collision CFL')
         call cflfile%write()
         ! Create LPT monitor
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp%VFmean,'VFp mean')
         call lptfile%add_column(lp%VFmax,'VFp max')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%Umin,'Particle Umin')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmin,'Particle Vmin')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmin,'Particle Wmin')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%add_column(lp%dmin,'Particle dmin')
         call lptfile%add_column(lp%dmax,'Particle dmax')
         call lptfile%write()
         ! Create timing monitor
         tfile=monitor(amroot=fs%cfg%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep number')
         call tfile%add_column(time%t,'Time')
         call tfile%add_column(wt_total%time,'Total [s]')
         call tfile%add_column(wt_vel%time,'Velocity [s]')
         call tfile%add_column(wt_vel%percent,'Velocity [%]')
         call tfile%add_column(wt_pres%time,'Pressure [s]')
         call tfile%add_column(wt_pres%percent,'Pressure [%]')
         call tfile%add_column(wt_lpt%time,'LPT [s]')
         call tfile%add_column(wt_lpt%percent,'LPT [%]')
         call tfile%add_column(wt_rest%time,'Rest [s]')
         call tfile%add_column(wt_rest%percent,'Rest [%]')
         call tfile%write()
         ! Create MSC monitor
         mscfile=monitor(fs%cfg%amRoot,'msc')
         call mscfile%add_column(time%n,'Timestep number')
         call mscfile%add_column(time%t,'Time')
         call mscfile%add_column(msc%rhomax,'RHOmax')
         do isc=1,msc%nscalar
            call mscfile%add_column(msc%SCmax(isc),msc%SCname(isc)//'_max')
            call mscfile%add_column(msc%SCmin(isc),msc%SCname(isc)//'_min')
         end do
         call mscfile%write()
         ! Create out flux accumulate file
         outfluxfile=monitor(fs%cfg%amRoot,'outflux')
         call outfluxfile%add_column(time%n,'Timestep number')
         call outfluxfile%add_column(time%t,'Time')
         do isc=1,msc%nscalar
            call outfluxfile%add_column(SCmassOutSum(isc),msc%SCname(isc)//"[kg]")
         end do
         call outfluxfile%write()
      end block create_monitor

   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      use mathtools, only: twoPi
      use parallel, only: parallel_time
      use multivdscalar_class, only: bcond
      implicit none
      real(WP) :: cfl
      integer :: isc,i,j,k,n
      type(bcond), pointer :: mybc

      ! Perform time integration
      do while (.not.time%done())

         ! Initial wallclock time
         wt_total%time_in=parallel_time()

         ! Increment time
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old density, velocity, and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         ! Remember old density
         msc%rhoold=msc%rho
         msc%SCold =msc%SC

         ! Remember old density and TCp
         tp%rhoold = tp%rho
         tp%SCold = tp%SC

         ! Get fluid stress
         wt_lpt%time_in=parallel_time()
         call fs%get_div_stress(resU,resV,resW)

         ! ================================= Particle solver ==========================
         call lp%collide(dt=time%dtmid)
         call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=msc%rho,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,&
             srcU=srcUlp,srcV=srcVlp,srcW=srcWlp)

         !> intial source term array for interphase transport
         srcMSC = 0.0_WP
         srcConti = 0.0_WP
         srcTp = 0.0_WP

         !< Update temperature
         !< The source term from interphase heat transfer is in unit of [W/m^3]
         call lp%updateTemp(U=fs%U,V=fs%V,W=fs%W,visc=fs%visc,rho=rho0,gCps=gCpfield,glambdas=gLambdafield,gTemp=gTfield,dt=time%dt,srcTp=srcTp)

         !< Reaction
         !< The MSC and Mass source term are in unit of [kg/m^3/s]
         !< The heat source term is for phase change, in unit of [W/m^3]
         call lp%react(time=time%t, dt=time%dt, scSRC=srcMSC, massSRC=srcConti,heatSRC=srcTp, gCp=gCp)

         wt_lpt%time=wt_lpt%time+parallel_time()-wt_lpt%time_in



         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            ! ============= MULTI SCALAR SOLVER =======================
            call applyLewisNumber()

            !< Reset metrics for bquick
            call msc%metric_reset()

            ! Build mid-time scalar
            msc%SC=0.5_WP*(msc%SC+msc%SCold)

            ! Get rhoSC time derivative
            call msc%get_drhoSCdt(resMSC,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Assemble explicit residual
            do isc=1,msc%nscalar
               resMSC(:,:,:,isc)=time%dt*resMSC(:,:,:,isc)&
               &   -2.0_WP*msc%rho*msc%SC(:,:,:,isc)+(msc%rho+msc%rhoold)*msc%SCold(:,:,:,isc)&
               &   +time%dt*srcMSC(:,:,:,isc)  !< Source term
               !< Get temperary solution for bquick
               mscTmp(:,:,:,isc)=2.0_WP*msc%SC(:,:,:,isc)-msc%SCold(:,:,:,isc)+resMSC(:,:,:,isc)/msc%rho
            end do

            !< For each species, if <0 or >1, set flag
            bqflag = .false.
            do isc=1,msc%nscalar
               do k=msc%cfg%kmino_,msc%cfg%kmaxo_
                  do j=msc%cfg%jmino_,msc%cfg%jmaxo_
                     do i=msc%cfg%imino_,msc%cfg%imaxo_
                        if ((mscTmp(i,j,k,isc).le.0.0_WP).or.(mscTmp(i,j,k,isc).ge.1.0_WP)) then
                           bqflag(i,j,k,isc)=.true.
                        end if
                     end do
                  end do
               end do
            end do

            ! Adjust metrics
            call msc%metric_adjust(mscTmp,bqflag)

            ! re-Assemble explicit residual
            call msc%get_drhoSCdt(resMSC,fs%rhoU,fs%rhoV,fs%rhoW)
            do isc=1,msc%nscalar
               resMSC(:,:,:,isc) = time%dt*resMSC(:,:,:,isc)&
               &   -2.0_WP*msc%rho*msc%SC(:,:,:,isc)+(msc%rho+msc%rhoold)*msc%SCold(:,:,:,isc)&
               &   +time%dt*srcMSC(:,:,:,isc)  !< Source term
            end do

            call msc%solve_implicit(time%dt,resMSC,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Apply these residuals, reconstruct the n+1 field
            msc%SC=2.0_WP*msc%SC-msc%SCold+resMSC

            ! Apply boundary conditions on the resulting field
            call msc%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               msc%SC(i,j,k,:)=inletComp
            end do
            call msc%apply_bcond(time%t,time%dt)

            call get_rho()
            call msc%get_max()
            call msc%get_int()

            ! ===================== Velocity Solver =======================

            wt_vel%time_in=parallel_time()

            ! Update density based on particle volume fraction and multi-vd
            fs%rho=0.5_WP*(msc%rho+msc%rhoold)

            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)

            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)

            ! Add momentum source term from lpt
            add_lpt_src: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcUlp(i-1:i,j,k))
                        resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcVlp(i,j-1:j,k))
                        resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcWlp(i,j,k-1:k))
                     end do
                  end do
               end do
            end block add_lpt_src

            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply other boundary conditions and update momentum
            call fs%apply_bcond(time%tmid,time%dtmid)
            call fs%rho_multiply()
            call fs%apply_bcond(time%tmid,time%dtmid)

            ! Reset Dirichlet BCs
            dirichlet_velocity: block
               use lowmach_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               call fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%rhoU(i,j,k)=inlet_velocity*dot_product(gRhos,inletComp)
                  fs%U(i,j,k)   =inlet_velocity
               end do
            end block dirichlet_velocity

            wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

            ! Solve Poisson equation
            wt_pres%time_in=parallel_time()

            call msc%get_drhodt(dt=time%dt,drhodt=resRHO)
            resRHO = resRHO - srcConti
            call fs%correct_mfr(drhodt=resRHO)
            call fs%get_div(drhodt=resRHO)


            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide
            wt_pres%time=wt_pres%time+parallel_time()-wt_pres%time_in

            ! ===================== Temperature Solver =====================
            ! Get Temperature rho and diffusivity
            tp%rho=msc%rho
            call get_tp_diff()

            ! Build mid-time temperature*Cp
            tp%SC = 0.5_WP*(tp%SC+tp%SCold)

            ! Assemble explicit residual
            call tp%get_drhoSCdt(resTp,fs%rhoU,fs%rhoV,fs%rhoW)
            resTp = time%dt*resTp-(2.0_WP*tp%rho*tp%SC-(tp%rho+tp%rhoold)*tp%SCold) &
            &        +time%dt*srcTp     !< Source term, srcTP should be in unit of [W/m^3]

            ! solve implicit
            call tp%solve_implicit(time%dt,resTp,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Advance temperature field
            tp%SC=2.0_WP*tp%SC-tp%SCold+resTp

            ! Apply boundary conditions
            call tp%apply_bcond(time%t, time%dt)
            dirichlet_temperature : block
               use vdscalar_class, only : bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               call tp%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  tp%SC(i,j,k)=(2*inletTemp-tp%SC(i+1,j,k)/dot_product(msc%SC(i+1,j,k,:),gCp))*dot_product(msc%SC(i,j,k,:),gCp)
               end do
            end block dirichlet_temperature

            ! Get temperature, averaged Cp and Lambda fields for output
            do i=tp%cfg%imino_,tp%cfg%imaxo_
               do j=tp%cfg%jmino_,tp%cfg%jmaxo_
                  do k=tp%cfg%kmino_,tp%cfg%kmaxo_
                     gCpfield(i,j,k)=dot_product(msc%SC(i,j,k,:),gCp)
                     gTfield(i,j,k)=tp%SC(i,j,k)/gCpfield(i,j,k)
                     gLambdafield(i,j,k)=dot_product(msc%SC(i,j,k,:),gLambda)
                  end do
               end do
            end do

            ! =======================================================


            ! Increment sub-iteration counter

            time%it=time%it+1

         end do

         ! ==================== Post Process ==================================



         !< check mass fraction sum
         call checkMassFractionSum()

         if (fs%cfg%imaxo_ .eq. fs%cfg%imaxo) then
            call getFluxOutSum(time%t,time%dt)
         end if

         ! Recompute interpolated velocity and divergence
         wt_vel%time_in=parallel_time()
         call fs%interp_vel(Ui,Vi,Wi)
         call msc%get_drhodt(dt=time%dt,drhodt=resRHO)
         resRHO = resRHO - srcConti
         call fs%get_div(drhodt=resRHO)
         wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
               integer :: i,j
               call lp%update_partmesh(pmesh)
               do i=1,lp%np_
                  pmesh%var(1,i)=lp%p(i)%d
                  pmesh%vec(:,1,i)=lp%p(i)%vel
                  pmesh%vec(:,2,i)=lp%p(i)%angVel
                  pmesh%vec(:,3,i)=lp%p(i)%Acol
                  pmesh%var(2,i)=lp%p(i)%T
                  pmesh%var(3,i)=lp%p(i)%mass
                  do j=1,speciesNum
                     pmesh%var(3+j,i)=lp%p(i)%composition(j)
                  end do
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call lp%get_max()
         call mfile%write()
         call cflfile%write()
         call lptfile%write()
         call mscfile%write()
         call outfluxfile%write()

         ! Monitor timing
         wt_total%time=parallel_time()-wt_total%time_in
         wt_vel%percent=wt_vel%time/wt_total%time*100.0_WP
         wt_pres%percent=wt_pres%time/wt_total%time*100.0_WP
         wt_lpt%percent=wt_lpt%time/wt_total%time*100.0_WP
         wt_rest%time=wt_total%time-wt_vel%time-wt_pres%time-wt_lpt%time
         wt_rest%percent=wt_rest%time/wt_total%time*100.0_WP
         call tfile%write()
         wt_total%time=0.0_WP; wt_total%percent=0.0_WP
         wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
         wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
         wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
         wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP
      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker

      ! Deallocate work arrays
      deallocate(resU,resV,resW,srcUlp,srcVlp,srcWlp,Ui,Vi,Wi,resRHO,resMSC,resTp,mscTmp,bqflag,srcMSC,srcConti,srcTp)

   end subroutine simulation_final

end module simulation
