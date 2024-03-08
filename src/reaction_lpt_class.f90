module reaction_lpt_class
   use dev_lpt_class,     only: lpt
   use config_class,   only: config
   use monitor_class,     only: monitor
   use string, only: str_medium
   use precision, only: WP
   use, intrinsic :: iso_c_binding
   use fsundials_context_mod
   use fcvode_mod                 ! Fortran interface to CVODE
   use fnvector_serial_mod        ! Fortran interface to serial N_Vector
   use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
   use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
   use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
   use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
   use fsundials_nvector_mod      ! Fortran interface to generic N_Vector
   use param, only: param_read

   implicit none

   PUBLIC :: reaction_lpt


   real(c_double)                 :: rtol=1.0d-5, atol=1.0d-10   ! global relative and absolute tolerance
   integer(c_int)                 :: ierr                        ! error flag from C functions

   !< Each prs element should match each particle in lpt class
   !< It will only be destroyed or created when particle is destroyed or created
   !< It is a linked list so that when delete or add new noodes, the address of composition of existing nodes will not change
   type :: prs
      integer(KIND=8) :: id
      real(WP), dimension(:), allocatable :: composition
      real(WP) :: t_cur
      integer, dimension(3) :: ind
      real(WP) :: Temp
      type(prs), pointer :: next => NULL()
      !< Below are CVode related variables
      type(c_ptr)                    :: ctx          ! SUNDIALS context
      type(c_ptr)                    :: cvode_mem    ! CVODE memory
      type(N_Vector),        pointer :: sunvec_y     ! sundials vector
      type(SUNMatrix),       pointer :: sunmat_A     ! sundials matrix
      type(SUNLinearSolver), pointer :: sunlinsol_LS ! sundials linear Solver
   end type prs  !< particle's reaction solving related


   type :: reaction_lpt
      class(config), POINTER :: cfg
      type(lpt), POINTER :: lp
      character(len=str_medium) :: name='UNNAMED_LPT_Reaction'     !< Solver name (default=UNNAMED_LPT)
      integer :: np_        !< local number of particles, should be synced with lpt class
      integer(c_long) :: compNum        !< number of compositions, it should be the same as equation numbers for ODE solver
      integer(kind=8), DIMENSION(:), ALLOCATABLE :: activeParticles       !< Remenber the active particle id
      type(c_funptr) :: RhsFunc, JacFunc

      type(prs), pointer :: head => NULL() !< head of the linked particle-list

   contains
      !< Below are public operations, can be called from outside
      procedure :: react
      procedure :: syncLpt
      procedure :: proceedReact
      procedure :: writeBackLpt

      !< Below are private operations, commonly should not be called from outside
      procedure :: addParticleNode
      procedure :: deleteParticleNode
      procedure :: getActiveParticles
      procedure :: initNodeCVode
   end type reaction_lpt


   interface reaction_lpt
      procedure constructor
   end interface reaction_lpt

contains

   function constructor(cfg, lp, RhsFunc, JacFunc) result(self)
      implicit none
      type(reaction_lpt) :: self
      type(lpt), target, INTENT(IN) :: lp
      type(c_funptr), INTENT(IN) :: RhsFunc, JacFunc
      class(config),target, INTENT(IN) :: cfg
      integer :: i
      self%cfg => cfg
      self%lp => lp
      self%RhsFunc = RhsFunc
      self%JacFunc = JacFunc
      self%np_ = lp%np_
      self%compNum = lp%compNum
      self%head => NULL()
      do i=1, lp%np_
         if (lp%p(i)%flag == 0) then !< only add active particles
            call self%addParticleNode(lp%p(i)%id, lp%p(i)%composition, 0.0_WP, lp%p(i)%ind, lp%p(i)%T)
            call self%initNodeCVode(lp%p(i)%id)
         end if
      end do
      call self%getActiveParticles()
   end function constructor


   subroutine react(self, time, dt)
      implicit none
      class(reaction_lpt), INTENT(INOUT) :: self
      real(WP), target, INTENT(IN) :: time, dt
      write(*,*) '[Reaction Solver] Start to solve ODE for particles in lpt...'
      call self%syncLpt(time)
      call self%proceedReact(time, dt)
      call self%writeBackLpt()
   end subroutine react


   subroutine syncLpt(self, time)
      implicit none
      class(reaction_lpt), INTENT(INOUT) :: self
      real(WP), INTENT(IN) :: time
      integer(kind=8), dimension(:), allocatable :: activelptParticles, newParticles, deadParticles
      integer :: i,j, newCount, deadCount, activeInLpt
      real(WP), dimension(:), allocatable :: tempComp
      integer, dimension(3) :: index
      real(WP) :: Temp
      type(prs), pointer :: current
      !< Get active lpt particles, put their id into activelptParticles
      activeInLpt = 0
      allocate(activelptParticles(self%lp%np_))
      activelptParticles = -1
      do i=1, self%lp%np_
         if (self%lp%p(i)%flag==0 .and. self%lp%p(i)%id>0) then
            activelptParticles(i) = self%lp%p(i)%id
            activeInLpt = activeInLpt + 1
         end if
      end do

      !< Get active reaction particles
      call self%getActiveParticles()

      !< Compare and get new and dead particles in reaction solver
      allocate(newParticles(activeInLpt))
      newCount = 0
      if (self%np_ == 0) then
         newParticles = activelptParticles
         newCount = activeInLpt
      else
         do i=1, activeInLpt
            if (all(self%activeParticles /= activelptParticles(i))) then
               newCount = newCount + 1
               newParticles(newCount) = activelptParticles(i)
            end if
         end do
      end if

      allocate(deadParticles(self%np_))
      deadCount = 0
      if (activeInLpt == 0) then
         deadParticles = self%activeParticles
         deadCount = self%np_
      else
         do i=1, self%np_
            if (all(activelptParticles /= self%activeParticles(i))) then
               deadCount = deadCount + 1
               deadParticles(deadCount) = self%activeParticles(i)
            end if
         end do
      end if
      

      !< Add new particles and delete dead particles
      allocate(tempComp(self%compNum))
      do i=1, newCount
         !< Find the composition of new particle with such id
         do j=1, self%lp%np_
            if (newParticles(i) == self%lp%p(j)%id) then
               tempComp = self%lp%p(j)%composition
               index = self%lp%p(j)%ind
               Temp = self%lp%p(j)%T
            end if
         end do
         !< Add new particle node with id and composition, inital time would be current time
         call self%addParticleNode(newParticles(i), tempComp, time, index, Temp)
         call self%initNodeCVode(newParticles(i))
      end do
      do i=1, deadCount
         call self%deleteParticleNode(deadParticles(i))
      end do

      !< Update active particles
      call self%getActiveParticles()

      !< Update particle location index
      current => self%head
      do while (associated(current))
         do i=1, self%lp%np_
            if (current%id == self%lp%p(i)%id) then
               current%ind = self%lp%p(i)%ind
               current%Temp = self%lp%p(i)%T
            end if
         end do
         current => current%next
      end do

   end subroutine syncLpt


   subroutine proceedReact(self, t, dt)
      implicit none
      class(reaction_lpt), INTENT(INOUT) :: self
      real(WP), target, INTENT(IN) :: dt, t
      real(WP), target :: T_local
      type(prs), pointer :: current
      real(WP)                 :: tcur(1)      ! returned current time
      current => self%head
      do while (associated(current))
         !< Set user data Temperature
         T_local = current%Temp
         ierr = FCVodeSetUserData(current%cvode_mem, c_loc(T_local))
         if (ierr/=0) print *, "[CVODE] Setting user data error"

         !< Solve ODE
         ierr = FCVode(current%cvode_mem, (t+dt), current%sunvec_y, tcur, CV_NORMAL)
         if (ierr/=0) print *, "[CVODE] CVode error"

         current => current%next
      end do
   end subroutine proceedReact


   subroutine writeBackLpt(self)
      implicit none
      class(reaction_lpt), INTENT(INOUT) :: self
      type(prs), pointer :: current
      integer :: i
      logical :: found
      current => self%head
      do while (associated(current))
         found = .false.
         do i=1, self%lp%np_
            !< write beck only when the particle id is the same and the particle is active
            if ((current%id == self%lp%p(i)%id) .and. (self%lp%p(i)%flag==0)) then
               self%lp%p(i)%composition = current%composition
               found = .true.
            end if
         end do
         if (.not. found) print *, '[Reaction][Error] No such particle id in lpt, writeBackLpt failed'
         current => current%next
      end do
   end subroutine writeBackLpt


!< Belowing are private operations

   subroutine addParticleNode(self, id, comp, t_cur, index, Temp)
      implicit none
      class(reaction_lpt), INTENT(INOUT) :: self
      integer(KIND=8), INTENT(IN) :: id
      real(WP), dimension(:), INTENT(IN) :: comp
      real(WP), INTENT(IN) :: t_cur, Temp
      integer, dimension(3),  INTENT(IN) :: index
      type(prs), pointer :: new_node
      type(prs), pointer :: current

      ! write(*,*) '[Reaction][Info] Add new particle node with id:', id
      allocate(new_node)
      new_node%id = id
      new_node%composition = comp
      new_node%t_cur = t_cur
      new_node%next => NULL()
      new_node%ind = index
      new_node%Temp = Temp

      if (associated(self%head)) then
         current => self%head
         do while (associated(current%next))
            current => current%next
         end do
         current%next => new_node
      else
         self%head => new_node
      end if
   end subroutine addParticleNode

   subroutine deleteParticleNode(self, id)
      implicit none
      class(reaction_lpt), INTENT(INOUT) :: self
      integer(KIND=8), INTENT(IN) :: id
      type(prs), pointer :: current
      type(prs), pointer :: previous
      logical :: found

      found = .false.
      if (associated(self%head)) then
         current => self%head
         previous => NULL()
         do while (associated(current))
            if (current%id == id) then
               found = .true.
               if (associated(previous)) then
                  previous%next => current%next
               else
                  self%head => current%next
               end if
               call FCVodeFree(current%cvode_mem)
               ierr = FSUNLinSolFree(current%sunlinsol_LS)
               call FSUNMatDestroy(current%sunmat_A)
               call FN_VDestroy(current%sunvec_y)
               ierr = FSUNContext_Free(current%ctx)
               deallocate(current%composition)
               deallocate(current)
               return
            end if
            previous => current
            current => current%next
         end do
         if (.not. found) print *, '[Reaction][Error] No such particle id to delete'
      else
         print *, '[Reaction][Error] No node to delete'
      end if
   end subroutine deleteParticleNode

   subroutine getActiveParticles(self)
      implicit none
      class(reaction_lpt), INTENT(INOUT) :: self
      integer :: i
      integer(kind = 8), dimension(:), allocatable :: activeParticles
      integer :: count
      type(prs), pointer :: current

      current => self%head
      count = 0
      do while (associated(current))
         count = count + 1
         current => current%next
      end do

      allocate(activeParticles(count))
      current => self%head
      do i = 1, count
         activeParticles(i) = current%id
         current => current%next
      end do
      call move_alloc(activeParticles, self%activeParticles)
      self%np_ = count
   end subroutine getActiveParticles

   subroutine initNodeCVode(self, id)
      implicit none
      class(reaction_lpt), INTENT(INOUT) :: self
      integer(kind=8), INTENT(IN) :: id
      type(prs), pointer :: current
      logical :: found

      found = .false.
      current => self%head
      do while (associated(current))
         if (current%id == id) then
            found = .true.
            ierr = FSUNContext_Create(c_null_ptr, current%ctx)
            if (ierr/=0) print *, "[CVODE] context create error"
            current%sunvec_y => FN_VMake_Serial(self%compNum, current%composition, current%ctx)
            current%sunmat_A => FSUNDenseMatrix(self%compNum, self%compNum, current%ctx)
            current%sunlinsol_LS => FSUNLinSol_Dense(current%sunvec_y, current%sunmat_A, current%ctx)
            if (.not. associated(current%sunvec_y)) print *, "[CVODE] sunvec_y create failed"
            if (.not. associated(current%sunmat_A)) print *, "[CVODE] sunmat_A create failed"
            if (.not. associated(current%sunlinsol_LS)) print *, "[CVODE] sunlinsol_LS create failed"
            current%cvode_mem = FCVodeCreate(CV_BDF, current%ctx)
            if (.not. c_associated(current%cvode_mem)) print *, "[CVODE] CVODE_memory create failed"
            ierr = FCVodeInit(current%cvode_mem, self%RhsFunc, current%t_cur, current%sunvec_y)
            if (ierr/=0) print *, "[CVODE] CVode Init error"
            ierr = FCVodeSStolerances(current%cvode_mem, rtol, atol)
            if (ierr/=0) print *, "[CVODE] CVode Set tolerance error"
            ierr = FCVodeSetLinearSolver(current%cvode_mem, current%sunlinsol_LS, current%sunmat_A)
            if (ierr/=0) print *, "[CVODE] CVode Set Linear Solver error"
            ierr = FCVodeSetJacFn(current%cvode_mem, self%JacFunc)
            if (ierr/=0) print *, "[CVODE] CVode Set Jacobian Func error"
            return
         end if
         current => current%next
      end do
      if (.not. found) print *, "[CVODE][Error] No such particle id, initNodeCVode failed"
   end subroutine initNodeCVode


end module reaction_lpt_class
