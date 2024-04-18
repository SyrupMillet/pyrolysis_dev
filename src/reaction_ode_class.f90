! ===========Usage Instruction==========================================

! The only thing you need to modify is in the ode_mod.

! neq is the number of euqations in the ode system.

! You also need to modify the RhsFn and JacFn
! They are based on your equations

! The reaction kinetic parameters are dynamically calculated
! So you need to supply activation energy and pre-exponential factor

! =====================================================================

module reaction_ode_mod
   use, intrinsic :: iso_c_binding

   implicit none
   ! number of equations
   integer(c_long), parameter :: neq = 5

   ! kinetic parameters
   real(c_double),dimension(5), parameter :: Ea = [206.2,295.9,206.1,325.2,63.0]
   real(c_double),dimension(5), parameter :: k0 = [1.32E+13/60,2.56E+19/60,5.40E+13/60,4.99E+18/60,1.86E+01/60]

contains
   ! ----------------------------------------------------------------
   ! RhsFn provides the right hand side function for the
   ! ODE: dy/dt = f(t,y)
   !
   ! Return values:
   !    0 = success,
   !    1 = recoverable error,
   !   -1 = non-recoverable error
   ! ----------------------------------------------------------------
   integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
      result(ierr) bind(C,name='RhsFn')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod

      !======= Declarations =========
      implicit none
      ! User data variables
      real(c_double), pointer :: T
      real(c_double) :: k1,k2,k3,k4,k5

      ! calling variables
      real(c_double), value :: tn        ! current time
      type(N_Vector)        :: sunvec_y  ! solution N_Vector
      type(N_Vector)        :: sunvec_f  ! rhs N_Vector
      type(c_ptr),    value :: user_data ! user-defined data

      ! pointers to data in SUNDIALS vectors
      real(c_double), pointer :: yvec(:)
      real(c_double), pointer :: fvec(:)

      !======= Internals ============

      ! get data arrays from SUNDIALS vectors
      yvec => FN_VGetArrayPointer(sunvec_y)
      if (.not. associated(yvec)) then
         write(*,*) 'ERROR: yvec = NULL'
      end if
      fvec => FN_VGetArrayPointer(sunvec_f)
      if (.not. associated(fvec)) then
         write(*,*) 'ERROR: fvec = NULL'
      end if

      call c_f_pointer(user_data, T)
      ! write(*,*) 'T = ', T

      ! calculate kinetic parameters based on T
      k1 = k0(1) * exp(-Ea(1)*1000.0/(8.3145*T))
      k2 = k0(2) * exp(-Ea(2)*1000.0/(8.3145*T))
      k3 = k0(3) * exp(-Ea(3)*1000.0/(8.3145*T))
      k4 = k0(4) * exp(-Ea(4)*1000.0/(8.3145*T))
      k5 = k0(5) * exp(-Ea(5)*1000.0/(8.3145*T))

      ! fill RHS vector
      ! [Nitrogen, Polymer, Gas, Liquid, Wax]
      fvec(1) = 0.0d0
      fvec(2) = -(k1+k2+k3)*yvec(2)
      fvec(3) = k1*yvec(2)+k4*yvec(5)
      fvec(4) = k2*yvec(2)+k5*yvec(5)
      fvec(5) = k3*yvec(2)-k4*yvec(5)-k5*yvec(5)

      ! return success
      ierr = 0
      return

   end function RhsFn

   ! ----------------------------------------------------------------
   ! JacFn: The Jacobian of the ODE hand side function J = df/dy
   !
   ! Return values:
   !    0 = success,
   !    1 = recoverable error,
   !   -1 = non-recoverable error
   ! ----------------------------------------------------------------
   integer(c_int) function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, &
      user_data, tmp1, tmp2, tmp3) &
      result(ierr) bind(C,name='JacFn')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod
      use fsunmatrix_dense_mod
      use fsundials_matrix_mod

      !======= Declarations =========
      implicit none
      ! User data variables
      real(c_double), pointer :: T
      real(c_double) :: k1,k2,k3,k4,k5

      ! calling variables
      real(c_double), value :: tn               ! current time
      type(N_Vector)        :: sunvec_y         ! current solution N_Vector
      type(N_Vector)        :: sunvec_f         ! current rhs N_Vector
      type(SUNMatrix)       :: sunmat_J         ! Jacobian SUNMatrix
      type(c_ptr), value    :: user_data        ! user-defined data
      type(N_Vector)        :: tmp1, tmp2, tmp3 ! workspace N_Vectors

      ! pointers to data in SUNDIALS vectors
      real(c_double), pointer :: yvec(:)

      ! pointer to data in SUNDIALS matrix
      real(c_double), pointer :: Jmat(:)

      call c_f_pointer(user_data, T)
      ! write(*,*) '[JacFunc]T = ', T
      !======= Internals ============

      ! get data array from SUNDIALS vector
      yvec => FN_VGetArrayPointer(sunvec_y)
      if (.not. associated(yvec)) then
         write(*,*) 'ERROR: yvec = NULL'
      end if
      ! get data arrays from SUNDIALS vectors
      Jmat => FSUNDenseMatrix_Data(sunmat_J)
      if (.not. associated(Jmat)) then
         write(*,*) 'ERROR: Jmat = NULL'
      end if

      ! calculate kinetic parameters based on T
      k1 = k0(1) * exp(-Ea(1)*1000.0/(8.3145*T))
      k2 = k0(2) * exp(-Ea(2)*1000.0/(8.3145*T))
      k3 = k0(3) * exp(-Ea(3)*1000.0/(8.3145*T))
      k4 = k0(4) * exp(-Ea(4)*1000.0/(8.3145*T))
      k5 = k0(5) * exp(-Ea(5)*1000.0/(8.3145*T))

      ! [Nitrogen, Polymer, Gas, Liquid, Wax]
      ! fvec(1) = 0.0d0
      ! fvec(2) = -(k1+k2+k3)*yvec(2)
      ! fvec(3) = k1*yvec(2)+k4*yvec(5)
      ! fvec(4) = k2*yvec(2)+k5*yvec(5)
      ! fvec(5) = k3*yvec(2)-k4*yvec(5)-k5*yvec(5)
      ! fill Jacobian matrix
      Jmat = [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,&
               0.0d0,-(k1+k2+k3),k1,k2,k3,&
               0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,&
               0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,&
               0.0d0,0.0d0,k4,k5,-k4-k5]
      ! return success
      ierr = 0
      return

   end function JacFn

end module reaction_ode_mod

! Don't modify below this line unless you know what you are doing
! ----------------------------------------------------------------

module reactionSolver_class
   use, intrinsic :: iso_c_binding
   use fsundials_context_mod
   use fcvode_mod                 ! Fortran interface to CVODE
   use fnvector_serial_mod        ! Fortran interface to serial N_Vector
   use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
   use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
   use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
   use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
   use fsundials_nvector_mod      ! Fortran interface to generic N_Vector
   use reaction_ode_mod                    ! ODE functions
   use param, only: param_read

   implicit none

   real(c_double)                 :: rtol=1.0d-5, atol=1.0d-10   ! relative and absolute tolerance
   integer(c_int)                 :: ierr         ! error flag from C functions

   type :: reactionSolver
      type(c_ptr)                    :: ctx          ! SUNDIALS context
      type(c_ptr)                    :: cvode_mem    ! CVODE memory
      type(N_Vector),        pointer :: sunvec_y     ! sundials vector
      type(SUNMatrix),       pointer :: sunmat_A     ! sundials matrix
      type(SUNLinearSolver), pointer :: sunlinsol_LS ! sundials linear solver

      real(c_double) :: t_init

   contains
      procedure :: proceedReaction
      procedure :: clean
   end type reactionSolver

   interface reactionSolver
      procedure constructor
   end interface reactionSolver

contains

   function constructor(composition) result(self)
      implicit none
      type(reactionSolver) :: self
      real(c_double), dimension(neq), intent(inout) :: composition
      ierr = FSUNContext_Create(c_null_ptr, self%ctx)

      self%t_init = 0.0d0

      ! create SUNDIALS N_Vector with initial values
      self%sunvec_y => FN_VMake_Serial(neq, composition, self%ctx)
      if (.not. associated(self%sunvec_y)) then
         print *, 'ERROR: sunvec = NULL'
         stop 1
      end if

      self%sunmat_A => FSUNDenseMatrix(neq, neq, self%ctx)
      if (.not. associated(self%sunmat_A)) then
         print *, 'ERROR: sunmat = NULL'
         stop 1
      end if

      ! create a linear solver
      self%sunlinsol_LS => FSUNLinSol_Dense(self%sunvec_y, self%sunmat_A, self%ctx)
      if (.not. associated(self%sunlinsol_LS)) then
         print *, 'ERROR: sunlinsol = NULL'
         stop 1
      end if

      ! create CVode memory
      self%cvode_mem = FCVodeCreate(CV_BDF, self%ctx)
      if (.not. c_associated(self%cvode_mem)) then
         print *, 'ERROR: cvode_mem = NULL'
         stop 1
      end if

      ! initialize CVode
      ierr = FCVodeInit(self%cvode_mem, c_funloc(RhsFn), self%t_init, self%sunvec_y)
      if (ierr /= 0) then
         print *, 'ERROR: FCVodeInit failed'
         stop 1
      end if
      ! set relative and absolute tolerances
      ierr = FCVodeSStolerances(self%cvode_mem, rtol, atol)
      if (ierr /= 0) then
         print *, 'ERROR: FCVodeInit failed'
         stop 1
      end if
      ! attach linear solver
      ierr = FCVodeSetLinearSolver(self%cvode_mem, self%sunlinsol_LS, self%sunmat_A)
      if (ierr /= 0) then
         print *, 'ERROR: FCVodeInit failed'
         stop 1
      end if
      ! set Jacobian routine
      ierr = FCVodeSetJacFn(self%cvode_mem, c_funloc(JacFn))
      if (ierr /= 0) then
         print *, 'ERROR: FCVodeInit failed'
         stop 1
      end if
   end function constructor

   subroutine proceedReaction(self,time, dt, T, comp)
      implicit none
      class(reactionSolver) :: self
      real(c_double) :: dt,time
      real(c_double), target :: T
      real(c_double), dimension(neq) :: comp
      real(c_double)                 :: tcur(1)      ! current time

      ! free old sunvecy and linear solver
      ierr = FSUNLinSolFree(self%sunlinsol_LS)
      call FN_VDestroy(self%sunvec_y)

      self%sunvec_y => FN_VMake_Serial(neq, comp, self%ctx)
      if (.not. associated(self%sunvec_y)) then
         print *, 'ERROR: sunvec = NULL'
         stop 1
      end if

      ! re-create linear solver
      self%sunlinsol_LS => FSUNLinSol_Dense(self%sunvec_y, self%sunmat_A, self%ctx)
      if (.not. associated(self%sunlinsol_LS)) then
         print *, 'ERROR: sunlinsol = NULL'
         stop 1
      end if

      ! re-set linear solver
      ierr = FCVodeSetLinearSolver(self%cvode_mem, self%sunlinsol_LS, self%sunmat_A)
      if (ierr /= 0) then
         print *, 'ERROR: FCVodeSetLinearSolver failed'
         stop 1
      end if

      ! set Jacobian routine
      ierr = FCVodeSetJacFn(self%cvode_mem, c_funloc(JacFn))
      if (ierr /= 0) then
         print *, 'ERROR: FCVodeSetJacFn failed'
         stop 1
      end if

      ! re-init cvode
      ierr = FCVodeReInit(self%cvode_mem, time, self%sunvec_y)
      if (ierr /= 0) then
         print *, 'ERROR: FCVodeReInit failed'
         stop 1
      end if   



      ! Assign Temperature userdata
      ierr = FCVodeSetUserData(self%cvode_mem, c_loc(T))
        if (ierr /= 0) then
         print *, 'ERROR: FCVodeSetUserData failed'
         stop 1
      end if

      ierr = FCVode(self%cvode_mem, (time+dt), self%sunvec_y, tcur, CV_NORMAL)
      if (ierr /= 0) then
         print *, 'ERROR: FCVode failed', ierr
         stop 1
      end if

   end subroutine proceedReaction


   subroutine clean(self)
      use fcvode_mod
      implicit none
      class(reactionSolver) :: self
      call FCVodeFree(self%cvode_mem)
      ierr = FSUNLinSolFree(self%sunlinsol_LS)
      call FSUNMatDestroy(self%sunmat_A)
      call FN_VDestroy(self%sunvec_y)
      ierr = FSUNContext_Free(self%ctx)
   end subroutine clean

end module reactionSolver_class
