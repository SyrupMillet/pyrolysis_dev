module react_ode_mod
   use, intrinsic :: iso_c_binding

   implicit none

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

end module react_ode_mod