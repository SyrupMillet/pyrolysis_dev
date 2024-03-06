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
    type :: prs 
        integer(KIND=8) :: id
        real(WP), dimension(:), allocatable :: composition
        real(WP) :: t_cur
        logical :: active
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
        character(len=str_medium) :: name='UNNAMED_LPT'     !< Solver name (default=UNNAMED_LPT)
        integer :: np_        !< local number of particles, should be synced with lpt class
        type(prs), dimension(:), allocatable :: p 
        integer(c_long) :: compNum        !< number of compositions, it should be the same as equation numbers for ODE solver
        integer(kind=8), DIMENSION(:), ALLOCATABLE :: activeParticles       !< Remenber the active particle id 
    contains
        procedure :: checkShare
        procedure :: proceedReact 
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
        self%compNum = lp%compNum
        self%np_ = lp%np_

        allocate(self%p(self%np_))
        allocate(self%activeParticles(self%np_))

        do i=i, self%np_
            allocate(self%p(i)%composition(self%compNum))
            self%p(i)%composition = lp%p(i)%composition
            self%p(i)%id = lp%p(i)%id
            self%p(i)%t_cur = 0.0_WP
            self%p(i)%active = .True. 
            self%activeParticles(i) = lp%p(i)%id 
            
            !< Below we began to initiate our CVode solver
            ierr = FSUNContext_Create(c_null_ptr, self%p(i)%ctx)
            if (ierr/=0) print *, "[CVODE] context create error"
            self%p(i)%sunvec_y => FN_VMake_Serial(self%compNum, self%p(i)%composition, self%p(i)%ctx)
            self%p(i)%sunmat_A => FSUNDenseMatrix(self%compNum, self%compNum, self%p(i)%ctx)
            self%p(i)%sunlinsol_LS => FSUNLinSol_Dense(self%p(i)%sunvec_y, self%p(i)%sunmat_A, self%p(i)%ctx)
            if (.not. associated(self%p(i)%sunvec_y)) print *, "[CVODE] sunvec_y create failed"
            if (.not. associated(self%p(i)%sunmat_A)) print *, "[CVODE] sunmat_A create failed"
            if (.not. associated(self%p(i)%sunlinsol_LS)) print *, "[CVODE] sunlinsol_LS create failed"
            self%p(i)%cvode_mem = FCVodeCreate(CV_BDF, self%p(i)%ctx)
            if (.not. c_associated(self%p(i)%cvode_mem)) print *, "[CVODE] CVODE_memory create failed"
            ierr = FCVodeInit(self%p(i)%cvode_mem, RhsFunc, 0.0_WP, self%p(i)%sunvec_y)
            if (ierr/=0) print *, "[CVODE] CVode Init error"
            ierr = FCVodeSStolerances(self%p(i)%cvode_mem, rtol, atol)
            if (ierr/=0) print *, "[CVODE] CVode Set tolerance error"
            ierr = FCVodeSetLinearSolver(self%p(i)%cvode_mem, self%p(i)%sunlinsol_LS, self%p(i)%sunmat_A)
            if (ierr/=0) print *, "[CVODE] CVode Set Linear Solver error"
            ierr = FCVodeSetJacFn(self%p(i)%cvode_mem, JacFunc)
            if (ierr/=0) print *, "[CVODE] CVode Set Jacobian Func error"

        end do 


    end function constructor 

    subroutine checkShare(self)
        implicit none
        class(reaction_lpt) :: self
        integer :: i, j  
        integer(kind=8), DIMENSION(:), ALLOCATABLE :: newParticles, deadParticles
        integer :: maxpnum , newNum, deadNum
        logical :: isin 
        maxpnum = self%np_ + self%lp%np_
        ALLOCATE(newParticles(maxpnum)) ; newParticles = -1
        ALLOCATE(deadParticles(maxpnum)) ; deadParticles = -1
        newNum = 0 
        deadNum = 0

        ! This part the complexibility is O(n^2), not efficient. 
        countlptP : do i=1, self%lp%np_
            if (.not. isInReact(self%activeParticles, self%np_, self%lp%p(i)%id)) then
                newParticles(newNum+1) = self%lp%p(i)%id
                newNum = newNum+1 
            end if
        end do countlptP

        countReactP : do i=1, self%np_
            isin = .false. 
            do j=1, self%lp%np_
                if (self%activeParticles(i) .eq. self%lp%p(j)%id) isin = .true. 
            end do
            if (.not. isin) then
                deadNum = deadNum+1 
                deadParticles(deadNum) = self%activeParticles(i) 
            end if 
        end do countReactP

        !< Now we have new particles needs to be added to reaction solver 
        !< And dead particles needs to be removed from reaction solver
        !< But if no new no dead, nothing needs to be done
        if ((newNum .ne. 0) .or. (deadNum .ne. 0)) then
            




        end if 

        DEALLOCATE(newParticles) ; DEALLOCATE(deadParticles)

    contains

        function isInReact(activeParticles, length, id) result(ll)
            implicit none
            integer(kind=8), dimension(:), intent(in) :: activeParticles
            integer, intent(in) :: length
            integer(kind=8) :: id
            integer :: i 
            logical :: ll 
            ll = .false.
            do i=1, length
                if (id .eq. activeParticles(i)) then
                    ll = .true. 
                end if 
            end do 
        end function isInReact

    end subroutine checkShare




    subroutine proceedReact()

    end subroutine proceedReact




end module reaction_lpt_class
