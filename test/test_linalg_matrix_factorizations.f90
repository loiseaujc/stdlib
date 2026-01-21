! Test QR factorization 
module test_linalg_matrix_factorizations
    use testdrive, only: error_type, check, new_unittest, unittest_type
    use stdlib_linalg_constants
    use stdlib_linalg_state, only: LINALG_VALUE_ERROR,linalg_state_type
    use stdlib_linalg, only: qrfact, qr, mnorm
    use stdlib_math, only: all_close
    use ieee_arithmetic, only: ieee_value,ieee_quiet_nan

    implicit none (type,external)
    
    public :: test_qr_factorization

    real(sp), parameter :: rel_tol_sp = sqrt(epsilon(1.0_sp))
    real(dp), parameter :: rel_tol_dp = sqrt(epsilon(1.0_dp))

    contains

    !> QR factorization tests
    subroutine test_qr_factorization(tests)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: tests(:)
        
        allocate(tests(0))
        
        call add_test(tests,new_unittest("qr_random_tall_matrix_s",test_qr_random_tall_matrix_s))
        call add_test(tests,new_unittest("qr_random_wide_matrix_s",test_qr_random_wide_matrix_s))
        call add_test(tests,new_unittest("qr_random_tall_matrix_d",test_qr_random_tall_matrix_d))
        call add_test(tests,new_unittest("qr_random_wide_matrix_d",test_qr_random_wide_matrix_d))
        call add_test(tests,new_unittest("qr_random_tall_matrix_c",test_qr_random_tall_matrix_c))
        call add_test(tests,new_unittest("qr_random_wide_matrix_c",test_qr_random_wide_matrix_c))
        call add_test(tests,new_unittest("qr_random_tall_matrix_z",test_qr_random_tall_matrix_z))
        call add_test(tests,new_unittest("qr_random_wide_matrix_z",test_qr_random_wide_matrix_z))
    end subroutine test_qr_factorization

    !> QR factorization of a random matrix
    subroutine test_qr_random_tall_matrix_s(error)
        use stdlib_linalg, only: hermitian, qr_rsp_type
        type(error_type), allocatable, intent(out) :: error
        integer(ilp), parameter :: m   = 15_ilp
        integer(ilp), parameter :: n   =  4_ilp
        integer(ilp), parameter :: k   = min(m,n)
        real(sp) :: a(m,n), q(m, k), r(k, n)
        real(sp) :: rea(m,n),ima(m,n)
        type(linalg_state_type) :: state
        type(qr_rsp_type) :: F
        
        call random_number(rea)
        a = rea
        
        ! Reference QR decomposition.
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a, q, r, err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Instantiate derived-type for the QR.
        F = qrfact(A)

        ! Check the R matrix with reference.
        call check(error, all_close(R, F%R(), rel_tol=rel_tol_sp))
        if (allocated(error)) return

        ! Check the Q matrix with reference.
        call check(error, all_close(Q, F%Q()))
        if (allocated(error)) return
        
    end subroutine test_qr_random_tall_matrix_s

    subroutine test_qr_random_wide_matrix_s(error)
        use stdlib_linalg, only: hermitian, qr_rsp_type
        type(error_type), allocatable, intent(out) :: error
        integer(ilp), parameter :: m   = 4_ilp
        integer(ilp), parameter :: n   = 15_ilp
        integer(ilp), parameter :: k   = min(m,n)
        real(sp) :: a(m, n), q(m, k), r(k, n)
        real(sp) :: rea(m, n),ima(m, n)
        type(linalg_state_type) :: state
        type(qr_rsp_type) :: F
        
        call random_number(rea)
        a = rea
        
        ! Reference QR decomposition.
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a, q, r, err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Instantiate derived-type for the QR.
        F = qrfact(A)

        ! Check the R matrix with reference.
        call check(error, all_close(R, F%R(), rel_tol=rel_tol_sp))
        if (allocated(error)) return

        ! Check the Q matrix with reference.
        call check(error, all_close(Q, F%Q()))
        if (allocated(error)) return
        
    end subroutine test_qr_random_wide_matrix_s
    subroutine test_qr_random_tall_matrix_d(error)
        use stdlib_linalg, only: hermitian, qr_rdp_type
        type(error_type), allocatable, intent(out) :: error
        integer(ilp), parameter :: m   = 15_ilp
        integer(ilp), parameter :: n   =  4_ilp
        integer(ilp), parameter :: k   = min(m,n)
        real(dp) :: a(m,n), q(m, k), r(k, n)
        real(dp) :: rea(m,n),ima(m,n)
        type(linalg_state_type) :: state
        type(qr_rdp_type) :: F
        
        call random_number(rea)
        a = rea
        
        ! Reference QR decomposition.
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a, q, r, err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Instantiate derived-type for the QR.
        F = qrfact(A)

        ! Check the R matrix with reference.
        call check(error, all_close(R, F%R(), rel_tol=rel_tol_dp))
        if (allocated(error)) return

        ! Check the Q matrix with reference.
        call check(error, all_close(Q, F%Q()))
        if (allocated(error)) return
        
    end subroutine test_qr_random_tall_matrix_d

    subroutine test_qr_random_wide_matrix_d(error)
        use stdlib_linalg, only: hermitian, qr_rdp_type
        type(error_type), allocatable, intent(out) :: error
        integer(ilp), parameter :: m   = 4_ilp
        integer(ilp), parameter :: n   = 15_ilp
        integer(ilp), parameter :: k   = min(m,n)
        real(dp) :: a(m, n), q(m, k), r(k, n)
        real(dp) :: rea(m, n),ima(m, n)
        type(linalg_state_type) :: state
        type(qr_rdp_type) :: F
        
        call random_number(rea)
        a = rea
        
        ! Reference QR decomposition.
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a, q, r, err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Instantiate derived-type for the QR.
        F = qrfact(A)

        ! Check the R matrix with reference.
        call check(error, all_close(R, F%R(), rel_tol=rel_tol_dp))
        if (allocated(error)) return

        ! Check the Q matrix with reference.
        call check(error, all_close(Q, F%Q()))
        if (allocated(error)) return
        
    end subroutine test_qr_random_wide_matrix_d
    subroutine test_qr_random_tall_matrix_c(error)
        use stdlib_linalg, only: hermitian, qr_csp_type
        type(error_type), allocatable, intent(out) :: error
        integer(ilp), parameter :: m   = 15_ilp
        integer(ilp), parameter :: n   =  4_ilp
        integer(ilp), parameter :: k   = min(m,n)
        complex(sp) :: a(m,n), q(m, k), r(k, n)
        real(sp) :: rea(m,n),ima(m,n)
        type(linalg_state_type) :: state
        type(qr_csp_type) :: F
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        
        ! Reference QR decomposition.
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a, q, r, err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Instantiate derived-type for the QR.
        F = qrfact(A)

        ! Check the R matrix with reference.
        call check(error, all_close(R, F%R(), rel_tol=rel_tol_sp))
        if (allocated(error)) return

        ! Check the Q matrix with reference.
        call check(error, all_close(Q, F%Q()))
        if (allocated(error)) return
        
    end subroutine test_qr_random_tall_matrix_c

    subroutine test_qr_random_wide_matrix_c(error)
        use stdlib_linalg, only: hermitian, qr_csp_type
        type(error_type), allocatable, intent(out) :: error
        integer(ilp), parameter :: m   = 4_ilp
        integer(ilp), parameter :: n   = 15_ilp
        integer(ilp), parameter :: k   = min(m,n)
        complex(sp) :: a(m, n), q(m, k), r(k, n)
        real(sp) :: rea(m, n),ima(m, n)
        type(linalg_state_type) :: state
        type(qr_csp_type) :: F
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=sp)
        
        ! Reference QR decomposition.
        q = ieee_value(0.0_sp,ieee_quiet_nan)
        r = ieee_value(0.0_sp,ieee_quiet_nan)
        call qr(a, q, r, err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Instantiate derived-type for the QR.
        F = qrfact(A)

        ! Check the R matrix with reference.
        call check(error, all_close(R, F%R(), rel_tol=rel_tol_sp))
        if (allocated(error)) return

        ! Check the Q matrix with reference.
        call check(error, all_close(Q, F%Q()))
        if (allocated(error)) return
        
    end subroutine test_qr_random_wide_matrix_c
    subroutine test_qr_random_tall_matrix_z(error)
        use stdlib_linalg, only: hermitian, qr_cdp_type
        type(error_type), allocatable, intent(out) :: error
        integer(ilp), parameter :: m   = 15_ilp
        integer(ilp), parameter :: n   =  4_ilp
        integer(ilp), parameter :: k   = min(m,n)
        complex(dp) :: a(m,n), q(m, k), r(k, n)
        real(dp) :: rea(m,n),ima(m,n)
        type(linalg_state_type) :: state
        type(qr_cdp_type) :: F
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        
        ! Reference QR decomposition.
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a, q, r, err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Instantiate derived-type for the QR.
        F = qrfact(A)

        ! Check the R matrix with reference.
        call check(error, all_close(R, F%R(), rel_tol=rel_tol_dp))
        if (allocated(error)) return

        ! Check the Q matrix with reference.
        call check(error, all_close(Q, F%Q()))
        if (allocated(error)) return
        
    end subroutine test_qr_random_tall_matrix_z

    subroutine test_qr_random_wide_matrix_z(error)
        use stdlib_linalg, only: hermitian, qr_cdp_type
        type(error_type), allocatable, intent(out) :: error
        integer(ilp), parameter :: m   = 4_ilp
        integer(ilp), parameter :: n   = 15_ilp
        integer(ilp), parameter :: k   = min(m,n)
        complex(dp) :: a(m, n), q(m, k), r(k, n)
        real(dp) :: rea(m, n),ima(m, n)
        type(linalg_state_type) :: state
        type(qr_cdp_type) :: F
        
        call random_number(rea)
        call random_number(ima)
        a = cmplx(rea,ima,kind=dp)
        
        ! Reference QR decomposition.
        q = ieee_value(0.0_dp,ieee_quiet_nan)
        r = ieee_value(0.0_dp,ieee_quiet_nan)
        call qr(a, q, r, err=state)
        call check(error,state%ok(),state%print())
        if (allocated(error)) return

        ! Instantiate derived-type for the QR.
        F = qrfact(A)

        ! Check the R matrix with reference.
        call check(error, all_close(R, F%R(), rel_tol=rel_tol_dp))
        if (allocated(error)) return

        ! Check the Q matrix with reference.
        call check(error, all_close(Q, F%Q()))
        if (allocated(error)) return
        
    end subroutine test_qr_random_wide_matrix_z

    ! gcc-15 bugfix utility
    subroutine add_test(tests,new_test)
        type(unittest_type), allocatable, intent(inout) :: tests(:)    
        type(unittest_type), intent(in) :: new_test
        
        integer :: n
        type(unittest_type), allocatable :: new_tests(:)
        
        if (allocated(tests)) then 
            n = size(tests)
        else
            n = 0
        end if
        
        allocate(new_tests(n+1))
        if (n>0) new_tests(1:n) = tests(1:n)
                 new_tests(1+n) = new_test
        call move_alloc(from=new_tests,to=tests)        
        
    end subroutine add_test

end module test_linalg_matrix_factorizations

program test_matrix_factorizations
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_linalg_matrix_factorizations, only : test_qr_factorization
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_matrix_factorizations: qr", test_qr_factorization) &
        ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if
end program test_matrix_factorizations
