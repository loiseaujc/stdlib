! Test Schur decomposition
module test_linalg_expm
    use testdrive, only: error_type, check, new_unittest, unittest_type
    use stdlib_linalg_constants
    use stdlib_linalg, only: expm, eye, norm, matrix_exp
    use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_ERROR, &
         LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR

    implicit none (type,external)
    
    public :: test_expm_computation

    contains

    !> schur decomposition tests
    subroutine test_expm_computation(tests)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: tests(:)
        
        allocate(tests(0))
        
        tests = [tests, new_unittest("expm_s",test_expm_s)]
        ! tests = [tests, new_unittest("Error-handling expm_s",test_error_handling_expm_s)]
        tests = [tests, new_unittest("expm_d",test_expm_d)]
        ! tests = [tests, new_unittest("Error-handling expm_d",test_error_handling_expm_d)]
        tests = [tests, new_unittest("expm_c",test_expm_c)]
        ! tests = [tests, new_unittest("Error-handling expm_c",test_error_handling_expm_c)]
        tests = [tests, new_unittest("expm_z",test_expm_z)]
        ! tests = [tests, new_unittest("Error-handling expm_z",test_error_handling_expm_z)]

    end subroutine test_expm_computation

    !> Matrix exponential with analytic expression.
    subroutine test_expm_s(error)
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer(ilp), parameter :: n = 5, m = 6
        ! Test matrix.
        real(sp) :: A(n, n), E(n, n), Eref(n, n)
        real(sp) :: err
        integer(ilp) :: i, j

        ! Initialize matrix.
        A = 0.0_sp
        do i = 1, n-1
            A(i, i+1) = m*1.0_sp
        enddo

        ! Reference with analytical exponential.
        Eref = eye(n, mold=1.0_sp)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        E = expm(A)

        ! Check result.
        err = norm(Eref - E, "inf")
        call check(error, err < (n**2)*epsilon(1.0_sp), "Analytical matrix exponential.")
        if (allocated(error)) return
        return        
    end subroutine test_expm_s    
    subroutine test_expm_d(error)
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer(ilp), parameter :: n = 5, m = 6
        ! Test matrix.
        real(dp) :: A(n, n), E(n, n), Eref(n, n)
        real(dp) :: err
        integer(ilp) :: i, j

        ! Initialize matrix.
        A = 0.0_dp
        do i = 1, n-1
            A(i, i+1) = m*1.0_dp
        enddo

        ! Reference with analytical exponential.
        Eref = eye(n, mold=1.0_dp)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        E = expm(A)

        ! Check result.
        err = norm(Eref - E, "inf")
        call check(error, err < (n**2)*epsilon(1.0_dp), "Analytical matrix exponential.")
        if (allocated(error)) return
        return        
    end subroutine test_expm_d    
    subroutine test_expm_c(error)
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer(ilp), parameter :: n = 5, m = 6
        ! Test matrix.
        complex(sp) :: A(n, n), E(n, n), Eref(n, n)
        real(sp) :: err
        integer(ilp) :: i, j

        ! Initialize matrix.
        A = 0.0_sp
        do i = 1, n-1
            A(i, i+1) = m*1.0_sp
        enddo

        ! Reference with analytical exponential.
        Eref = eye(n, mold=1.0_sp)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        E = expm(A)

        ! Check result.
        err = norm(Eref - E, "inf")
        call check(error, err < (n**2)*epsilon(1.0_sp), "Analytical matrix exponential.")
        if (allocated(error)) return
        return        
    end subroutine test_expm_c    
    subroutine test_expm_z(error)
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer(ilp), parameter :: n = 5, m = 6
        ! Test matrix.
        complex(dp) :: A(n, n), E(n, n), Eref(n, n)
        real(dp) :: err
        integer(ilp) :: i, j

        ! Initialize matrix.
        A = 0.0_dp
        do i = 1, n-1
            A(i, i+1) = m*1.0_dp
        enddo

        ! Reference with analytical exponential.
        Eref = eye(n, mold=1.0_dp)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        E = expm(A)

        ! Check result.
        err = norm(Eref - E, "inf")
        call check(error, err < (n**2)*epsilon(1.0_dp), "Analytical matrix exponential.")
        if (allocated(error)) return
        return        
    end subroutine test_expm_z    

    !> Test error handler.
    subroutine test_error_handling_expm_s(error)
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer(ilp), parameter :: n = 5, m = 6
        ! Test matrix.
        real(sp) :: A(n, n), E(n, n)
        type(linalg_state_type) :: err
        integer(ilp) :: i

        ! Initialize matrix.
        A = 0.0_sp
        do i = 1, n-1
            A(i, i+1) = m*1.0_sp
        enddo

        ! Compute matrix exponential.
        call matrix_exp(A, E, order=-1, err=err)
        ! Check result.
        call check(error, err%error(), "Negative Pade order")
        if (allocated(error)) return

        ! Compute matrix exponential.
        call matrix_exp(A(:n, :n-1), E, err=err)
        ! Check result.
        call check(error, err%error(), "Invalid matrix size")
        if (allocated(error)) return

        return        
    end subroutine test_error_handling_expm_s    
    subroutine test_error_handling_expm_d(error)
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer(ilp), parameter :: n = 5, m = 6
        ! Test matrix.
        real(dp) :: A(n, n), E(n, n)
        type(linalg_state_type) :: err
        integer(ilp) :: i

        ! Initialize matrix.
        A = 0.0_dp
        do i = 1, n-1
            A(i, i+1) = m*1.0_dp
        enddo

        ! Compute matrix exponential.
        call matrix_exp(A, E, order=-1, err=err)
        ! Check result.
        call check(error, err%error(), "Negative Pade order")
        if (allocated(error)) return

        ! Compute matrix exponential.
        call matrix_exp(A(:n, :n-1), E, err=err)
        ! Check result.
        call check(error, err%error(), "Invalid matrix size")
        if (allocated(error)) return

        return        
    end subroutine test_error_handling_expm_d    
    subroutine test_error_handling_expm_c(error)
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer(ilp), parameter :: n = 5, m = 6
        ! Test matrix.
        complex(sp) :: A(n, n), E(n, n)
        type(linalg_state_type) :: err
        integer(ilp) :: i

        ! Initialize matrix.
        A = 0.0_sp
        do i = 1, n-1
            A(i, i+1) = m*1.0_sp
        enddo

        ! Compute matrix exponential.
        call matrix_exp(A, E, order=-1, err=err)
        ! Check result.
        call check(error, err%error(), "Negative Pade order")
        if (allocated(error)) return

        ! Compute matrix exponential.
        call matrix_exp(A(:n, :n-1), E, err=err)
        ! Check result.
        call check(error, err%error(), "Invalid matrix size")
        if (allocated(error)) return

        return        
    end subroutine test_error_handling_expm_c    
    subroutine test_error_handling_expm_z(error)
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer(ilp), parameter :: n = 5, m = 6
        ! Test matrix.
        complex(dp) :: A(n, n), E(n, n)
        type(linalg_state_type) :: err
        integer(ilp) :: i

        ! Initialize matrix.
        A = 0.0_dp
        do i = 1, n-1
            A(i, i+1) = m*1.0_dp
        enddo

        ! Compute matrix exponential.
        call matrix_exp(A, E, order=-1, err=err)
        ! Check result.
        call check(error, err%error(), "Negative Pade order")
        if (allocated(error)) return

        ! Compute matrix exponential.
        call matrix_exp(A(:n, :n-1), E, err=err)
        ! Check result.
        call check(error, err%error(), "Invalid matrix size")
        if (allocated(error)) return

        return        
    end subroutine test_error_handling_expm_z    

end module test_linalg_expm

program test_expm
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_linalg_expm, only : test_expm_computation
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg_expm", test_expm_computation) &
        ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if
end program test_expm
