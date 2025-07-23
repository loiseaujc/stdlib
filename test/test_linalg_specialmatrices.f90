module test_specialmatrices
    use testdrive, only : new_unittest, unittest_type, error_type, check, skip_test
    use stdlib_kinds
    use stdlib_constants
    use stdlib_linalg, only: hermitian
    use stdlib_linalg_state, only: linalg_state_type
    use stdlib_math, only: all_close
    use stdlib_specialmatrices
    implicit none

contains


    !> Collect all exported unit tests
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            ! new_unittest('tridiagonal arithmetic', test_tridiagonal_arithmetic), &
            new_unittest('tridiagonal spmv kernel', test_tridiagonal_spmv) &
            ! new_unittest('tridiagonal error handling', test_tridiagonal_error_handling), &
            ! new_unittest('symtridiagonal arithmetic', test_symtridiagonal_arithmetic), &
            ! new_unittest('symtridiagonal spmv kernel', test_symtridiagonal_spmv), &
            ! ! new_unittest('hermtridiagonal arithmetic', test_hermtridiagonal_arithmetic), &
            ! new_unittest('hermtridiagonal spmv kernel', test_hermtridiagonal_spmv) &
        ]
    end subroutine

    !----------------------------------------
    !-----                              -----
    !-----     TRIDIAGONAL MATRICES     -----
    !-----                              -----
    !----------------------------------------

    subroutine test_tridiagonal_arithmetic(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(tridiagonal_sp_type) :: A, B, C
            type(symtridiagonal_sp_type) :: D
            real(sp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :), Dmat(:, :)
            real(sp), allocatable :: dl(:), dv(:), du(:), ev(:)
            real(sp), parameter :: alpha = 2.0_dp
            integer :: i, j

            ! Initialize A and B matrices.
            allocate(dl(n-1), dv(n), du(n-1), ev(n-1))
            call random_number(dl) ; call random_number(dv) ; call random_number(du)
            A = tridiagonal(dl, dv, du) ; Amat = dense(A)
            
            call random_number(dl) ; call random_number(dv) ; call random_number(du)
            B = tridiagonal(dl, dv, du) ; Bmat = dense(B)
            
            call random_number(dv) ; call random_number(ev)
            D = symtridiagonal(dv, ev) ; Dmat = dense(D)

       
            ! Matrix addition.
            C = A + B ; Cmat = dense(C) ! Tridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return

            C = A + D ; Cmat = dense(C) ! Tridiag + SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Dmat), .true.)
            if (allocated(error)) return

            C = D + A ; Cmat = dense(C) ! SymTridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Dmat), .true.)
            if (allocated(error)) return


            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C) ! Tridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return

            C = A - D ; Cmat = dense(C) ! Tridiag - SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Dmat), .true.)
            if (allocated(error)) return

            C = D - A ; Cmat = dense(C) ! SymTridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Dmat-Amat), .true.)
            if (allocated(error)) return



            ! Matrix scalar multiplication
            C = alpha * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = A *alpha ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(tridiagonal_dp_type) :: A, B, C
            type(symtridiagonal_dp_type) :: D
            real(dp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :), Dmat(:, :)
            real(dp), allocatable :: dl(:), dv(:), du(:), ev(:)
            real(dp), parameter :: alpha = 2.0_dp
            integer :: i, j

            ! Initialize A and B matrices.
            allocate(dl(n-1), dv(n), du(n-1), ev(n-1))
            call random_number(dl) ; call random_number(dv) ; call random_number(du)
            A = tridiagonal(dl, dv, du) ; Amat = dense(A)
            
            call random_number(dl) ; call random_number(dv) ; call random_number(du)
            B = tridiagonal(dl, dv, du) ; Bmat = dense(B)
            
            call random_number(dv) ; call random_number(ev)
            D = symtridiagonal(dv, ev) ; Dmat = dense(D)

       
            ! Matrix addition.
            C = A + B ; Cmat = dense(C) ! Tridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return

            C = A + D ; Cmat = dense(C) ! Tridiag + SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Dmat), .true.)
            if (allocated(error)) return

            C = D + A ; Cmat = dense(C) ! SymTridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Dmat), .true.)
            if (allocated(error)) return


            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C) ! Tridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return

            C = A - D ; Cmat = dense(C) ! Tridiag - SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Dmat), .true.)
            if (allocated(error)) return

            C = D - A ; Cmat = dense(C) ! SymTridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Dmat-Amat), .true.)
            if (allocated(error)) return



            ! Matrix scalar multiplication
            C = alpha * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = A *alpha ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(tridiagonal_csp_type) :: A, B, C
            type(symtridiagonal_csp_type) :: D
            complex(sp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :), Dmat(:, :)
            complex(sp), allocatable :: dl(:), dv(:), du(:), ev(:)
            complex(sp), parameter :: alpha = 2.0_dp
            integer :: i, j
            type(hermtridiagonal_csp_type) :: E
            complex(sp), allocatable :: Emat(:, :)
            real(wp), allocatable :: data(:, :)

            ! Initialize A and B matrices.
            allocate(dl(n-1), dv(n), du(n-1), ev(n-1))
            allocate(data(n, 2))
            call random_number(data) ; dl%re = data(:n-1, 1) ; dl%im = data(:n-1, 2)
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; du%re = data(:n-1, 1) ; du%im = data(:n-1, 2)
            A = tridiagonal(dl, dv, du) ; Amat = dense(A)
            
            call random_number(data) ; dl%re = data(:n-1, 1) ; dl%im = data(:n-1, 2)
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; du%re = data(:n-1, 1) ; du%im = data(:n-1, 2)
            B = tridiagonal(dl, dv, du) ; Bmat = dense(B)
            
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            D = symtridiagonal(dv, ev) ; Dmat = dense(D)

            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            E = hermtridiagonal(dv, ev) ; Emat = dense(E)
       
            ! Matrix addition.
            C = A + B ; Cmat = dense(C) ! Tridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return

            C = A + D ; Cmat = dense(C) ! Tridiag + SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Dmat), .true.)
            if (allocated(error)) return

            C = D + A ; Cmat = dense(C) ! SymTridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Dmat), .true.)
            if (allocated(error)) return

            C = A + E ; Cmat = dense(C) ! Tridiag + HermTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Emat), .true.)
            if (allocated(error)) return

            C = E + A ; Cmat = dense(C) ! HermTridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Emat), .true.)
            if (allocated(error)) return

            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C) ! Tridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return

            C = A - D ; Cmat = dense(C) ! Tridiag - SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Dmat), .true.)
            if (allocated(error)) return

            C = D - A ; Cmat = dense(C) ! SymTridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Dmat-Amat), .true.)
            if (allocated(error)) return

            C = A - E ; Cmat = dense(C) ! Tridiag - HermTridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Emat), .true.)
            if (allocated(error)) return

            C = E - A ; Cmat = dense(C) ! HermTridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Emat-Amat), .true.)
            if (allocated(error)) return


            ! Matrix scalar multiplication
            C = alpha * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = A *alpha ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(tridiagonal_cdp_type) :: A, B, C
            type(symtridiagonal_cdp_type) :: D
            complex(dp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :), Dmat(:, :)
            complex(dp), allocatable :: dl(:), dv(:), du(:), ev(:)
            complex(dp), parameter :: alpha = 2.0_dp
            integer :: i, j
            type(hermtridiagonal_cdp_type) :: E
            complex(dp), allocatable :: Emat(:, :)
            real(wp), allocatable :: data(:, :)

            ! Initialize A and B matrices.
            allocate(dl(n-1), dv(n), du(n-1), ev(n-1))
            allocate(data(n, 2))
            call random_number(data) ; dl%re = data(:n-1, 1) ; dl%im = data(:n-1, 2)
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; du%re = data(:n-1, 1) ; du%im = data(:n-1, 2)
            A = tridiagonal(dl, dv, du) ; Amat = dense(A)
            
            call random_number(data) ; dl%re = data(:n-1, 1) ; dl%im = data(:n-1, 2)
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; du%re = data(:n-1, 1) ; du%im = data(:n-1, 2)
            B = tridiagonal(dl, dv, du) ; Bmat = dense(B)
            
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            D = symtridiagonal(dv, ev) ; Dmat = dense(D)

            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            E = hermtridiagonal(dv, ev) ; Emat = dense(E)
       
            ! Matrix addition.
            C = A + B ; Cmat = dense(C) ! Tridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return

            C = A + D ; Cmat = dense(C) ! Tridiag + SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Dmat), .true.)
            if (allocated(error)) return

            C = D + A ; Cmat = dense(C) ! SymTridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Dmat), .true.)
            if (allocated(error)) return

            C = A + E ; Cmat = dense(C) ! Tridiag + HermTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Emat), .true.)
            if (allocated(error)) return

            C = E + A ; Cmat = dense(C) ! HermTridiag + Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Emat), .true.)
            if (allocated(error)) return

            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C) ! Tridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return

            C = A - D ; Cmat = dense(C) ! Tridiag - SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Dmat), .true.)
            if (allocated(error)) return

            C = D - A ; Cmat = dense(C) ! SymTridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Dmat-Amat), .true.)
            if (allocated(error)) return

            C = A - E ; Cmat = dense(C) ! Tridiag - HermTridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Emat), .true.)
            if (allocated(error)) return

            C = E - A ; Cmat = dense(C) ! HermTridiag - Tridiag = Tridiag
            call check(error, all_close(Cmat, Emat-Amat), .true.)
            if (allocated(error)) return


            ! Matrix scalar multiplication
            C = alpha * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = A *alpha ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return
        end block
    end subroutine

    subroutine test_tridiagonal_spmv(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(tridiagonal_sp_type) :: A
            real(sp), allocatable :: Amat(:,:), dl(:), dv(:), du(:)
            real(sp), allocatable :: x(:)
            real(sp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(dl(n-1), dv(n), du(n-1))
            call random_number(dl) ; call random_number(dv) ; call random_number(du)
            A = tridiagonal(dl, dv, du) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source = 0.0_wp)   ; call random_number(x)
            allocate(y1(n), source = 0.0_wp)  ; allocate(y2(n), source=0.0_wp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(tridiagonal_dp_type) :: A
            real(dp), allocatable :: Amat(:,:), dl(:), dv(:), du(:)
            real(dp), allocatable :: x(:)
            real(dp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(dl(n-1), dv(n), du(n-1))
            call random_number(dl) ; call random_number(dv) ; call random_number(du)
            A = tridiagonal(dl, dv, du) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source = 0.0_wp)   ; call random_number(x)
            allocate(y1(n), source = 0.0_wp)  ; allocate(y2(n), source=0.0_wp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

        end block
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(tridiagonal_csp_type) :: A
            complex(sp), allocatable :: Amat(:,:), dl(:), dv(:), du(:)
            real(wp), allocatable :: data(:, :)
            complex(sp), allocatable :: x(:)
            complex(sp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(dl(n-1), dv(n), du(n-1))
            allocate(data(n, 2), source=0.0_wp)
            call random_number(data) ; dl = cmplx(data(:n-1, 1), data(:n-1, 2), kind=wp)
            call random_number(data) ; dv = cmplx(data(:, 1), data(:, 2), kind=wp)
            call random_number(data) ; du = cmplx(data(:n-1, 1), data(:n-1, 2), kind=wp)
            A = tridiagonal(dl, dv, du) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source=zero_csp)
            call random_number(data) ; x = cmplx(data(:, 1), data(:, 2), kind=wp)
            allocate(y1(n), source = zero_csp)  ; allocate(y2(n), source=zero_csp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.H @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(hermitian(Amat), x) ; call spmv(A, x, y2, op="H")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(tridiagonal_cdp_type) :: A
            complex(dp), allocatable :: Amat(:,:), dl(:), dv(:), du(:)
            real(wp), allocatable :: data(:, :)
            complex(dp), allocatable :: x(:)
            complex(dp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(dl(n-1), dv(n), du(n-1))
            allocate(data(n, 2), source=0.0_wp)
            call random_number(data) ; dl = cmplx(data(:n-1, 1), data(:n-1, 2), kind=wp)
            call random_number(data) ; dv = cmplx(data(:, 1), data(:, 2), kind=wp)
            call random_number(data) ; du = cmplx(data(:n-1, 1), data(:n-1, 2), kind=wp)
            A = tridiagonal(dl, dv, du) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source=zero_cdp)
            call random_number(data) ; x = cmplx(data(:, 1), data(:, 2), kind=wp)
            allocate(y1(n), source = zero_cdp)  ; allocate(y2(n), source=zero_cdp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.H @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(hermitian(Amat), x) ; call spmv(A, x, y2, op="H")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return
        end block
    end subroutine

    subroutine test_tridiagonal_error_handling(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(tridiagonal_sp_type) :: A
            real(sp), allocatable :: dl(:), dv(:), du(:)
            type(linalg_state_type) :: state
            integer :: i

            !> Test constructor from arrays.
            dl = [(1.0_wp, i = 1, n-2)] ; du = dl
            dv = [(2.0_wp, i = 1, n)]
            A = tridiagonal(dl, dv, du, state)
            call check(error, state%ok(), .false.)
            if (allocated(error)) return

            !> Test contructor from constants.
            A = tridiagonal(dl(1), dv(1), du(1), -n, state)
            call check(error, state%ok(), .false.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(tridiagonal_dp_type) :: A
            real(dp), allocatable :: dl(:), dv(:), du(:)
            type(linalg_state_type) :: state
            integer :: i

            !> Test constructor from arrays.
            dl = [(1.0_wp, i = 1, n-2)] ; du = dl
            dv = [(2.0_wp, i = 1, n)]
            A = tridiagonal(dl, dv, du, state)
            call check(error, state%ok(), .false.)
            if (allocated(error)) return

            !> Test contructor from constants.
            A = tridiagonal(dl(1), dv(1), du(1), -n, state)
            call check(error, state%ok(), .false.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(tridiagonal_csp_type) :: A
            complex(sp), allocatable :: dl(:), dv(:), du(:)
            type(linalg_state_type) :: state
            integer :: i

            !> Test constructor from arrays.
            dl = [(1.0_wp, i = 1, n-2)] ; du = dl
            dv = [(2.0_wp, i = 1, n)]
            A = tridiagonal(dl, dv, du, state)
            call check(error, state%ok(), .false.)
            if (allocated(error)) return

            !> Test contructor from constants.
            A = tridiagonal(dl(1), dv(1), du(1), -n, state)
            call check(error, state%ok(), .false.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(tridiagonal_cdp_type) :: A
            complex(dp), allocatable :: dl(:), dv(:), du(:)
            type(linalg_state_type) :: state
            integer :: i

            !> Test constructor from arrays.
            dl = [(1.0_wp, i = 1, n-2)] ; du = dl
            dv = [(2.0_wp, i = 1, n)]
            A = tridiagonal(dl, dv, du, state)
            call check(error, state%ok(), .false.)
            if (allocated(error)) return

            !> Test contructor from constants.
            A = tridiagonal(dl(1), dv(1), du(1), -n, state)
            call check(error, state%ok(), .false.)
            if (allocated(error)) return
        end block
    end subroutine

    !--------------------------------------------------
    !-----                                        -----
    !-----     SYMMETRIC TRIDIAGONAL MATRICES     -----
    !-----                                        -----
    !--------------------------------------------------

    subroutine test_symtridiagonal_arithmetic(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(symtridiagonal_sp_type) :: A, B, C
            type(tridiagonal_sp_type) :: D
            real(sp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :)
            real(sp), allocatable :: dv(:), ev(:)
            real(sp), parameter :: alpha = 2.0_dp
            integer :: i, j

            ! Initialize A and B matrices.
            allocate(ev(n-1), dv(n))
            call random_number(ev) ; call random_number(dv)
            A = symtridiagonal(dv, ev) ; Amat = dense(A)

            call random_number(ev) ; call random_number(dv)
            B = symtridiagonal(dv, ev) ; Bmat = dense(B)


            ! Matrix addition.
            C = A + B ; Cmat = dense(C) ! SymTridiag + SymTridiag = SymTridiag
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return


            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C) ! SymTridiag - SymTridiag = SymTridiag
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return


            ! Matrix scalar multiplication
            C = alpha * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = A *alpha ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(symtridiagonal_dp_type) :: A, B, C
            type(tridiagonal_dp_type) :: D
            real(dp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :)
            real(dp), allocatable :: dv(:), ev(:)
            real(dp), parameter :: alpha = 2.0_dp
            integer :: i, j

            ! Initialize A and B matrices.
            allocate(ev(n-1), dv(n))
            call random_number(ev) ; call random_number(dv)
            A = symtridiagonal(dv, ev) ; Amat = dense(A)

            call random_number(ev) ; call random_number(dv)
            B = symtridiagonal(dv, ev) ; Bmat = dense(B)


            ! Matrix addition.
            C = A + B ; Cmat = dense(C) ! SymTridiag + SymTridiag = SymTridiag
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return


            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C) ! SymTridiag - SymTridiag = SymTridiag
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return


            ! Matrix scalar multiplication
            C = alpha * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = A *alpha ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(symtridiagonal_csp_type) :: A, B, C
            type(tridiagonal_csp_type) :: D
            complex(sp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :)
            complex(sp), allocatable :: dv(:), ev(:)
            complex(sp), parameter :: alpha = 2.0_dp
            integer :: i, j
            type(hermtridiagonal_csp_type) :: E
            complex(sp), allocatable :: Emat(:, :)
            real(wp), allocatable :: data(:, :)

            ! Initialize A and B matrices.
            allocate(ev(n-1), dv(n))
            allocate(data(n, 2))
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            A = symtridiagonal(dv, ev) ; Amat = dense(A)

            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            B = symtridiagonal(dv, ev) ; Bmat = dense(B)

            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            E = hermtridiagonal(dv, ev) ; Emat = dense(E)

            ! Matrix addition.
            C = A + B ; Cmat = dense(C) ! SymTridiag + SymTridiag = SymTridiag
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return

            D = A + E ; Cmat = dense(D) ! SymTridiag + HermTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Emat), .true.)
            if (allocated(error)) return

            D = E + A ; Cmat = dense(D) ! HermTridiag + SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Emat), .true.)
            if (allocated(error)) return

            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C) ! SymTridiag - SymTridiag = SymTridiag
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return

            D = A - E ; Cmat = dense(D) ! SymTridiag - HermTridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Emat), .true.)
            if (allocated(error)) return

            D = E - A ; Cmat = dense(D) ! HermTridiag - SymTridiag = Tridiag
            call check(error, all_close(Cmat, Emat-Amat), .true.)
            if (allocated(error)) return

            ! Matrix scalar multiplication
            C = alpha * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = A *alpha ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(symtridiagonal_cdp_type) :: A, B, C
            type(tridiagonal_cdp_type) :: D
            complex(dp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :)
            complex(dp), allocatable :: dv(:), ev(:)
            complex(dp), parameter :: alpha = 2.0_dp
            integer :: i, j
            type(hermtridiagonal_cdp_type) :: E
            complex(dp), allocatable :: Emat(:, :)
            real(wp), allocatable :: data(:, :)

            ! Initialize A and B matrices.
            allocate(ev(n-1), dv(n))
            allocate(data(n, 2))
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            A = symtridiagonal(dv, ev) ; Amat = dense(A)

            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            B = symtridiagonal(dv, ev) ; Bmat = dense(B)

            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            E = hermtridiagonal(dv, ev) ; Emat = dense(E)

            ! Matrix addition.
            C = A + B ; Cmat = dense(C) ! SymTridiag + SymTridiag = SymTridiag
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return

            D = A + E ; Cmat = dense(D) ! SymTridiag + HermTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Emat), .true.)
            if (allocated(error)) return

            D = E + A ; Cmat = dense(D) ! HermTridiag + SymTridiag = Tridiag
            call check(error, all_close(Cmat, Amat+Emat), .true.)
            if (allocated(error)) return

            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C) ! SymTridiag - SymTridiag = SymTridiag
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return

            D = A - E ; Cmat = dense(D) ! SymTridiag - HermTridiag = Tridiag
            call check(error, all_close(Cmat, Amat-Emat), .true.)
            if (allocated(error)) return

            D = E - A ; Cmat = dense(D) ! HermTridiag - SymTridiag = Tridiag
            call check(error, all_close(Cmat, Emat-Amat), .true.)
            if (allocated(error)) return

            ! Matrix scalar multiplication
            C = alpha * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = A *alpha ; Cmat = dense(C)
            call check(error, all_close(Cmat, alpha * Amat), .true.)
            if (allocated(error)) return
        end block
    end subroutine

    subroutine test_symtridiagonal_spmv(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(symtridiagonal_sp_type) :: A
            real(sp), allocatable :: Amat(:,:), dv(:), ev(:)
            real(sp), allocatable :: x(:)
            real(sp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(ev(n-1), dv(n))
            call random_number(dv) ; call random_number(ev)
            A = symtridiagonal(dv, ev) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source = 0.0_wp)   ; call random_number(x)
            allocate(y1(n), source = 0.0_wp)  ; allocate(y2(n), source=0.0_wp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(symtridiagonal_dp_type) :: A
            real(dp), allocatable :: Amat(:,:), dv(:), ev(:)
            real(dp), allocatable :: x(:)
            real(dp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(ev(n-1), dv(n))
            call random_number(dv) ; call random_number(ev)
            A = symtridiagonal(dv, ev) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source = 0.0_wp)   ; call random_number(x)
            allocate(y1(n), source = 0.0_wp)  ; allocate(y2(n), source=0.0_wp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

        end block
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(symtridiagonal_csp_type) :: A
            complex(sp), allocatable :: Amat(:,:), dv(:), ev(:)
            real(wp), allocatable :: data(:, :)
            complex(sp), allocatable :: x(:)
            complex(sp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(ev(n-1), dv(n))
            allocate(data(n, 2), source=0.0_wp)
            call random_number(data) ; dv%re = data(:, 1) ; dv%im = data(:, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            A = symtridiagonal(dv, ev) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source=zero_csp)
            call random_number(data) ; x%re = data(:, 1) ; x%im = data(:, 2)
            allocate(y1(n), source = zero_csp)  ; allocate(y2(n), source=zero_csp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.H @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(hermitian(Amat), x) ; call spmv(A, x, y2, op="H")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(symtridiagonal_cdp_type) :: A
            complex(dp), allocatable :: Amat(:,:), dv(:), ev(:)
            real(wp), allocatable :: data(:, :)
            complex(dp), allocatable :: x(:)
            complex(dp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(ev(n-1), dv(n))
            allocate(data(n, 2), source=0.0_wp)
            call random_number(data) ; dv%re = data(:, 1) ; dv%im = data(:, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            A = symtridiagonal(dv, ev) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source=zero_cdp)
            call random_number(data) ; x%re = data(:, 1) ; x%im = data(:, 2)
            allocate(y1(n), source = zero_cdp)  ; allocate(y2(n), source=zero_cdp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.H @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(hermitian(Amat), x) ; call spmv(A, x, y2, op="H")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return
        end block
    end subroutine

    !--------------------------------------
    !-----                            -----
    !-----     HERMITIAN MATRICES     -----
    !-----                            -----
    !--------------------------------------

    subroutine test_hermtridiagonal_arithmetic(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(hermtridiagonal_csp_type) :: A, B, C
            type(tridiagonal_csp_type) :: D
            complex(sp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :), Dmat(:, :)
            complex(sp), allocatable :: dv(:), ev(:)
            complex(sp), parameter :: alpha = cmplx(2.0_wp, 2.0_wp, kind=wp)
            real(sp), parameter :: beta = 2.0_wp
            integer :: i, j
            real(wp), allocatable :: data(:, :)

            ! Initialize A and B matrices.
            allocate(ev(n-1), dv(n))
            allocate(data(n, 2))
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            A = hermtridiagonal(dv, ev) ; Amat = dense(A)

            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            B = hermtridiagonal(dv, ev) ; Bmat = dense(B)

            ! Matrix addition.
            C = A + B ; Cmat = dense(C)
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return

            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C)
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return

            ! Matrix scalar multiplication
            D = alpha * A ; Dmat = dense(D)
            call check(error, all_close(Dmat, alpha * Amat), .true.)
            if (allocated(error)) return

            D = A * alpha ; Dmat = dense(D)
            call check(error, all_close(Dmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = beta * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, beta * Amat), .true.)
            if (allocated(error)) return

            C = A * beta ; Cmat = dense(C)
            call check(error, all_close(Cmat, beta * Amat), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(hermtridiagonal_cdp_type) :: A, B, C
            type(tridiagonal_cdp_type) :: D
            complex(dp), allocatable :: Amat(:, :), Bmat(:, :), Cmat(:, :), Dmat(:, :)
            complex(dp), allocatable :: dv(:), ev(:)
            complex(dp), parameter :: alpha = cmplx(2.0_wp, 2.0_wp, kind=wp)
            real(dp), parameter :: beta = 2.0_wp
            integer :: i, j
            real(wp), allocatable :: data(:, :)

            ! Initialize A and B matrices.
            allocate(ev(n-1), dv(n))
            allocate(data(n, 2))
            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            A = hermtridiagonal(dv, ev) ; Amat = dense(A)

            call random_number(data) ; dv%re = data(:n, 1) ; dv%im = data(:n, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            B = hermtridiagonal(dv, ev) ; Bmat = dense(B)

            ! Matrix addition.
            C = A + B ; Cmat = dense(C)
            call check(error, all_close(Cmat, Amat+Bmat), .true.)
            if (allocated(error)) return

            ! Matrix subtraction.
            C = A - B ; Cmat = dense(C)
            call check(error, all_close(Cmat, Amat-Bmat), .true.)
            if (allocated(error)) return

            ! Matrix scalar multiplication
            D = alpha * A ; Dmat = dense(D)
            call check(error, all_close(Dmat, alpha * Amat), .true.)
            if (allocated(error)) return

            D = A * alpha ; Dmat = dense(D)
            call check(error, all_close(Dmat, alpha * Amat), .true.)
            if (allocated(error)) return

            C = beta * A ; Cmat = dense(C)
            call check(error, all_close(Cmat, beta * Amat), .true.)
            if (allocated(error)) return

            C = A * beta ; Cmat = dense(C)
            call check(error, all_close(Cmat, beta * Amat), .true.)
            if (allocated(error)) return
        end block
    end subroutine

    subroutine test_hermtridiagonal_spmv(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            integer, parameter :: n = 5
            type(hermtridiagonal_csp_type) :: A
            complex(sp), allocatable :: Amat(:,:), dv(:), ev(:)
            real(wp), allocatable :: data(:, :)
            complex(sp), allocatable :: x(:)
            complex(sp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(ev(n-1), dv(n))
            allocate(data(n, 2), source=0.0_wp)
            call random_number(data) ; dv%re = data(:, 1) ; dv%im = data(:, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            A = hermtridiagonal(dv, ev) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source=zero_csp)
            call random_number(data) ; x%re = data(:, 1) ; x%im = data(:, 2)
            allocate(y1(n), source = zero_csp)  ; allocate(y2(n), source=zero_csp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.H @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(hermitian(Amat), x) ; call spmv(A, x, y2, op="H")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            integer, parameter :: n = 5
            type(hermtridiagonal_cdp_type) :: A
            complex(dp), allocatable :: Amat(:,:), dv(:), ev(:)
            real(wp), allocatable :: data(:, :)
            complex(dp), allocatable :: x(:)
            complex(dp), allocatable :: y1(:), y2(:)

            ! Initialize matrix.
            allocate(ev(n-1), dv(n))
            allocate(data(n, 2), source=0.0_wp)
            call random_number(data) ; dv%re = data(:, 1) ; dv%im = data(:, 2)
            call random_number(data) ; ev%re = data(:n-1, 1) ; ev%im = data(:n-1, 2)
            A = hermtridiagonal(dv, ev) ; Amat = dense(A)

            ! Random vectors.
            allocate(x(n), source=zero_cdp)
            call random_number(data) ; x%re = data(:, 1) ; x%im = data(:, 2)
            allocate(y1(n), source = zero_cdp)  ; allocate(y2(n), source=zero_cdp)

            ! Test y = A @ x
            y1 = matmul(Amat, x) ; call spmv(A, x, y2)
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.T @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(transpose(Amat), x) ; call spmv(A, x, y2, op="T")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return

            ! Test y = A.H @ x
            y1 = 0.0_wp ; y2 = 0.0_wp
            y1 = matmul(hermitian(Amat), x) ; call spmv(A, x, y2, op="H")
            call check(error, all_close(y1, y2), .true.)
            if (allocated(error)) return
        end block
    end subroutine


end module


program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_specialmatrices, only : collect_suite
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("Tridiagonal matrices", collect_suite) &
        ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if
end program
