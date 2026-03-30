! Cholesky factorization of a matrix, based on LAPACK *POTRF functions
submodule (stdlib_linalg) stdlib_linalg_matrix_factorizations
    use stdlib_linalg_lapack, only: geqrt, gemqrt
    use stdlib_linalg_lapack_aux, only: handle_geqrt_info, handle_gemqrt_info
    implicit none

    character(len=*), parameter :: this = 'matrix_factorizations'

    contains

    !-----------------------------------------------------
    !-----     DERIVED-TYPE FOR QR FACTORIZATION     -----
    !-----------------------------------------------------

    pure module function stdlib_qrfact_s(A) result(F)
        real(sp), intent(in) :: A(:, :)
        !> Matrix to be factorized.
        type(qr_rsp_type) :: F

        !----- Internal variables -----
        integer(ilp), parameter :: min_nb = 36 ! Minimum block size.
        integer(ilp) :: m, n, nb, info
        real(sp), allocatable :: work(:)
        type(linalg_state_type) :: err

        !> Problem dimensions.
        m  = size(A, 1)         ! Number of rows.
        n  = size(A, 2)         ! Number of columns.
        nb = min(m, n, min_nb)  ! Block size.

        !> Copy matrix.
        F%data = A
        F%ldt = nb

        !> Allocate data.
        allocate(F%t(F%ldt, min(m, n)), F%work(nb*max(m, n)))

        !> QR factorization.
        call geqrt(m, n, nb, F%data, m, F%t, F%ldt, F%work, info)
        call handle_geqrt_info(this, info, m, n, nb, F%ldt, err)
        call linalg_error_handling(err)
    end function stdlib_qrfact_s

    module function get_qfactor_rsp(self) result(Q)
        class(qr_rsp_type), intent(inout) :: self
        !> Compact representation of the QR factorization.
        real(sp), allocatable :: Q(:, :)
        !> Orthogonal matrix.

        !----- Internal variables -----
        integer(ilp) :: i, j, k, info
        character(len=1), parameter :: side='L', trans='N'
        real(sp), parameter :: zero = 0.0_sp
        type(linalg_state_type) :: err
        
        associate( m => size(self%data, 1), n => size(self%data, 2), &
                   nb => self%ldt, ldv => size(self%data, 1), ldt => self%ldt, ldc => size(self%data, 1))
        !> Compute the Q matrix.
        k = min(m, n)
        Q = eye(m, k, mold=zero)
        call gemqrt(side, trans, m, k, k, &
                    nb, self%data(:, :k), ldv, self%t, ldt, Q, ldc, self%work, info)
        call handle_gemqrt_info(this, info, side, trans, m, n, k, nb, ldv, ldt, ldc, err)
        call linalg_error_handling(err)
        end associate
    end function get_qfactor_rsp

    pure module function get_rfactor_rsp(self) result(R)
        class(qr_rsp_type), intent(in) :: self
        !> Compact representation of the QR factorization.
        real(sp), allocatable :: R(:, :)
        !> Upper-triangular matrix.

        !----- Internal variables -----
        integer(ilp) :: i, j, k, m, n
        real(sp), parameter :: zero =  0.0_sp

        !> Matrix size.
        m = size(self%data, 1)
        n = size(self%data, 2)
        k = min(m, n)

        !> Allocate R matrix.
        allocate(R(k, n), source=zero)

        !> Fill-in matrix.
        do concurrent(i=1:k, j=1:n, i <= j)
            R(i, j) = self%data(i, j)
        enddo
    end function get_rfactor_rsp
    pure module function stdlib_qrfact_d(A) result(F)
        real(dp), intent(in) :: A(:, :)
        !> Matrix to be factorized.
        type(qr_rdp_type) :: F

        !----- Internal variables -----
        integer(ilp), parameter :: min_nb = 36 ! Minimum block size.
        integer(ilp) :: m, n, nb, info
        real(dp), allocatable :: work(:)
        type(linalg_state_type) :: err

        !> Problem dimensions.
        m  = size(A, 1)         ! Number of rows.
        n  = size(A, 2)         ! Number of columns.
        nb = min(m, n, min_nb)  ! Block size.

        !> Copy matrix.
        F%data = A
        F%ldt = nb

        !> Allocate data.
        allocate(F%t(F%ldt, min(m, n)), F%work(nb*max(m, n)))

        !> QR factorization.
        call geqrt(m, n, nb, F%data, m, F%t, F%ldt, F%work, info)
        call handle_geqrt_info(this, info, m, n, nb, F%ldt, err)
        call linalg_error_handling(err)
    end function stdlib_qrfact_d

    module function get_qfactor_rdp(self) result(Q)
        class(qr_rdp_type), intent(inout) :: self
        !> Compact representation of the QR factorization.
        real(dp), allocatable :: Q(:, :)
        !> Orthogonal matrix.

        !----- Internal variables -----
        integer(ilp) :: i, j, k, info
        character(len=1), parameter :: side='L', trans='N'
        real(dp), parameter :: zero = 0.0_dp
        type(linalg_state_type) :: err
        
        associate( m => size(self%data, 1), n => size(self%data, 2), &
                   nb => self%ldt, ldv => size(self%data, 1), ldt => self%ldt, ldc => size(self%data, 1))
        !> Compute the Q matrix.
        k = min(m, n)
        Q = eye(m, k, mold=zero)
        call gemqrt(side, trans, m, k, k, &
                    nb, self%data(:, :k), ldv, self%t, ldt, Q, ldc, self%work, info)
        call handle_gemqrt_info(this, info, side, trans, m, n, k, nb, ldv, ldt, ldc, err)
        call linalg_error_handling(err)
        end associate
    end function get_qfactor_rdp

    pure module function get_rfactor_rdp(self) result(R)
        class(qr_rdp_type), intent(in) :: self
        !> Compact representation of the QR factorization.
        real(dp), allocatable :: R(:, :)
        !> Upper-triangular matrix.

        !----- Internal variables -----
        integer(ilp) :: i, j, k, m, n
        real(dp), parameter :: zero =  0.0_dp

        !> Matrix size.
        m = size(self%data, 1)
        n = size(self%data, 2)
        k = min(m, n)

        !> Allocate R matrix.
        allocate(R(k, n), source=zero)

        !> Fill-in matrix.
        do concurrent(i=1:k, j=1:n, i <= j)
            R(i, j) = self%data(i, j)
        enddo
    end function get_rfactor_rdp
    pure module function stdlib_qrfact_c(A) result(F)
        complex(sp), intent(in) :: A(:, :)
        !> Matrix to be factorized.
        type(qr_csp_type) :: F

        !----- Internal variables -----
        integer(ilp), parameter :: min_nb = 36 ! Minimum block size.
        integer(ilp) :: m, n, nb, info
        complex(sp), allocatable :: work(:)
        type(linalg_state_type) :: err

        !> Problem dimensions.
        m  = size(A, 1)         ! Number of rows.
        n  = size(A, 2)         ! Number of columns.
        nb = min(m, n, min_nb)  ! Block size.

        !> Copy matrix.
        F%data = A
        F%ldt = nb

        !> Allocate data.
        allocate(F%t(F%ldt, min(m, n)), F%work(nb*max(m, n)))

        !> QR factorization.
        call geqrt(m, n, nb, F%data, m, F%t, F%ldt, F%work, info)
        call handle_geqrt_info(this, info, m, n, nb, F%ldt, err)
        call linalg_error_handling(err)
    end function stdlib_qrfact_c

    module function get_qfactor_csp(self) result(Q)
        class(qr_csp_type), intent(inout) :: self
        !> Compact representation of the QR factorization.
        complex(sp), allocatable :: Q(:, :)
        !> Orthogonal matrix.

        !----- Internal variables -----
        integer(ilp) :: i, j, k, info
        character(len=1), parameter :: side='L', trans='N'
        complex(sp), parameter :: zero = cmplx(0.0_sp, 0.0_sp, kind=sp)
        type(linalg_state_type) :: err
        
        associate( m => size(self%data, 1), n => size(self%data, 2), &
                   nb => self%ldt, ldv => size(self%data, 1), ldt => self%ldt, ldc => size(self%data, 1))
        !> Compute the Q matrix.
        k = min(m, n)
        Q = eye(m, k, mold=zero)
        call gemqrt(side, trans, m, k, k, &
                    nb, self%data(:, :k), ldv, self%t, ldt, Q, ldc, self%work, info)
        call handle_gemqrt_info(this, info, side, trans, m, n, k, nb, ldv, ldt, ldc, err)
        call linalg_error_handling(err)
        end associate
    end function get_qfactor_csp

    pure module function get_rfactor_csp(self) result(R)
        class(qr_csp_type), intent(in) :: self
        !> Compact representation of the QR factorization.
        complex(sp), allocatable :: R(:, :)
        !> Upper-triangular matrix.

        !----- Internal variables -----
        integer(ilp) :: i, j, k, m, n
        complex(sp), parameter :: zero = cmplx(0.0_sp, 0.0_sp, kind=sp) 

        !> Matrix size.
        m = size(self%data, 1)
        n = size(self%data, 2)
        k = min(m, n)

        !> Allocate R matrix.
        allocate(R(k, n), source=zero)

        !> Fill-in matrix.
        do concurrent(i=1:k, j=1:n, i <= j)
            R(i, j) = self%data(i, j)
        enddo
    end function get_rfactor_csp
    pure module function stdlib_qrfact_z(A) result(F)
        complex(dp), intent(in) :: A(:, :)
        !> Matrix to be factorized.
        type(qr_cdp_type) :: F

        !----- Internal variables -----
        integer(ilp), parameter :: min_nb = 36 ! Minimum block size.
        integer(ilp) :: m, n, nb, info
        complex(dp), allocatable :: work(:)
        type(linalg_state_type) :: err

        !> Problem dimensions.
        m  = size(A, 1)         ! Number of rows.
        n  = size(A, 2)         ! Number of columns.
        nb = min(m, n, min_nb)  ! Block size.

        !> Copy matrix.
        F%data = A
        F%ldt = nb

        !> Allocate data.
        allocate(F%t(F%ldt, min(m, n)), F%work(nb*max(m, n)))

        !> QR factorization.
        call geqrt(m, n, nb, F%data, m, F%t, F%ldt, F%work, info)
        call handle_geqrt_info(this, info, m, n, nb, F%ldt, err)
        call linalg_error_handling(err)
    end function stdlib_qrfact_z

    module function get_qfactor_cdp(self) result(Q)
        class(qr_cdp_type), intent(inout) :: self
        !> Compact representation of the QR factorization.
        complex(dp), allocatable :: Q(:, :)
        !> Orthogonal matrix.

        !----- Internal variables -----
        integer(ilp) :: i, j, k, info
        character(len=1), parameter :: side='L', trans='N'
        complex(dp), parameter :: zero = cmplx(0.0_dp, 0.0_dp, kind=dp)
        type(linalg_state_type) :: err
        
        associate( m => size(self%data, 1), n => size(self%data, 2), &
                   nb => self%ldt, ldv => size(self%data, 1), ldt => self%ldt, ldc => size(self%data, 1))
        !> Compute the Q matrix.
        k = min(m, n)
        Q = eye(m, k, mold=zero)
        call gemqrt(side, trans, m, k, k, &
                    nb, self%data(:, :k), ldv, self%t, ldt, Q, ldc, self%work, info)
        call handle_gemqrt_info(this, info, side, trans, m, n, k, nb, ldv, ldt, ldc, err)
        call linalg_error_handling(err)
        end associate
    end function get_qfactor_cdp

    pure module function get_rfactor_cdp(self) result(R)
        class(qr_cdp_type), intent(in) :: self
        !> Compact representation of the QR factorization.
        complex(dp), allocatable :: R(:, :)
        !> Upper-triangular matrix.

        !----- Internal variables -----
        integer(ilp) :: i, j, k, m, n
        complex(dp), parameter :: zero = cmplx(0.0_dp, 0.0_dp, kind=dp) 

        !> Matrix size.
        m = size(self%data, 1)
        n = size(self%data, 2)
        k = min(m, n)

        !> Allocate R matrix.
        allocate(R(k, n), source=zero)

        !> Fill-in matrix.
        do concurrent(i=1:k, j=1:n, i <= j)
            R(i, j) = self%data(i, j)
        enddo
    end function get_rfactor_cdp
end submodule stdlib_linalg_matrix_factorizations

