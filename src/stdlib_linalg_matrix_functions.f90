submodule (stdlib_linalg) stdlib_linalg_matrix_functions
    use stdlib_linalg_constants
    use stdlib_linalg_blas, only: gemm
    use stdlib_linalg_lapack, only: gesv
    use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_ERROR, &
         LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR
    implicit none

    real(sp), parameter :: zero_s = 0._sp
    real(sp), parameter :: one_s  = 1._sp
    real(dp), parameter :: zero_d = 0._dp
    real(dp), parameter :: one_d  = 1._dp
    complex(sp), parameter :: zero_c = (0._sp, 0._sp)
    complex(sp), parameter :: one_c  = (1._sp, 0._sp)
    complex(dp), parameter :: zero_z = (0._dp, 0._dp)
    complex(dp), parameter :: one_z  = (1._dp, 0._dp)

contains

    module function stdlib_expm_s(A, order, err) result(E)
        !> Input matrix A(n, n).
        real(sp), intent(in) :: A(:, :)
        !> [optional] Order of the Pade approximation.
        integer(ilp), optional, intent(in) :: order
        !> [optional] State return flag.
        type(linalg_state_type), optional, intent(out) :: err
        !> Exponential of the input matrix E = exp(A).
        real(sp), allocatable :: E(:, :)

        ! Internal variables.
        real(sp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(sp)        :: a_norm, c
        integer(ilp)        :: m, n, ee, k, s, order_, i, j
        logical(lk)         :: p
        character(len=*), parameter :: this = "expm"
        type(linalg_state_type)     :: err0

        ! Deal with optional args.
        order_ = 10 ; if (present(order)) order_ = order

        ! Problem's dimension.
        m = size(A, 1) ; n = size(A, 2)

        if (m /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Invalid matrix size A=',[m, n])
            call linalg_error_handling(err0, err)
            return
        else if (order_ < 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, 'Order of Pade approximation &
                                    needs to be positive, order=', order_)
            call linalg_error_handling(err0, err)
            return
        endif

        ! Compute the L-infinity norm.
        a_norm = mnorm(A, "inf")

        ! Determine scaling factor for the matrix.
        ee = int(log(a_norm) / log(2.0_sp)) + 1
        s  = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A/2.0_sp**s ; X = A2

        ! First step of the Pade approximation.
        c = 0.5_sp
        allocate (E, source=A2) ; allocate (Q, source=A2)
        do concurrent(i=1:n, j=1:n)
            E(i, j) =  c*E(i, j) ; if (i == j) E(i, j) = 1.0_sp + E(i, j) ! E = I + c*A2
            Q(i, j) = -c*Q(i, j) ; if (i == j) Q(i, j) = 1.0_sp + Q(i, j) ! Q = I - c*A2
        enddo

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, order_
            c = c * (order_ - k + 1) / (k * (2*order_ - k + 1))
            X = matmul(A2, X)
            do concurrent(i=1:n, j=1:n)
                E(i, j) = E(i, j) + c*X(i, j)       ! E = E + c*X
            enddo
            if (p) then
                do concurrent(i=1:n, j=1:n)
                    Q(i, j) = Q(i, j) + c*X(i, j)   ! Q = Q + c*X
                enddo
            else
                do concurrent(i=1:n, j=1:n)
                    Q(i, j) = Q(i, j) - c*X(i, j)   ! Q = Q - c*X
                enddo
            endif
            p = .not. p
        enddo

        block
            integer(ilp) :: ipiv(n), info
            call gesv(n, n, Q, n, ipiv, E, n, info) ! E = inv(Q) @ E
            call handle_gesv_info(info, n, n, n, err0)
            call linalg_error_handling(err0, err)
        end block

        ! Matrix squaring.
        block
            real(sp) :: E_tmp(n, n)
            do k = 1, s
                E_tmp = E
                call gemm("N", "N", n, n, n, one_s, E_tmp, n, E_tmp, n, zero_s, E, n)
            enddo
        end block
        return
    contains
        elemental subroutine handle_gesv_info(info,lda,n,nrhs,err)
            integer(ilp), intent(in) :: info,lda,n,nrhs
            type(linalg_state_type), intent(out) :: err
            ! Process output
            select case (info)
            case (0)
            ! Success
            case (-1)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid matrix size a=',[lda,n])
            case (-7)
                err = linalg_state_type(this,LINALG_ERROR,'invalid matrix size a=',[lda,n])
            case (1:)
                err = linalg_state_type(this,LINALG_ERROR,'singular matrix')
            case default
                err = linalg_state_type(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select
        end subroutine handle_gesv_info
    end function stdlib_expm_s
    module function stdlib_expm_d(A, order, err) result(E)
        !> Input matrix A(n, n).
        real(dp), intent(in) :: A(:, :)
        !> [optional] Order of the Pade approximation.
        integer(ilp), optional, intent(in) :: order
        !> [optional] State return flag.
        type(linalg_state_type), optional, intent(out) :: err
        !> Exponential of the input matrix E = exp(A).
        real(dp), allocatable :: E(:, :)

        ! Internal variables.
        real(dp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(dp)        :: a_norm, c
        integer(ilp)        :: m, n, ee, k, s, order_, i, j
        logical(lk)         :: p
        character(len=*), parameter :: this = "expm"
        type(linalg_state_type)     :: err0

        ! Deal with optional args.
        order_ = 10 ; if (present(order)) order_ = order

        ! Problem's dimension.
        m = size(A, 1) ; n = size(A, 2)

        if (m /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Invalid matrix size A=',[m, n])
            call linalg_error_handling(err0, err)
            return
        else if (order_ < 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, 'Order of Pade approximation &
                                    needs to be positive, order=', order_)
            call linalg_error_handling(err0, err)
            return
        endif

        ! Compute the L-infinity norm.
        a_norm = mnorm(A, "inf")

        ! Determine scaling factor for the matrix.
        ee = int(log(a_norm) / log(2.0_dp)) + 1
        s  = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A/2.0_dp**s ; X = A2

        ! First step of the Pade approximation.
        c = 0.5_dp
        allocate (E, source=A2) ; allocate (Q, source=A2)
        do concurrent(i=1:n, j=1:n)
            E(i, j) =  c*E(i, j) ; if (i == j) E(i, j) = 1.0_dp + E(i, j) ! E = I + c*A2
            Q(i, j) = -c*Q(i, j) ; if (i == j) Q(i, j) = 1.0_dp + Q(i, j) ! Q = I - c*A2
        enddo

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, order_
            c = c * (order_ - k + 1) / (k * (2*order_ - k + 1))
            X = matmul(A2, X)
            do concurrent(i=1:n, j=1:n)
                E(i, j) = E(i, j) + c*X(i, j)       ! E = E + c*X
            enddo
            if (p) then
                do concurrent(i=1:n, j=1:n)
                    Q(i, j) = Q(i, j) + c*X(i, j)   ! Q = Q + c*X
                enddo
            else
                do concurrent(i=1:n, j=1:n)
                    Q(i, j) = Q(i, j) - c*X(i, j)   ! Q = Q - c*X
                enddo
            endif
            p = .not. p
        enddo

        block
            integer(ilp) :: ipiv(n), info
            call gesv(n, n, Q, n, ipiv, E, n, info) ! E = inv(Q) @ E
            call handle_gesv_info(info, n, n, n, err0)
            call linalg_error_handling(err0, err)
        end block

        ! Matrix squaring.
        block
            real(dp) :: E_tmp(n, n)
            do k = 1, s
                E_tmp = E
                call gemm("N", "N", n, n, n, one_d, E_tmp, n, E_tmp, n, zero_d, E, n)
            enddo
        end block
        return
    contains
        elemental subroutine handle_gesv_info(info,lda,n,nrhs,err)
            integer(ilp), intent(in) :: info,lda,n,nrhs
            type(linalg_state_type), intent(out) :: err
            ! Process output
            select case (info)
            case (0)
            ! Success
            case (-1)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid matrix size a=',[lda,n])
            case (-7)
                err = linalg_state_type(this,LINALG_ERROR,'invalid matrix size a=',[lda,n])
            case (1:)
                err = linalg_state_type(this,LINALG_ERROR,'singular matrix')
            case default
                err = linalg_state_type(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select
        end subroutine handle_gesv_info
    end function stdlib_expm_d
    module function stdlib_expm_c(A, order, err) result(E)
        !> Input matrix A(n, n).
        complex(sp), intent(in) :: A(:, :)
        !> [optional] Order of the Pade approximation.
        integer(ilp), optional, intent(in) :: order
        !> [optional] State return flag.
        type(linalg_state_type), optional, intent(out) :: err
        !> Exponential of the input matrix E = exp(A).
        complex(sp), allocatable :: E(:, :)

        ! Internal variables.
        complex(sp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(sp)        :: a_norm, c
        integer(ilp)        :: m, n, ee, k, s, order_, i, j
        logical(lk)         :: p
        character(len=*), parameter :: this = "expm"
        type(linalg_state_type)     :: err0

        ! Deal with optional args.
        order_ = 10 ; if (present(order)) order_ = order

        ! Problem's dimension.
        m = size(A, 1) ; n = size(A, 2)

        if (m /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Invalid matrix size A=',[m, n])
            call linalg_error_handling(err0, err)
            return
        else if (order_ < 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, 'Order of Pade approximation &
                                    needs to be positive, order=', order_)
            call linalg_error_handling(err0, err)
            return
        endif

        ! Compute the L-infinity norm.
        a_norm = mnorm(A, "inf")

        ! Determine scaling factor for the matrix.
        ee = int(log(a_norm) / log(2.0_sp)) + 1
        s  = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A/2.0_sp**s ; X = A2

        ! First step of the Pade approximation.
        c = 0.5_sp
        allocate (E, source=A2) ; allocate (Q, source=A2)
        do concurrent(i=1:n, j=1:n)
            E(i, j) =  c*E(i, j) ; if (i == j) E(i, j) = 1.0_sp + E(i, j) ! E = I + c*A2
            Q(i, j) = -c*Q(i, j) ; if (i == j) Q(i, j) = 1.0_sp + Q(i, j) ! Q = I - c*A2
        enddo

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, order_
            c = c * (order_ - k + 1) / (k * (2*order_ - k + 1))
            X = matmul(A2, X)
            do concurrent(i=1:n, j=1:n)
                E(i, j) = E(i, j) + c*X(i, j)       ! E = E + c*X
            enddo
            if (p) then
                do concurrent(i=1:n, j=1:n)
                    Q(i, j) = Q(i, j) + c*X(i, j)   ! Q = Q + c*X
                enddo
            else
                do concurrent(i=1:n, j=1:n)
                    Q(i, j) = Q(i, j) - c*X(i, j)   ! Q = Q - c*X
                enddo
            endif
            p = .not. p
        enddo

        block
            integer(ilp) :: ipiv(n), info
            call gesv(n, n, Q, n, ipiv, E, n, info) ! E = inv(Q) @ E
            call handle_gesv_info(info, n, n, n, err0)
            call linalg_error_handling(err0, err)
        end block

        ! Matrix squaring.
        block
            complex(sp) :: E_tmp(n, n)
            do k = 1, s
                E_tmp = E
                call gemm("N", "N", n, n, n, one_c, E_tmp, n, E_tmp, n, zero_c, E, n)
            enddo
        end block
        return
    contains
        elemental subroutine handle_gesv_info(info,lda,n,nrhs,err)
            integer(ilp), intent(in) :: info,lda,n,nrhs
            type(linalg_state_type), intent(out) :: err
            ! Process output
            select case (info)
            case (0)
            ! Success
            case (-1)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid matrix size a=',[lda,n])
            case (-7)
                err = linalg_state_type(this,LINALG_ERROR,'invalid matrix size a=',[lda,n])
            case (1:)
                err = linalg_state_type(this,LINALG_ERROR,'singular matrix')
            case default
                err = linalg_state_type(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select
        end subroutine handle_gesv_info
    end function stdlib_expm_c
    module function stdlib_expm_z(A, order, err) result(E)
        !> Input matrix A(n, n).
        complex(dp), intent(in) :: A(:, :)
        !> [optional] Order of the Pade approximation.
        integer(ilp), optional, intent(in) :: order
        !> [optional] State return flag.
        type(linalg_state_type), optional, intent(out) :: err
        !> Exponential of the input matrix E = exp(A).
        complex(dp), allocatable :: E(:, :)

        ! Internal variables.
        complex(dp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(dp)        :: a_norm, c
        integer(ilp)        :: m, n, ee, k, s, order_, i, j
        logical(lk)         :: p
        character(len=*), parameter :: this = "expm"
        type(linalg_state_type)     :: err0

        ! Deal with optional args.
        order_ = 10 ; if (present(order)) order_ = order

        ! Problem's dimension.
        m = size(A, 1) ; n = size(A, 2)

        if (m /= n) then
            err0 = linalg_state_type(this,LINALG_VALUE_ERROR,'Invalid matrix size A=',[m, n])
            call linalg_error_handling(err0, err)
            return
        else if (order_ < 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, 'Order of Pade approximation &
                                    needs to be positive, order=', order_)
            call linalg_error_handling(err0, err)
            return
        endif

        ! Compute the L-infinity norm.
        a_norm = mnorm(A, "inf")

        ! Determine scaling factor for the matrix.
        ee = int(log(a_norm) / log(2.0_dp)) + 1
        s  = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A/2.0_dp**s ; X = A2

        ! First step of the Pade approximation.
        c = 0.5_dp
        allocate (E, source=A2) ; allocate (Q, source=A2)
        do concurrent(i=1:n, j=1:n)
            E(i, j) =  c*E(i, j) ; if (i == j) E(i, j) = 1.0_dp + E(i, j) ! E = I + c*A2
            Q(i, j) = -c*Q(i, j) ; if (i == j) Q(i, j) = 1.0_dp + Q(i, j) ! Q = I - c*A2
        enddo

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, order_
            c = c * (order_ - k + 1) / (k * (2*order_ - k + 1))
            X = matmul(A2, X)
            do concurrent(i=1:n, j=1:n)
                E(i, j) = E(i, j) + c*X(i, j)       ! E = E + c*X
            enddo
            if (p) then
                do concurrent(i=1:n, j=1:n)
                    Q(i, j) = Q(i, j) + c*X(i, j)   ! Q = Q + c*X
                enddo
            else
                do concurrent(i=1:n, j=1:n)
                    Q(i, j) = Q(i, j) - c*X(i, j)   ! Q = Q - c*X
                enddo
            endif
            p = .not. p
        enddo

        block
            integer(ilp) :: ipiv(n), info
            call gesv(n, n, Q, n, ipiv, E, n, info) ! E = inv(Q) @ E
            call handle_gesv_info(info, n, n, n, err0)
            call linalg_error_handling(err0, err)
        end block

        ! Matrix squaring.
        block
            complex(dp) :: E_tmp(n, n)
            do k = 1, s
                E_tmp = E
                call gemm("N", "N", n, n, n, one_z, E_tmp, n, E_tmp, n, zero_z, E, n)
            enddo
        end block
        return
    contains
        elemental subroutine handle_gesv_info(info,lda,n,nrhs,err)
            integer(ilp), intent(in) :: info,lda,n,nrhs
            type(linalg_state_type), intent(out) :: err
            ! Process output
            select case (info)
            case (0)
            ! Success
            case (-1)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err = linalg_state_type(this,LINALG_VALUE_ERROR,'invalid matrix size a=',[lda,n])
            case (-7)
                err = linalg_state_type(this,LINALG_ERROR,'invalid matrix size a=',[lda,n])
            case (1:)
                err = linalg_state_type(this,LINALG_ERROR,'singular matrix')
            case default
                err = linalg_state_type(this,LINALG_INTERNAL_ERROR,'catastrophic error')
            end select
        end subroutine handle_gesv_info
    end function stdlib_expm_z

end submodule stdlib_linalg_matrix_functions
