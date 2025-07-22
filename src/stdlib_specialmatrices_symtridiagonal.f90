submodule (stdlib_specialmatrices) symtridiagonal_matrices
    use stdlib_linalg_lapack, only: lagtm

    character(len=*), parameter :: this = "symtridiagonal matrices"
    contains

    !--------------------------------
    !-----                      -----
    !-----     CONSTRUCTORS     -----
    !-----                      -----
    !--------------------------------

    pure module function initialize_symtridiagonal_pure_sp(dv, ev) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        real(sp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(symtridiagonal_sp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if (size(ev, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector ev does not have the correct length.")
            call linalg_error_handling(err0)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dv = dv ; A%ev = ev
    end function

    pure module function initialize_constant_symtridiagonal_pure_sp(dv, ev, n) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        real(sp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(symtridiagonal_sp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: i
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        ! Matrix elements.
        A%dv = [(dv, i = 1, n-1)]
        A%ev = [(ev, i = 1, n-1)]
    end function

    module function initialize_symtridiagonal_impure_sp(dv, ev, err) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        real(sp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(symtridiagonal_sp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        if (size(ev, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector ev does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dv = dv ; A%ev = ev
    end function

    module function initialize_constant_symtridiagonal_impure_sp(dv, ev, n, err) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        real(sp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(symtridiagonal_sp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: i
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        ! Matrix elements.
        A%dv = [(dv, i = 1, n-1)]
        A%ev = [(ev, i = 1, n)]
    end function
    pure module function initialize_symtridiagonal_pure_dp(dv, ev) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        real(dp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(symtridiagonal_dp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if (size(ev, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector ev does not have the correct length.")
            call linalg_error_handling(err0)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dv = dv ; A%ev = ev
    end function

    pure module function initialize_constant_symtridiagonal_pure_dp(dv, ev, n) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        real(dp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(symtridiagonal_dp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: i
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        ! Matrix elements.
        A%dv = [(dv, i = 1, n-1)]
        A%ev = [(ev, i = 1, n-1)]
    end function

    module function initialize_symtridiagonal_impure_dp(dv, ev, err) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        real(dp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(symtridiagonal_dp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        if (size(ev, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector ev does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dv = dv ; A%ev = ev
    end function

    module function initialize_constant_symtridiagonal_impure_dp(dv, ev, n, err) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        real(dp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(symtridiagonal_dp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: i
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        ! Matrix elements.
        A%dv = [(dv, i = 1, n-1)]
        A%ev = [(ev, i = 1, n)]
    end function
    pure module function initialize_symtridiagonal_pure_csp(dv, ev) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(sp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(symtridiagonal_csp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if (size(ev, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector ev does not have the correct length.")
            call linalg_error_handling(err0)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dv = dv ; A%ev = ev
    end function

    pure module function initialize_constant_symtridiagonal_pure_csp(dv, ev, n) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        complex(sp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(symtridiagonal_csp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: i
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        ! Matrix elements.
        A%dv = [(dv, i = 1, n-1)]
        A%ev = [(ev, i = 1, n-1)]
    end function

    module function initialize_symtridiagonal_impure_csp(dv, ev, err) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(sp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(symtridiagonal_csp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        if (size(ev, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector ev does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dv = dv ; A%ev = ev
    end function

    module function initialize_constant_symtridiagonal_impure_csp(dv, ev, n, err) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        complex(sp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(symtridiagonal_csp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: i
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        ! Matrix elements.
        A%dv = [(dv, i = 1, n-1)]
        A%ev = [(ev, i = 1, n)]
    end function
    pure module function initialize_symtridiagonal_pure_cdp(dv, ev) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(dp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(symtridiagonal_cdp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if (size(ev, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector ev does not have the correct length.")
            call linalg_error_handling(err0)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dv = dv ; A%ev = ev
    end function

    pure module function initialize_constant_symtridiagonal_pure_cdp(dv, ev, n) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        complex(dp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(symtridiagonal_cdp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: i
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        ! Matrix elements.
        A%dv = [(dv, i = 1, n-1)]
        A%ev = [(ev, i = 1, n-1)]
    end function

    module function initialize_symtridiagonal_impure_cdp(dv, ev, err) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(dp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(symtridiagonal_cdp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        if (size(ev, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector ev does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dv = dv ; A%ev = ev
    end function

    module function initialize_constant_symtridiagonal_impure_cdp(dv, ev, n, err) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        complex(dp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(symtridiagonal_cdp_type) :: A
        !! Corresponding symtridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: i
        type(linalg_state_type) :: err0

        ! Description of the matrix.
        A%n = n
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        ! Matrix elements.
        A%dv = [(dv, i = 1, n-1)]
        A%ev = [(ev, i = 1, n)]
    end function

    !-----------------------------------------
    !-----                               -----
    !-----     MATRIX-VECTOR PRODUCT     -----
    !-----                               -----
    !-----------------------------------------

    !! spmv_tridiag
    module subroutine spmv_symtridiag_1d_sp(A, x, y, alpha, beta, op)
        type(symtridiagonal_sp_type), intent(in) :: A
        real(sp), intent(in), contiguous, target :: x(:)
        real(sp), intent(inout), contiguous, target :: y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(sp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_
        real(sp), pointer :: xmat(:, :), ymat(:, :)

        ! Deal with optional arguments.
        alpha_ = 1.0_sp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_sp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_sp
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call lagtm(op_, n, nrhs, alpha_, A%ev, A%dv, A%ev, xmat, ldx, beta_, ymat, ldy)
    end subroutine
    module subroutine spmv_symtridiag_2d_sp(A, x, y, alpha, beta, op)
        type(symtridiagonal_sp_type), intent(in) :: A
        real(sp), intent(in), contiguous, target :: x(:,:)
        real(sp), intent(inout), contiguous, target :: y(:,:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(sp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_

        ! Deal with optional arguments.
        alpha_ = 1.0_sp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_sp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_sp
        nrhs =  size(x, dim=2, kind=ilp) 

        call lagtm(op_, n, nrhs, alpha_, A%ev, A%dv, A%ev, x, ldx, beta_, y, ldy)
    end subroutine
    module subroutine spmv_symtridiag_1d_dp(A, x, y, alpha, beta, op)
        type(symtridiagonal_dp_type), intent(in) :: A
        real(dp), intent(in), contiguous, target :: x(:)
        real(dp), intent(inout), contiguous, target :: y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(dp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_
        real(dp), pointer :: xmat(:, :), ymat(:, :)

        ! Deal with optional arguments.
        alpha_ = 1.0_dp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_dp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_dp
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call lagtm(op_, n, nrhs, alpha_, A%ev, A%dv, A%ev, xmat, ldx, beta_, ymat, ldy)
    end subroutine
    module subroutine spmv_symtridiag_2d_dp(A, x, y, alpha, beta, op)
        type(symtridiagonal_dp_type), intent(in) :: A
        real(dp), intent(in), contiguous, target :: x(:,:)
        real(dp), intent(inout), contiguous, target :: y(:,:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(dp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_

        ! Deal with optional arguments.
        alpha_ = 1.0_dp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_dp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_dp
        nrhs =  size(x, dim=2, kind=ilp) 

        call lagtm(op_, n, nrhs, alpha_, A%ev, A%dv, A%ev, x, ldx, beta_, y, ldy)
    end subroutine
    module subroutine spmv_symtridiag_1d_csp(A, x, y, alpha, beta, op)
        type(symtridiagonal_csp_type), intent(in) :: A
        complex(sp), intent(in), contiguous, target :: x(:)
        complex(sp), intent(inout), contiguous, target :: y(:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(sp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_
        complex(sp), pointer :: xmat(:, :), ymat(:, :)

        ! Deal with optional arguments.
        alpha_ = 1.0_sp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_sp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_sp
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call lagtm(op_, n, nrhs, alpha_, A%ev, A%dv, A%ev, xmat, ldx, beta_, ymat, ldy)
    end subroutine
    module subroutine spmv_symtridiag_2d_csp(A, x, y, alpha, beta, op)
        type(symtridiagonal_csp_type), intent(in) :: A
        complex(sp), intent(in), contiguous, target :: x(:,:)
        complex(sp), intent(inout), contiguous, target :: y(:,:)
        real(sp), intent(in), optional :: alpha
        real(sp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(sp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_

        ! Deal with optional arguments.
        alpha_ = 1.0_sp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_sp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_sp
        nrhs =  size(x, dim=2, kind=ilp) 

        call lagtm(op_, n, nrhs, alpha_, A%ev, A%dv, A%ev, x, ldx, beta_, y, ldy)
    end subroutine
    module subroutine spmv_symtridiag_1d_cdp(A, x, y, alpha, beta, op)
        type(symtridiagonal_cdp_type), intent(in) :: A
        complex(dp), intent(in), contiguous, target :: x(:)
        complex(dp), intent(inout), contiguous, target :: y(:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(dp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_
        complex(dp), pointer :: xmat(:, :), ymat(:, :)

        ! Deal with optional arguments.
        alpha_ = 1.0_dp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_dp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_dp
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call lagtm(op_, n, nrhs, alpha_, A%ev, A%dv, A%ev, xmat, ldx, beta_, ymat, ldy)
    end subroutine
    module subroutine spmv_symtridiag_2d_cdp(A, x, y, alpha, beta, op)
        type(symtridiagonal_cdp_type), intent(in) :: A
        complex(dp), intent(in), contiguous, target :: x(:,:)
        complex(dp), intent(inout), contiguous, target :: y(:,:)
        real(dp), intent(in), optional :: alpha
        real(dp), intent(in), optional :: beta
        character(1), intent(in), optional :: op

        ! Internal variables.
        real(dp) :: alpha_, beta_
        integer(ilp) :: n, nrhs, ldx, ldy
        character(1) :: op_

        ! Deal with optional arguments.
        alpha_ = 1.0_dp ; if (present(alpha)) alpha_ = alpha
        beta_  = 0.0_dp ; if (present(beta))  beta_  = beta
        op_    = "N"        ; if (present(op))    op_    = op

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_dp
        nrhs =  size(x, dim=2, kind=ilp) 

        call lagtm(op_, n, nrhs, alpha_, A%ev, A%dv, A%ev, x, ldx, beta_, y, ldy)
    end subroutine

    !-------------------------------------
    !-----                           -----
    !-----     UTILITY FUNCTIONS     -----
    !-----                           -----
    !-------------------------------------

    pure module function symtridiagonal_to_dense_sp(A) result(B)
        !! Convert a `symtridiagonal` matrix to its dense representation.
        type(symtridiagonal_sp_type), intent(in) :: A
        !! Input symtridiagonal matrix.
        real(sp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
        allocate(B(n, n), source=zero_sp)
        B(1, 1) = A%dv(1) ; B(1, 2) = A%ev(1)
        do concurrent (i=2:n-1)
            B(i, i-1) = A%ev(i-1)
            B(i, i) = A%dv(i)
            B(i, i+1) = A%ev(i)
        enddo
        B(n, n-1) = A%ev(n-1) ; B(n, n) = A%dv(n)
        end associate
    end function
    pure module function symtridiagonal_to_dense_dp(A) result(B)
        !! Convert a `symtridiagonal` matrix to its dense representation.
        type(symtridiagonal_dp_type), intent(in) :: A
        !! Input symtridiagonal matrix.
        real(dp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
        allocate(B(n, n), source=zero_dp)
        B(1, 1) = A%dv(1) ; B(1, 2) = A%ev(1)
        do concurrent (i=2:n-1)
            B(i, i-1) = A%ev(i-1)
            B(i, i) = A%dv(i)
            B(i, i+1) = A%ev(i)
        enddo
        B(n, n-1) = A%ev(n-1) ; B(n, n) = A%dv(n)
        end associate
    end function
    pure module function symtridiagonal_to_dense_csp(A) result(B)
        !! Convert a `symtridiagonal` matrix to its dense representation.
        type(symtridiagonal_csp_type), intent(in) :: A
        !! Input symtridiagonal matrix.
        complex(sp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
        allocate(B(n, n), source=zero_csp)
        B(1, 1) = A%dv(1) ; B(1, 2) = A%ev(1)
        do concurrent (i=2:n-1)
            B(i, i-1) = A%ev(i-1)
            B(i, i) = A%dv(i)
            B(i, i+1) = A%ev(i)
        enddo
        B(n, n-1) = A%ev(n-1) ; B(n, n) = A%dv(n)
        end associate
    end function
    pure module function symtridiagonal_to_dense_cdp(A) result(B)
        !! Convert a `symtridiagonal` matrix to its dense representation.
        type(symtridiagonal_cdp_type), intent(in) :: A
        !! Input symtridiagonal matrix.
        complex(dp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
        allocate(B(n, n), source=zero_cdp)
        B(1, 1) = A%dv(1) ; B(1, 2) = A%ev(1)
        do concurrent (i=2:n-1)
            B(i, i-1) = A%ev(i-1)
            B(i, i) = A%dv(i)
            B(i, i+1) = A%ev(i)
        enddo
        B(n, n-1) = A%ev(n-1) ; B(n, n) = A%dv(n)
        end associate
    end function

    pure module function transpose_symtridiagonal_sp(A) result(B)
        type(symtridiagonal_sp_type), intent(in) :: A
        !! Input matrix.
        type(symtridiagonal_sp_type) :: B
        B = A
    end function
    pure module function transpose_symtridiagonal_dp(A) result(B)
        type(symtridiagonal_dp_type), intent(in) :: A
        !! Input matrix.
        type(symtridiagonal_dp_type) :: B
        B = A
    end function
    pure module function transpose_symtridiagonal_csp(A) result(B)
        type(symtridiagonal_csp_type), intent(in) :: A
        !! Input matrix.
        type(symtridiagonal_csp_type) :: B
        B = A
    end function
    pure module function transpose_symtridiagonal_cdp(A) result(B)
        type(symtridiagonal_cdp_type), intent(in) :: A
        !! Input matrix.
        type(symtridiagonal_cdp_type) :: B
        B = A
    end function

    pure module function hermitian_symtridiagonal_sp(A) result(B)
        type(symtridiagonal_sp_type), intent(in) :: A
        !! Input matrix.
        type(symtridiagonal_sp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
    end function
    pure module function hermitian_symtridiagonal_dp(A) result(B)
        type(symtridiagonal_dp_type), intent(in) :: A
        !! Input matrix.
        type(symtridiagonal_dp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
    end function
    pure module function hermitian_symtridiagonal_csp(A) result(B)
        type(symtridiagonal_csp_type), intent(in) :: A
        !! Input matrix.
        type(symtridiagonal_csp_type) :: B
        B = symtridiagonal(A%dv, conjg(A%ev))
    end function
    pure module function hermitian_symtridiagonal_cdp(A) result(B)
        type(symtridiagonal_cdp_type), intent(in) :: A
        !! Input matrix.
        type(symtridiagonal_cdp_type) :: B
        B = symtridiagonal(A%dv, conjg(A%ev))
    end function

    pure module function scalar_multiplication_symtridiagonal_sp(alpha, A) result(B)
        real(sp), intent(in) :: alpha
        type(symtridiagonal_sp_type), intent(in) :: A
        type(symtridiagonal_sp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
        B%dv = alpha*B%dv; B%ev = alpha*B%ev
    end function

    pure module function scalar_multiplication_bis_symtridiagonal_sp(A, alpha) result(B)
        type(symtridiagonal_sp_type), intent(in) :: A
        real(sp), intent(in) :: alpha
        type(symtridiagonal_sp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
        B%dv = alpha*B%dv; B%ev = alpha*B%ev
    end function
    pure module function scalar_multiplication_symtridiagonal_dp(alpha, A) result(B)
        real(dp), intent(in) :: alpha
        type(symtridiagonal_dp_type), intent(in) :: A
        type(symtridiagonal_dp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
        B%dv = alpha*B%dv; B%ev = alpha*B%ev
    end function

    pure module function scalar_multiplication_bis_symtridiagonal_dp(A, alpha) result(B)
        type(symtridiagonal_dp_type), intent(in) :: A
        real(dp), intent(in) :: alpha
        type(symtridiagonal_dp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
        B%dv = alpha*B%dv; B%ev = alpha*B%ev
    end function
    pure module function scalar_multiplication_symtridiagonal_csp(alpha, A) result(B)
        complex(sp), intent(in) :: alpha
        type(symtridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
        B%dv = alpha*B%dv; B%ev = alpha*B%ev
    end function

    pure module function scalar_multiplication_bis_symtridiagonal_csp(A, alpha) result(B)
        type(symtridiagonal_csp_type), intent(in) :: A
        complex(sp), intent(in) :: alpha
        type(symtridiagonal_csp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
        B%dv = alpha*B%dv; B%ev = alpha*B%ev
    end function
    pure module function scalar_multiplication_symtridiagonal_cdp(alpha, A) result(B)
        complex(dp), intent(in) :: alpha
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
        B%dv = alpha*B%dv; B%ev = alpha*B%ev
    end function

    pure module function scalar_multiplication_bis_symtridiagonal_cdp(A, alpha) result(B)
        type(symtridiagonal_cdp_type), intent(in) :: A
        complex(dp), intent(in) :: alpha
        type(symtridiagonal_cdp_type) :: B
        B = symtridiagonal(A%dv, A%ev)
        B%dv = alpha*B%dv; B%ev = alpha*B%ev
    end function

    pure module function matrix_add_symtridiagonal_sp(A, B) result(C)
        type(symtridiagonal_sp_type), intent(in) :: A
        type(symtridiagonal_sp_type), intent(in) :: B
        type(symtridiagonal_sp_type) :: C
        C = symtridiagonal(A%dv, A%ev)
        C%dv = C%dv + B%dv; C%ev = C%ev + B%ev
    end function

    pure module function matrix_sub_symtridiagonal_sp(A, B) result(C)
        type(symtridiagonal_sp_type), intent(in) :: A
        type(symtridiagonal_sp_type), intent(in) :: B
        type(symtridiagonal_sp_type) :: C
        C = symtridiagonal(A%dv, A%ev)
        C%dv = C%dv - B%dv; C%ev = C%ev - B%ev
    end function
    pure module function matrix_add_symtridiagonal_dp(A, B) result(C)
        type(symtridiagonal_dp_type), intent(in) :: A
        type(symtridiagonal_dp_type), intent(in) :: B
        type(symtridiagonal_dp_type) :: C
        C = symtridiagonal(A%dv, A%ev)
        C%dv = C%dv + B%dv; C%ev = C%ev + B%ev
    end function

    pure module function matrix_sub_symtridiagonal_dp(A, B) result(C)
        type(symtridiagonal_dp_type), intent(in) :: A
        type(symtridiagonal_dp_type), intent(in) :: B
        type(symtridiagonal_dp_type) :: C
        C = symtridiagonal(A%dv, A%ev)
        C%dv = C%dv - B%dv; C%ev = C%ev - B%ev
    end function
    pure module function matrix_add_symtridiagonal_csp(A, B) result(C)
        type(symtridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type), intent(in) :: B
        type(symtridiagonal_csp_type) :: C
        C = symtridiagonal(A%dv, A%ev)
        C%dv = C%dv + B%dv; C%ev = C%ev + B%ev
    end function

    pure module function matrix_sub_symtridiagonal_csp(A, B) result(C)
        type(symtridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type), intent(in) :: B
        type(symtridiagonal_csp_type) :: C
        C = symtridiagonal(A%dv, A%ev)
        C%dv = C%dv - B%dv; C%ev = C%ev - B%ev
    end function
    pure module function matrix_add_symtridiagonal_cdp(A, B) result(C)
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type), intent(in) :: B
        type(symtridiagonal_cdp_type) :: C
        C = symtridiagonal(A%dv, A%ev)
        C%dv = C%dv + B%dv; C%ev = C%ev + B%ev
    end function

    pure module function matrix_sub_symtridiagonal_cdp(A, B) result(C)
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type), intent(in) :: B
        type(symtridiagonal_cdp_type) :: C
        C = symtridiagonal(A%dv, A%ev)
        C%dv = C%dv - B%dv; C%ev = C%ev - B%ev
    end function

end submodule
