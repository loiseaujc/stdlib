submodule (stdlib_specialmatrices) tridiagonal_matrices
    use stdlib_linalg_lapack, only: lagtm

    character(len=*), parameter :: this = "tridiagonal matrices"
    contains

    !--------------------------------
    !-----                      -----
    !-----     CONSTRUCTORS     -----
    !-----                      -----
    !--------------------------------

    ! ----- Tridiagonal matrices -----

    pure module function initialize_tridiagonal_pure_sp(dl, dv, du) result(A)
        !! Construct a `tridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        real(sp), intent(in) :: dl(:), dv(:), du(:)
        !! tridiagonal matrix elements.
        type(tridiagonal_sp_type) :: A
        !! Corresponding tridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if (size(dl, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector dl does not have the correct length.")
            call linalg_error_handling(err0)
        endif
        if (size(du, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dl = dl ; A%dv = dv ; A%du = du
    end function

    pure module function initialize_constant_tridiagonal_pure_sp(dl, dv, du, n) result(A)
        !! Construct a `tridiagonal` matrix with constant elements.
        real(sp), intent(in) :: dl, dv, du
        !! tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(tridiagonal_sp_type) :: A
        !! Corresponding tridiagonal matrix.

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
        A%dl = [(dl, i = 1, n-1)]
        A%dv = [(dv, i = 1, n)]
        A%du = [(du, i = 1, n-1)]
    end function

    module function initialize_tridiagonal_impure_sp(dl, dv, du, err) result(A)
        !! Construct a `tridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        real(sp), intent(in) :: dl(:), dv(:), du(:)
        !! tridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(tridiagonal_sp_type) :: A
        !! Corresponding tridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        if (size(dl, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector dl does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif
        if (size(du, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dl = dl ; A%dv = dv ; A%du = du
    end function

    module function initialize_constant_tridiagonal_impure_sp(dl, dv, du, n, err) result(A)
        !! Construct a `tridiagonal` matrix with constant elements.
        real(sp), intent(in) :: dl, dv, du
        !! tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(tridiagonal_sp_type) :: A
        !! Corresponding tridiagonal matrix.

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
        A%dl = [(dl, i = 1, n-1)]
        A%dv = [(dv, i = 1, n)]
        A%du = [(du, i = 1, n-1)]
    end function
    pure module function initialize_tridiagonal_pure_dp(dl, dv, du) result(A)
        !! Construct a `tridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        real(dp), intent(in) :: dl(:), dv(:), du(:)
        !! tridiagonal matrix elements.
        type(tridiagonal_dp_type) :: A
        !! Corresponding tridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if (size(dl, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector dl does not have the correct length.")
            call linalg_error_handling(err0)
        endif
        if (size(du, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dl = dl ; A%dv = dv ; A%du = du
    end function

    pure module function initialize_constant_tridiagonal_pure_dp(dl, dv, du, n) result(A)
        !! Construct a `tridiagonal` matrix with constant elements.
        real(dp), intent(in) :: dl, dv, du
        !! tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(tridiagonal_dp_type) :: A
        !! Corresponding tridiagonal matrix.

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
        A%dl = [(dl, i = 1, n-1)]
        A%dv = [(dv, i = 1, n)]
        A%du = [(du, i = 1, n-1)]
    end function

    module function initialize_tridiagonal_impure_dp(dl, dv, du, err) result(A)
        !! Construct a `tridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        real(dp), intent(in) :: dl(:), dv(:), du(:)
        !! tridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(tridiagonal_dp_type) :: A
        !! Corresponding tridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        if (size(dl, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector dl does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif
        if (size(du, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dl = dl ; A%dv = dv ; A%du = du
    end function

    module function initialize_constant_tridiagonal_impure_dp(dl, dv, du, n, err) result(A)
        !! Construct a `tridiagonal` matrix with constant elements.
        real(dp), intent(in) :: dl, dv, du
        !! tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(tridiagonal_dp_type) :: A
        !! Corresponding tridiagonal matrix.

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
        A%dl = [(dl, i = 1, n-1)]
        A%dv = [(dv, i = 1, n)]
        A%du = [(du, i = 1, n-1)]
    end function
    pure module function initialize_tridiagonal_pure_csp(dl, dv, du) result(A)
        !! Construct a `tridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(sp), intent(in) :: dl(:), dv(:), du(:)
        !! tridiagonal matrix elements.
        type(tridiagonal_csp_type) :: A
        !! Corresponding tridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if (size(dl, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector dl does not have the correct length.")
            call linalg_error_handling(err0)
        endif
        if (size(du, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dl = dl ; A%dv = dv ; A%du = du
    end function

    pure module function initialize_constant_tridiagonal_pure_csp(dl, dv, du, n) result(A)
        !! Construct a `tridiagonal` matrix with constant elements.
        complex(sp), intent(in) :: dl, dv, du
        !! tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(tridiagonal_csp_type) :: A
        !! Corresponding tridiagonal matrix.

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
        A%dl = [(dl, i = 1, n-1)]
        A%dv = [(dv, i = 1, n)]
        A%du = [(du, i = 1, n-1)]
    end function

    module function initialize_tridiagonal_impure_csp(dl, dv, du, err) result(A)
        !! Construct a `tridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(sp), intent(in) :: dl(:), dv(:), du(:)
        !! tridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(tridiagonal_csp_type) :: A
        !! Corresponding tridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        if (size(dl, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector dl does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif
        if (size(du, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dl = dl ; A%dv = dv ; A%du = du
    end function

    module function initialize_constant_tridiagonal_impure_csp(dl, dv, du, n, err) result(A)
        !! Construct a `tridiagonal` matrix with constant elements.
        complex(sp), intent(in) :: dl, dv, du
        !! tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(tridiagonal_csp_type) :: A
        !! Corresponding tridiagonal matrix.

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
        A%dl = [(dl, i = 1, n-1)]
        A%dv = [(dv, i = 1, n)]
        A%du = [(du, i = 1, n-1)]
    end function
    pure module function initialize_tridiagonal_pure_cdp(dl, dv, du) result(A)
        !! Construct a `tridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(dp), intent(in) :: dl(:), dv(:), du(:)
        !! tridiagonal matrix elements.
        type(tridiagonal_cdp_type) :: A
        !! Corresponding tridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0)
        endif
        if (size(dl, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector dl does not have the correct length.")
            call linalg_error_handling(err0)
        endif
        if (size(du, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dl = dl ; A%dv = dv ; A%du = du
    end function

    pure module function initialize_constant_tridiagonal_pure_cdp(dl, dv, du, n) result(A)
        !! Construct a `tridiagonal` matrix with constant elements.
        complex(dp), intent(in) :: dl, dv, du
        !! tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(tridiagonal_cdp_type) :: A
        !! Corresponding tridiagonal matrix.

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
        A%dl = [(dl, i = 1, n-1)]
        A%dv = [(dv, i = 1, n)]
        A%du = [(du, i = 1, n-1)]
    end function

    module function initialize_tridiagonal_impure_cdp(dl, dv, du, err) result(A)
        !! Construct a `tridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(dp), intent(in) :: dl(:), dv(:), du(:)
        !! tridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(tridiagonal_cdp_type) :: A
        !! Corresponding tridiagonal matrix.

        ! Internal variables.
        integer(ilp) :: n
        type(linalg_state_type) :: err0

        ! Sanity check.
        n = size(dv, kind=ilp)
        if (n <= 0) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Matrix size needs to be positive, n = ", n, ".")
            call linalg_error_handling(err0, err)
        endif
        if (size(dl, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector dl does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif
        if (size(du, kind=ilp) /= n-1) then
            err0 = linalg_state_type(this, LINALG_VALUE_ERROR, "Vector du does not have the correct length.")
            call linalg_error_handling(err0, err)
        endif

        ! Description of the matrix.
        A%n = n
        ! Matrix elements.
        A%dl = dl ; A%dv = dv ; A%du = du
    end function

    module function initialize_constant_tridiagonal_impure_cdp(dl, dv, du, n, err) result(A)
        !! Construct a `tridiagonal` matrix with constant elements.
        complex(dp), intent(in) :: dl, dv, du
        !! tridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(tridiagonal_cdp_type) :: A
        !! Corresponding tridiagonal matrix.

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
        A%dl = [(dl, i = 1, n-1)]
        A%dv = [(dv, i = 1, n)]
        A%du = [(du, i = 1, n-1)]
    end function

    !----- Symmetric Tridiagonal matrices -----

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
        A%dl = ev ; A%dv = dv ; A%du = ev
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
        A%dl = [(ev, i = 1, n-1)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n-1)]
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
        A%dl = ev ; A%dv = dv ; A%du = ev
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
        A%dl = [(ev, i = 1, n)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n)]
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
        A%dl = ev ; A%dv = dv ; A%du = ev
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
        A%dl = [(ev, i = 1, n-1)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n-1)]
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
        A%dl = ev ; A%dv = dv ; A%du = ev
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
        A%dl = [(ev, i = 1, n)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n)]
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
        A%dl = ev ; A%dv = dv ; A%du = ev
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
        A%dl = [(ev, i = 1, n-1)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n-1)]
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
        A%dl = ev ; A%dv = dv ; A%du = ev
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
        A%dl = [(ev, i = 1, n)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n)]
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
        A%dl = ev ; A%dv = dv ; A%du = ev
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
        A%dl = [(ev, i = 1, n-1)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n-1)]
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
        A%dl = ev ; A%dv = dv ; A%du = ev
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
        A%dl = [(ev, i = 1, n)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n)]
    end function

    !----- Hermitian Tridiagonal matrices -----

    pure module function initialize_hermtridiagonal_pure_csp(dv, ev) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(sp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(hermtridiagonal_csp_type) :: A
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
        A%dl = conjg(ev) ; A%dv = dv ; A%du = ev
    end function

    pure module function initialize_constant_hermtridiagonal_pure_csp(dv, ev, n) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        complex(sp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(hermtridiagonal_csp_type) :: A
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
        A%dl = [(conjg(ev), i = 1, n-1)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n-1)]
    end function

    module function initialize_hermtridiagonal_impure_csp(dv, ev, err) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(sp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(hermtridiagonal_csp_type) :: A
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
        A%dl = conjg(ev) ; A%dv = dv ; A%du = ev
    end function

    module function initialize_constant_hermtridiagonal_impure_csp(dv, ev, n, err) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        complex(sp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(hermtridiagonal_csp_type) :: A
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
        A%dl = [(conjg(ev), i = 1, n)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n)]
    end function
    pure module function initialize_hermtridiagonal_pure_cdp(dv, ev) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(dp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(hermtridiagonal_cdp_type) :: A
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
        A%dl = conjg(ev) ; A%dv = dv ; A%du = ev
    end function

    pure module function initialize_constant_hermtridiagonal_pure_cdp(dv, ev, n) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        complex(dp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(hermtridiagonal_cdp_type) :: A
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
        A%dl = [(conjg(ev), i = 1, n-1)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n-1)]
    end function

    module function initialize_hermtridiagonal_impure_cdp(dv, ev, err) result(A)
        !! Construct a `symtridiagonal` matrix from the rank-1 arrays
        !! `dl`, `dv` and `du`.
        complex(dp), intent(in) :: dv(:), ev(:)
        !! symtridiagonal matrix elements.
        type(linalg_state_type), intent(out) :: err
        !! Error handling.
        type(hermtridiagonal_cdp_type) :: A
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
        A%dl = conjg(ev) ; A%dv = dv ; A%du = ev
    end function

    module function initialize_constant_hermtridiagonal_impure_cdp(dv, ev, n, err) result(A)
        !! Construct a `symtridiagonal` matrix with constant elements.
        complex(dp), intent(in) :: dv, ev
        !! symtridiagonal matrix elements.
        integer(ilp), intent(in) :: n
        !! Matrix dimension.
        type(linalg_state_type), intent(out) :: err
        !! Error handling
        type(hermtridiagonal_cdp_type) :: A
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
        A%dl = [(conjg(ev), i = 1, n)]
        A%dv = [(dv, i = 1, n-1)]
        A%du = [(ev, i = 1, n)]
    end function

    !-----------------------------------------
    !-----                               -----
    !-----     MATRIX-VECTOR PRODUCT     -----
    !-----                               -----
    !-----------------------------------------

    !! spmv_tridiag
    module subroutine spmv_tridiag_1d_sp(A, x, y, alpha, beta, op)
        class(tridiagonal_sp_type), intent(in) :: A
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
        if (op_ == "H") op_ = "C"

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_sp
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call lagtm(op_, n, nrhs, alpha_, A%dl, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
        nullify(xmat, ymat)
    end subroutine
    module subroutine spmv_tridiag_2d_sp(A, x, y, alpha, beta, op)
        class(tridiagonal_sp_type), intent(in) :: A
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
        if (op_ == "H") op_ = "C"

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_sp
        nrhs =  size(x, dim=2, kind=ilp) 

        call lagtm(op_, n, nrhs, alpha_, A%dl, A%dv, A%du, x, ldx, beta_, y, ldy)
    end subroutine
    module subroutine spmv_tridiag_1d_dp(A, x, y, alpha, beta, op)
        class(tridiagonal_dp_type), intent(in) :: A
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
        if (op_ == "H") op_ = "C"

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_dp
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call lagtm(op_, n, nrhs, alpha_, A%dl, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
        nullify(xmat, ymat)
    end subroutine
    module subroutine spmv_tridiag_2d_dp(A, x, y, alpha, beta, op)
        class(tridiagonal_dp_type), intent(in) :: A
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
        if (op_ == "H") op_ = "C"

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_dp
        nrhs =  size(x, dim=2, kind=ilp) 

        call lagtm(op_, n, nrhs, alpha_, A%dl, A%dv, A%du, x, ldx, beta_, y, ldy)
    end subroutine
    module subroutine spmv_tridiag_1d_csp(A, x, y, alpha, beta, op)
        class(tridiagonal_csp_type), intent(in) :: A
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
        if (op_ == "H") op_ = "C"

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_sp
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call lagtm(op_, n, nrhs, alpha_, A%dl, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
        nullify(xmat, ymat)
    end subroutine
    module subroutine spmv_tridiag_2d_csp(A, x, y, alpha, beta, op)
        class(tridiagonal_csp_type), intent(in) :: A
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
        if (op_ == "H") op_ = "C"

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_sp
        nrhs =  size(x, dim=2, kind=ilp) 

        call lagtm(op_, n, nrhs, alpha_, A%dl, A%dv, A%du, x, ldx, beta_, y, ldy)
    end subroutine
    module subroutine spmv_tridiag_1d_cdp(A, x, y, alpha, beta, op)
        class(tridiagonal_cdp_type), intent(in) :: A
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
        if (op_ == "H") op_ = "C"

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_dp
        nrhs =  1 

        ! Pointer trick.
        xmat(1:n, 1:nrhs) => x ; ymat(1:n, 1:nrhs) => y
        call lagtm(op_, n, nrhs, alpha_, A%dl, A%dv, A%du, xmat, ldx, beta_, ymat, ldy)
        nullify(xmat, ymat)
    end subroutine
    module subroutine spmv_tridiag_2d_cdp(A, x, y, alpha, beta, op)
        class(tridiagonal_cdp_type), intent(in) :: A
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
        if (op_ == "H") op_ = "C"

        ! Prepare Lapack arguments.
        n = A%n ; ldx = n ; ldy = n ; y = 0.0_dp
        nrhs =  size(x, dim=2, kind=ilp) 

        call lagtm(op_, n, nrhs, alpha_, A%dl, A%dv, A%du, x, ldx, beta_, y, ldy)
    end subroutine

    !-------------------------------------
    !-----                           -----
    !-----     UTILITY FUNCTIONS     -----
    !-----                           -----
    !-------------------------------------

    !----- To dense matrix -----

    pure module function tridiagonal_to_dense_sp(A) result(B)
        !! Convert a `tridiagonal` matrix to its dense representation.
        class(tridiagonal_sp_type), intent(in) :: A
        !! Input tridiagonal matrix.
        real(sp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
        allocate(B(n, n), source=zero_sp)
        B(1, 1) = A%dv(1) ; B(1, 2) = A%du(1)
        do concurrent (i=2:n-1)
            B(i, i-1) = A%dl(i-1)
            B(i, i) = A%dv(i)
            B(i, i+1) = A%du(i)
        enddo
        B(n, n-1) = A%dl(n-1) ; B(n, n) = A%dv(n)
        end associate
    end function
    pure module function tridiagonal_to_dense_dp(A) result(B)
        !! Convert a `tridiagonal` matrix to its dense representation.
        class(tridiagonal_dp_type), intent(in) :: A
        !! Input tridiagonal matrix.
        real(dp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
        allocate(B(n, n), source=zero_dp)
        B(1, 1) = A%dv(1) ; B(1, 2) = A%du(1)
        do concurrent (i=2:n-1)
            B(i, i-1) = A%dl(i-1)
            B(i, i) = A%dv(i)
            B(i, i+1) = A%du(i)
        enddo
        B(n, n-1) = A%dl(n-1) ; B(n, n) = A%dv(n)
        end associate
    end function
    pure module function tridiagonal_to_dense_csp(A) result(B)
        !! Convert a `tridiagonal` matrix to its dense representation.
        class(tridiagonal_csp_type), intent(in) :: A
        !! Input tridiagonal matrix.
        complex(sp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
        allocate(B(n, n), source=zero_csp)
        B(1, 1) = A%dv(1) ; B(1, 2) = A%du(1)
        do concurrent (i=2:n-1)
            B(i, i-1) = A%dl(i-1)
            B(i, i) = A%dv(i)
            B(i, i+1) = A%du(i)
        enddo
        B(n, n-1) = A%dl(n-1) ; B(n, n) = A%dv(n)
        end associate
    end function
    pure module function tridiagonal_to_dense_cdp(A) result(B)
        !! Convert a `tridiagonal` matrix to its dense representation.
        class(tridiagonal_cdp_type), intent(in) :: A
        !! Input tridiagonal matrix.
        complex(dp), allocatable :: B(:, :)
        !! Corresponding dense matrix.

        ! Internal variables.
        integer(ilp) :: i

        associate (n => A%n)
        allocate(B(n, n), source=zero_cdp)
        B(1, 1) = A%dv(1) ; B(1, 2) = A%du(1)
        do concurrent (i=2:n-1)
            B(i, i-1) = A%dl(i-1)
            B(i, i) = A%dv(i)
            B(i, i+1) = A%du(i)
        enddo
        B(n, n-1) = A%dl(n-1) ; B(n, n) = A%dv(n)
        end associate
    end function

    !----- Matrix transposition -----

    pure module function transpose_tridiagonal_sp(A) result(B)
        type(tridiagonal_sp_type), intent(in) :: A
        !! Input matrix.
        type(tridiagonal_sp_type) :: B
        B = tridiagonal(A%du, A%dv, A%dl)
    end function

    pure module function transpose_symtridiagonal_sp(A) result(B)
        type(symtridiagonal_sp_type), intent(in) :: A
        type(symtridiagonal_sp_type) :: B
        B = symtridiagonal(A%dv, A%du)
    end function
        
    pure module function transpose_tridiagonal_dp(A) result(B)
        type(tridiagonal_dp_type), intent(in) :: A
        !! Input matrix.
        type(tridiagonal_dp_type) :: B
        B = tridiagonal(A%du, A%dv, A%dl)
    end function

    pure module function transpose_symtridiagonal_dp(A) result(B)
        type(symtridiagonal_dp_type), intent(in) :: A
        type(symtridiagonal_dp_type) :: B
        B = symtridiagonal(A%dv, A%du)
    end function
        
    pure module function transpose_tridiagonal_csp(A) result(B)
        type(tridiagonal_csp_type), intent(in) :: A
        !! Input matrix.
        type(tridiagonal_csp_type) :: B
        B = tridiagonal(A%du, A%dv, A%dl)
    end function

    pure module function transpose_symtridiagonal_csp(A) result(B)
        type(symtridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type) :: B
        B = symtridiagonal(A%dv, A%du)
    end function
        
    pure module function transpose_hermtridiagonal_csp(A) result(B)
        type(hermtridiagonal_csp_type), intent(in) :: A
        type(hermtridiagonal_csp_type) :: B
        B = hermtridiagonal(A%dv, A%dl)
    end function
    pure module function transpose_tridiagonal_cdp(A) result(B)
        type(tridiagonal_cdp_type), intent(in) :: A
        !! Input matrix.
        type(tridiagonal_cdp_type) :: B
        B = tridiagonal(A%du, A%dv, A%dl)
    end function

    pure module function transpose_symtridiagonal_cdp(A) result(B)
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type) :: B
        B = symtridiagonal(A%dv, A%du)
    end function
        
    pure module function transpose_hermtridiagonal_cdp(A) result(B)
        type(hermtridiagonal_cdp_type), intent(in) :: A
        type(hermtridiagonal_cdp_type) :: B
        B = hermtridiagonal(A%dv, A%dl)
    end function

    !----- Hermitian operator -----

    pure module function hermitian_tridiagonal_sp(A) result(B)
        type(tridiagonal_sp_type), intent(in) :: A
        !! Input matrix.
        type(tridiagonal_sp_type) :: B
        B = tridiagonal(A%du, A%dv, A%dl)
    end function
       
    pure module function hermitian_symtridiagonal_sp(A) result(B)
        type(symtridiagonal_sp_type), intent(in) :: A
        type(symtridiagonal_sp_type) :: B
        B = A
    end function
        
    pure module function hermitian_tridiagonal_dp(A) result(B)
        type(tridiagonal_dp_type), intent(in) :: A
        !! Input matrix.
        type(tridiagonal_dp_type) :: B
        B = tridiagonal(A%du, A%dv, A%dl)
    end function
       
    pure module function hermitian_symtridiagonal_dp(A) result(B)
        type(symtridiagonal_dp_type), intent(in) :: A
        type(symtridiagonal_dp_type) :: B
        B = A
    end function
        
    pure module function hermitian_tridiagonal_csp(A) result(B)
        type(tridiagonal_csp_type), intent(in) :: A
        !! Input matrix.
        type(tridiagonal_csp_type) :: B
        B = tridiagonal(conjg(A%du), conjg(A%dv), conjg(A%dl))
    end function
       
    pure module function hermitian_symtridiagonal_csp(A) result(B)
        type(symtridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type) :: B
        B = symtridiagonal(conjg(A%dv), conjg(A%du))
    end function
        
    pure module function hermitian_hermtridiagonal_csp(A) result(B)
        type(hermtridiagonal_csp_type), intent(in) :: A
        type(hermtridiagonal_csp_type) :: B
        B = A
    end function
    pure module function hermitian_tridiagonal_cdp(A) result(B)
        type(tridiagonal_cdp_type), intent(in) :: A
        !! Input matrix.
        type(tridiagonal_cdp_type) :: B
        B = tridiagonal(conjg(A%du), conjg(A%dv), conjg(A%dl))
    end function
       
    pure module function hermitian_symtridiagonal_cdp(A) result(B)
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type) :: B
        B = symtridiagonal(conjg(A%dv), conjg(A%du))
    end function
        
    pure module function hermitian_hermtridiagonal_cdp(A) result(B)
        type(hermtridiagonal_cdp_type), intent(in) :: A
        type(hermtridiagonal_cdp_type) :: B
        B = A
    end function

    !----- Scalar multiplication -----

    pure module function scalar_multiplication_tridiagonal_sp(alpha, A) result(B)
        real(sp), intent(in) :: alpha
        type(tridiagonal_sp_type), intent(in) :: A
        type(tridiagonal_sp_type) :: B
        B = tridiagonal(A%dl, A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function

    pure module function scalar_multiplication_bis_tridiagonal_sp(A, alpha) result(B)
        type(tridiagonal_sp_type), intent(in) :: A
        real(sp), intent(in) :: alpha
        type(tridiagonal_sp_type) :: B
        B = scalar_multiplication_tridiagonal_sp(alpha, A)
    end function

    pure module function scalar_multiplication_symtridiagonal_sp(alpha, A) result(B)
        real(sp), intent(in) :: alpha
        type(symtridiagonal_sp_type), intent(in) :: A
        type(symtridiagonal_sp_type) :: B
        B = symtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = B%dl
    end function

    pure module function scalar_multiplication_bis_symtridiagonal_sp(A, alpha) result(B)
        type(symtridiagonal_sp_type), intent(in) :: A
        real(sp), intent(in) :: alpha
        type(symtridiagonal_sp_type) :: B
        B = symtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = B%dl
    end function

    pure module function scalar_multiplication_tridiagonal_dp(alpha, A) result(B)
        real(dp), intent(in) :: alpha
        type(tridiagonal_dp_type), intent(in) :: A
        type(tridiagonal_dp_type) :: B
        B = tridiagonal(A%dl, A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function

    pure module function scalar_multiplication_bis_tridiagonal_dp(A, alpha) result(B)
        type(tridiagonal_dp_type), intent(in) :: A
        real(dp), intent(in) :: alpha
        type(tridiagonal_dp_type) :: B
        B = scalar_multiplication_tridiagonal_dp(alpha, A)
    end function

    pure module function scalar_multiplication_symtridiagonal_dp(alpha, A) result(B)
        real(dp), intent(in) :: alpha
        type(symtridiagonal_dp_type), intent(in) :: A
        type(symtridiagonal_dp_type) :: B
        B = symtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = B%dl
    end function

    pure module function scalar_multiplication_bis_symtridiagonal_dp(A, alpha) result(B)
        type(symtridiagonal_dp_type), intent(in) :: A
        real(dp), intent(in) :: alpha
        type(symtridiagonal_dp_type) :: B
        B = symtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = B%dl
    end function

    pure module function scalar_multiplication_tridiagonal_csp(alpha, A) result(B)
        complex(sp), intent(in) :: alpha
        type(tridiagonal_csp_type), intent(in) :: A
        type(tridiagonal_csp_type) :: B
        B = tridiagonal(A%dl, A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function

    pure module function scalar_multiplication_bis_tridiagonal_csp(A, alpha) result(B)
        type(tridiagonal_csp_type), intent(in) :: A
        complex(sp), intent(in) :: alpha
        type(tridiagonal_csp_type) :: B
        B = scalar_multiplication_tridiagonal_csp(alpha, A)
    end function

    pure module function scalar_multiplication_symtridiagonal_csp(alpha, A) result(B)
        complex(sp), intent(in) :: alpha
        type(symtridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type) :: B
        B = symtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = B%dl
    end function

    pure module function scalar_multiplication_bis_symtridiagonal_csp(A, alpha) result(B)
        type(symtridiagonal_csp_type), intent(in) :: A
        complex(sp), intent(in) :: alpha
        type(symtridiagonal_csp_type) :: B
        B = symtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = B%dl
    end function

     pure module function scalar_multiplication_hermtridiagonal_csp(alpha, A) result(B)
        complex(sp), intent(in) :: alpha
        type(hermtridiagonal_csp_type), intent(in) :: A
        type(tridiagonal_csp_type) :: B
        B = tridiagonal(A%dl, A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function

    pure module function scalar_multiplication_bis_hermtridiagonal_csp(A, alpha) result(B)
        type(hermtridiagonal_csp_type), intent(in) :: A
        complex(sp), intent(in) :: alpha
        type(tridiagonal_csp_type) :: B
        B = tridiagonal(A%dl, A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function
    
    pure module function real_scalar_multiplication_hermtridiagonal_csp(alpha, A) result(B)
        real(sp), intent(in) :: alpha
        type(hermtridiagonal_csp_type), intent(in) :: A
        type(hermtridiagonal_csp_type) :: B
        B = hermtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function

    pure module function real_scalar_multiplication_bis_hermtridiagonal_csp(A, alpha) result(B)
        type(hermtridiagonal_csp_type), intent(in) :: A
        real(sp), intent(in) :: alpha
        type(hermtridiagonal_csp_type) :: B
        B = hermtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function
    pure module function scalar_multiplication_tridiagonal_cdp(alpha, A) result(B)
        complex(dp), intent(in) :: alpha
        type(tridiagonal_cdp_type), intent(in) :: A
        type(tridiagonal_cdp_type) :: B
        B = tridiagonal(A%dl, A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function

    pure module function scalar_multiplication_bis_tridiagonal_cdp(A, alpha) result(B)
        type(tridiagonal_cdp_type), intent(in) :: A
        complex(dp), intent(in) :: alpha
        type(tridiagonal_cdp_type) :: B
        B = scalar_multiplication_tridiagonal_cdp(alpha, A)
    end function

    pure module function scalar_multiplication_symtridiagonal_cdp(alpha, A) result(B)
        complex(dp), intent(in) :: alpha
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type) :: B
        B = symtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = B%dl
    end function

    pure module function scalar_multiplication_bis_symtridiagonal_cdp(A, alpha) result(B)
        type(symtridiagonal_cdp_type), intent(in) :: A
        complex(dp), intent(in) :: alpha
        type(symtridiagonal_cdp_type) :: B
        B = symtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = B%dl
    end function

     pure module function scalar_multiplication_hermtridiagonal_cdp(alpha, A) result(B)
        complex(dp), intent(in) :: alpha
        type(hermtridiagonal_cdp_type), intent(in) :: A
        type(tridiagonal_cdp_type) :: B
        B = tridiagonal(A%dl, A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function

    pure module function scalar_multiplication_bis_hermtridiagonal_cdp(A, alpha) result(B)
        type(hermtridiagonal_cdp_type), intent(in) :: A
        complex(dp), intent(in) :: alpha
        type(tridiagonal_cdp_type) :: B
        B = tridiagonal(A%dl, A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function
    
    pure module function real_scalar_multiplication_hermtridiagonal_cdp(alpha, A) result(B)
        real(dp), intent(in) :: alpha
        type(hermtridiagonal_cdp_type), intent(in) :: A
        type(hermtridiagonal_cdp_type) :: B
        B = hermtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function

    pure module function real_scalar_multiplication_bis_hermtridiagonal_cdp(A, alpha) result(B)
        type(hermtridiagonal_cdp_type), intent(in) :: A
        real(dp), intent(in) :: alpha
        type(hermtridiagonal_cdp_type) :: B
        B = hermtridiagonal(A%dv, A%du)
        B%dl = alpha*B%dl; B%dv = alpha*B%dv; B%du = alpha*B%du
    end function

    !----- Matrix addition -----

    ! Tridiag + Tridiag = Tridiag
    pure module function matrix_add_tridiag_tridiag_sp(A, B) result(C)
        type(tridiagonal_sp_type), intent(in) :: A, B
        type(tridiagonal_sp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! Tridiag + SymTridiag = Tridiag
    pure module function matrix_add_tridiag_symtridiag_sp(A, B) result(C)
        type(tridiagonal_sp_type), intent(in) :: A
        type(symtridiagonal_sp_type), intent(in) :: B
        type(tridiagonal_sp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + Tridiag = Tridiag
    pure module function matrix_add_symtridiag_tridiag_sp(A, B) result(C)
        type(symtridiagonal_sp_type), intent(in) :: A
        type(tridiagonal_sp_type), intent(in) :: B
        type(tridiagonal_sp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + SymTridiag = SymTridiag
    pure module function matrix_add_symtridiag_symtridiag_sp(A, B) result(C)
        type(symtridiagonal_sp_type), intent(in) :: A, B
        type(symtridiagonal_sp_type) :: C
        C = symtridiagonal(A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%dl
    end function

    ! Tridiag + Tridiag = Tridiag
    pure module function matrix_add_tridiag_tridiag_dp(A, B) result(C)
        type(tridiagonal_dp_type), intent(in) :: A, B
        type(tridiagonal_dp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! Tridiag + SymTridiag = Tridiag
    pure module function matrix_add_tridiag_symtridiag_dp(A, B) result(C)
        type(tridiagonal_dp_type), intent(in) :: A
        type(symtridiagonal_dp_type), intent(in) :: B
        type(tridiagonal_dp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + Tridiag = Tridiag
    pure module function matrix_add_symtridiag_tridiag_dp(A, B) result(C)
        type(symtridiagonal_dp_type), intent(in) :: A
        type(tridiagonal_dp_type), intent(in) :: B
        type(tridiagonal_dp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + SymTridiag = SymTridiag
    pure module function matrix_add_symtridiag_symtridiag_dp(A, B) result(C)
        type(symtridiagonal_dp_type), intent(in) :: A, B
        type(symtridiagonal_dp_type) :: C
        C = symtridiagonal(A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%dl
    end function

    ! Tridiag + Tridiag = Tridiag
    pure module function matrix_add_tridiag_tridiag_csp(A, B) result(C)
        type(tridiagonal_csp_type), intent(in) :: A, B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! Tridiag + SymTridiag = Tridiag
    pure module function matrix_add_tridiag_symtridiag_csp(A, B) result(C)
        type(tridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + Tridiag = Tridiag
    pure module function matrix_add_symtridiag_tridiag_csp(A, B) result(C)
        type(symtridiagonal_csp_type), intent(in) :: A
        type(tridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + SymTridiag = SymTridiag
    pure module function matrix_add_symtridiag_symtridiag_csp(A, B) result(C)
        type(symtridiagonal_csp_type), intent(in) :: A, B
        type(symtridiagonal_csp_type) :: C
        C = symtridiagonal(A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%dl
    end function

    ! Tridiag + HermTridiag = Tridiag
    pure module function matrix_add_tridiag_hermtridiag_csp(A, B) result(C)
        type(tridiagonal_csp_type), intent(in) :: A
        type(hermtridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! HermTridiag + Tridiag = Tridiag
    pure module function matrix_add_hermtridiag_tridiag_csp(A, B) result(C)
        type(hermtridiagonal_csp_type), intent(in) :: A
        type(tridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + HermTridiag = Tridiag
    pure module function matrix_add_symtridiag_hermtridiag_csp(A, B) result(C)
        type(symtridiagonal_csp_type), intent(in) :: A
        type(hermtridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! HermTridiag + SymTridiag = Tridiag
    pure module function matrix_add_hermtridiag_symtridiag_csp(A, B) result(C)
        type(hermtridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! HermTridiag + HermTridiag = HermTridiag
    pure module function matrix_add_hermtridiag_hermtridiag_csp(A, B) result(C)
        type(hermtridiagonal_csp_type), intent(in) :: A, B
        type(hermtridiagonal_csp_type) :: C
        C = hermtridiagonal(A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = conjg(C%dl)
    end function
    ! Tridiag + Tridiag = Tridiag
    pure module function matrix_add_tridiag_tridiag_cdp(A, B) result(C)
        type(tridiagonal_cdp_type), intent(in) :: A, B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! Tridiag + SymTridiag = Tridiag
    pure module function matrix_add_tridiag_symtridiag_cdp(A, B) result(C)
        type(tridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + Tridiag = Tridiag
    pure module function matrix_add_symtridiag_tridiag_cdp(A, B) result(C)
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(tridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + SymTridiag = SymTridiag
    pure module function matrix_add_symtridiag_symtridiag_cdp(A, B) result(C)
        type(symtridiagonal_cdp_type), intent(in) :: A, B
        type(symtridiagonal_cdp_type) :: C
        C = symtridiagonal(A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%dl
    end function

    ! Tridiag + HermTridiag = Tridiag
    pure module function matrix_add_tridiag_hermtridiag_cdp(A, B) result(C)
        type(tridiagonal_cdp_type), intent(in) :: A
        type(hermtridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! HermTridiag + Tridiag = Tridiag
    pure module function matrix_add_hermtridiag_tridiag_cdp(A, B) result(C)
        type(hermtridiagonal_cdp_type), intent(in) :: A
        type(tridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! SymTridiag + HermTridiag = Tridiag
    pure module function matrix_add_symtridiag_hermtridiag_cdp(A, B) result(C)
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(hermtridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! HermTridiag + SymTridiag = Tridiag
    pure module function matrix_add_hermtridiag_symtridiag_cdp(A, B) result(C)
        type(hermtridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = C%du + B%du
    end function

    ! HermTridiag + HermTridiag = HermTridiag
    pure module function matrix_add_hermtridiag_hermtridiag_cdp(A, B) result(C)
        type(hermtridiagonal_cdp_type), intent(in) :: A, B
        type(hermtridiagonal_cdp_type) :: C
        C = hermtridiagonal(A%dv, A%du)
        C%dl = C%dl + B%dl ; C%dv = C%dv + B%dv ; C%du = conjg(C%dl)
    end function

    !----- Matrix subtraction -----

    ! Tridiag - Tridiag = Tridiag
    pure module function matrix_sub_tridiag_tridiag_sp(A, B) result(C)
        type(tridiagonal_sp_type), intent(in) :: A, B
        type(tridiagonal_sp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! Tridiag - SymTridiag = Tridiag
    pure module function matrix_sub_tridiag_symtridiag_sp(A, B) result(C)
        type(tridiagonal_sp_type), intent(in) :: A
        type(symtridiagonal_sp_type), intent(in) :: B
        type(tridiagonal_sp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - Tridiag = Tridiag
    pure module function matrix_sub_symtridiag_tridiag_sp(A, B) result(C)
        type(symtridiagonal_sp_type), intent(in) :: A
        type(tridiagonal_sp_type), intent(in) :: B
        type(tridiagonal_sp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - SymTridiag = SymTridiag
    pure module function matrix_sub_symtridiag_symtridiag_sp(A, B) result(C)
        type(symtridiagonal_sp_type), intent(in) :: A, B
        type(symtridiagonal_sp_type) :: C
        C = symtridiagonal(A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%dl
    end function

    ! Tridiag - Tridiag = Tridiag
    pure module function matrix_sub_tridiag_tridiag_dp(A, B) result(C)
        type(tridiagonal_dp_type), intent(in) :: A, B
        type(tridiagonal_dp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! Tridiag - SymTridiag = Tridiag
    pure module function matrix_sub_tridiag_symtridiag_dp(A, B) result(C)
        type(tridiagonal_dp_type), intent(in) :: A
        type(symtridiagonal_dp_type), intent(in) :: B
        type(tridiagonal_dp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - Tridiag = Tridiag
    pure module function matrix_sub_symtridiag_tridiag_dp(A, B) result(C)
        type(symtridiagonal_dp_type), intent(in) :: A
        type(tridiagonal_dp_type), intent(in) :: B
        type(tridiagonal_dp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - SymTridiag = SymTridiag
    pure module function matrix_sub_symtridiag_symtridiag_dp(A, B) result(C)
        type(symtridiagonal_dp_type), intent(in) :: A, B
        type(symtridiagonal_dp_type) :: C
        C = symtridiagonal(A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%dl
    end function

    ! Tridiag - Tridiag = Tridiag
    pure module function matrix_sub_tridiag_tridiag_csp(A, B) result(C)
        type(tridiagonal_csp_type), intent(in) :: A, B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! Tridiag - SymTridiag = Tridiag
    pure module function matrix_sub_tridiag_symtridiag_csp(A, B) result(C)
        type(tridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - Tridiag = Tridiag
    pure module function matrix_sub_symtridiag_tridiag_csp(A, B) result(C)
        type(symtridiagonal_csp_type), intent(in) :: A
        type(tridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - SymTridiag = SymTridiag
    pure module function matrix_sub_symtridiag_symtridiag_csp(A, B) result(C)
        type(symtridiagonal_csp_type), intent(in) :: A, B
        type(symtridiagonal_csp_type) :: C
        C = symtridiagonal(A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%dl
    end function

    ! Tridiag - HermTridiag = Tridiag
    pure module function matrix_sub_tridiag_hermtridiag_csp(A, B) result(C)
        type(tridiagonal_csp_type), intent(in) :: A
        type(hermtridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! HermTridiag - Tridiag = Tridiag
    pure module function matrix_sub_hermtridiag_tridiag_csp(A, B) result(C)
        type(hermtridiagonal_csp_type), intent(in) :: A
        type(tridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - HermTridiag = Tridiag
    pure module function matrix_sub_symtridiag_hermtridiag_csp(A, B) result(C)
        type(symtridiagonal_csp_type), intent(in) :: A
        type(hermtridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! HermTridiag - SymTridiag = Tridiag
    pure module function matrix_sub_hermtridiag_symtridiag_csp(A, B) result(C)
        type(hermtridiagonal_csp_type), intent(in) :: A
        type(symtridiagonal_csp_type), intent(in) :: B
        type(tridiagonal_csp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! HermTridiag - HermTridiag = HermTridiag
    pure module function matrix_sub_hermtridiag_hermtridiag_csp(A, B) result(C)
        type(hermtridiagonal_csp_type), intent(in) :: A, B
        type(hermtridiagonal_csp_type) :: C
        C = hermtridiagonal(A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = conjg(C%dl)
    end function
    ! Tridiag - Tridiag = Tridiag
    pure module function matrix_sub_tridiag_tridiag_cdp(A, B) result(C)
        type(tridiagonal_cdp_type), intent(in) :: A, B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! Tridiag - SymTridiag = Tridiag
    pure module function matrix_sub_tridiag_symtridiag_cdp(A, B) result(C)
        type(tridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - Tridiag = Tridiag
    pure module function matrix_sub_symtridiag_tridiag_cdp(A, B) result(C)
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(tridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - SymTridiag = SymTridiag
    pure module function matrix_sub_symtridiag_symtridiag_cdp(A, B) result(C)
        type(symtridiagonal_cdp_type), intent(in) :: A, B
        type(symtridiagonal_cdp_type) :: C
        C = symtridiagonal(A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%dl
    end function

    ! Tridiag - HermTridiag = Tridiag
    pure module function matrix_sub_tridiag_hermtridiag_cdp(A, B) result(C)
        type(tridiagonal_cdp_type), intent(in) :: A
        type(hermtridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! HermTridiag - Tridiag = Tridiag
    pure module function matrix_sub_hermtridiag_tridiag_cdp(A, B) result(C)
        type(hermtridiagonal_cdp_type), intent(in) :: A
        type(tridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! SymTridiag - HermTridiag = Tridiag
    pure module function matrix_sub_symtridiag_hermtridiag_cdp(A, B) result(C)
        type(symtridiagonal_cdp_type), intent(in) :: A
        type(hermtridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! HermTridiag - SymTridiag = Tridiag
    pure module function matrix_sub_hermtridiag_symtridiag_cdp(A, B) result(C)
        type(hermtridiagonal_cdp_type), intent(in) :: A
        type(symtridiagonal_cdp_type), intent(in) :: B
        type(tridiagonal_cdp_type) :: C
        C = tridiagonal(A%dl, A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = C%du - B%du
    end function

    ! HermTridiag - HermTridiag = HermTridiag
    pure module function matrix_sub_hermtridiag_hermtridiag_cdp(A, B) result(C)
        type(hermtridiagonal_cdp_type), intent(in) :: A, B
        type(hermtridiagonal_cdp_type) :: C
        C = hermtridiagonal(A%dv, A%du)
        C%dl = C%dl - B%dl ; C%dv = C%dv - B%dv ; C%du = conjg(C%dl)
    end function

end submodule
