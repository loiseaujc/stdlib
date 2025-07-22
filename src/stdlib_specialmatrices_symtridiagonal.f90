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

end submodule
