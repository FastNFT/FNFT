    UINT K = D;
    UINT M = D;
    REAL main_spec[K];
    REAL aux_spec[M];

    fnft_kdvp_opts_t opts = fnft_kdvp_default_opts();
    opts.grid_spacing = 0.001;
    INT ret_code = fnft_kdvp(D, q, X, E, &K, main_spec, &M, aux_spec, NULL, &opts);

    if (K != 3) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    if (M != 1) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    main_spec[1] = main_spec[2]; // remove the signs
    main_spec[2] = main_spec[4];
    REAL err = misc_rel_err_real(K, main_spec, main_spec_exact);
#ifdef DEBUG
    printf("err = %g, err_bnds[0] = %g\n", err, err_bnds[0]);
#endif
    if (!(err <= err_bnds[0])) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    err = misc_rel_err_real(M, aux_spec, aux_spec_exact);
#ifdef DEBUG
    printf("err = %g, err_bnds[1] = %g\n", err, err_bnds[1]);
#endif
    if (!(err <= err_bnds[1])) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    return ret_code;

