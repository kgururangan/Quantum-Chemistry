function [Nelec, NFZ_core, NFZ_vir, Nact_h, Nact_p, run] = check_run_config(Nelec,NFZ_core,NFZ_vir,Nact_h,Nact_p,run)



    assert(2*NFZ_core < Nelec, 'Error: 2 * Number of frozen occupied orbitals should be less than number of electrons!')

    if isempty(run.flag_RHF)
        disp('Warning: Restriced HF flag not specified! Assuming closed-shell RHF...')
        run.flag_RHF = 1;
    end

    if isempty(run.diis_size)
        disp('Warning: DIIS size not specified. Using default diis_size = 3...')
        run.diis_size = 3;
    end

    if isempty(run.cc_tol)
        disp('Warning: CC convergence tolerance not specified. Using default 1e-8...')
        run.cc_tol = 1e-8;
    end

    if isempty(run.cc_maxit)
        disp('Warning: Maximum number of iterations for CC calculation not specified! Using default maxit = 100 ...')
        run.cc_maxit = 100;
    end

    if isempty(run.eom_maxit)
        disp('Warning: Maximum number of iterations for EOMCC calculation not specified! Using default maxit = 1000 ...')
        run.eom_maxit = 1000;
    end

    if run.left_cc == 1
        if isempty(run.lcc_maxit)
            fprintf('Warning: Maximum number of iterations for EOMCC calculation not specified! Using default maxit = %d ...\n',50*run.cc_maxit)
            run.lcc_maxit = 50*run.cc_maxit;
        end
    end

    if strcmp(run.calc_type, 'CCSDT3')

        if isempty(Nact_h) || isempty(Nact_p)
            disp('Warning: Number of active hole or particle orbitals not specified! Full space assumed active...')
            run.Nact_h = 1e4;
            run.Nact_p = 1e4;
        end

    end

    if strcmp(run.calc_type,'EOMCCSD')

        assert(~isempty(run.nroot),'Error: Number of states not specified for EOM calculation!')

        if isempty(run.eom_tol)
            disp('Warning: EOM convergence tolerance not specified. Using default 1e-6...')
            run.eom_tol = 1e-6;
        end

        if isempty(run.max_nvec_per_root)
            disp('Warning: EOM Maximum number of vectors per root not specified. Using default max_nvec_per_root = 4...')
            run.max_nvec_per_root = 4;
        end

	if isempty(run.nvec_per_root)
            disp('Warning: EOM number of vectors per root not specified. Using default nvec_per_root = 1...')
            run.nvec_per_root = 1;
        end
    end

    if ~strcmp(run.calc_type,'CCSDT3')
        Nact_h = 10000;
        Nact_p = 10000;
    end

end

