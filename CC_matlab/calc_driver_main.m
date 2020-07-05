function [] = calc_driver_main(input_file)

    format long

    % Path to source .m files
    SOURCE_DIR = '/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CC_matlab';

    % Add source .m files to path
    addpath(genpath(SOURCE_DIR));
        
    fprintf('==========================================================================================================================\n')
    fprintf('================================================== Coupled-Cluster MATLAB ================================================\n')
    fprintf('==========================================================================================================================\n')
    fprintf('\nRepository online: <https://github.com/kgururangan/Quantum-Chemistry.git>\n')
    fprintf('Author: Karthik Gururangan (gururang@msu.edu)\n\n\n')

    % Parse input to get system information and integral locations
    [Nelec, nfzc, nfzv, nact_h, nact_p, onebody_path, twobody_path, run] = parse_input(input_file);
    % Load in onebody and twobody integrals
    fprintf('Onebody file: %s\n',onebody_path)
    fprintf('Twobody file: %s\n',twobody_path)
    fprintf('Loading integrals...')
    [e1int, e2int, Vnuc, Norb] = load_integrals(onebody_path, twobody_path);
    fprintf(' completed!\n')

    %
    if run.flag_RHF == 1
        Nocc_a = (Nelec/2);
        Nocc_b = (Nelec/2);
    end

    hf_energy = calculate_hf_energy(e1int,e2int,Nocc_a,Nocc_b);

    filename = strip(split(input_file,'.')); proj_name = filename{1};
	
    fprintf('\nProject Name: %s\n',proj_name)
    fprintf('---------------------------------\n')	
    fprintf('Working directory: %s\n\n',pwd')

    fprintf('System Information:\n')
    fprintf('---------------------\n')
    fprintf('   Number of electrons = %d\n',Nelec)
    fprintf('   Number of spatial orbitals = %d\n',Norb)
    fprintf('   Number of frozen core spatial orbitals = %d\n',nfzc)
    fprintf('   Number of frozen virtual spatial orbitals = %d\n',nfzv)
    fprintf('   Number of occupied orbitals (alpha) = %d\n',Nocc_a-nfzc)
    fprintf('   Number of occupied orbitals (beta) = %d\n',Nocc_b-nfzc)
    fprintf('   Number of unoccupied orbitals (alpha) = %d\n',Norb-Nocc_a-nfzv)
    fprintf('   Number of unoccupied orbitals (beta) = %d\n',Norb-Nocc_b-nfzv)
    fprintf('   Number of active occupied spatial orbitals = %d\n',nact_h)
    fprintf('   Number of active unoccupied spatial orbitals = %d\n',nact_p)

    if run.left_cc == 1
	    fprintf('\nCalculation Type = %s + Left-CC\n',run.calc_type)
    else
	    fprintf('\nCalculation Type = %s\n',run.calc_type)
    end

    sys_cc = build_system_cc(e1int,e2int,Vnuc,Nocc_a,Nocc_b,nfzc,nfzv,nact_h,nact_p);
    sys_ucc = build_system_ucc(e1int,e2int,Nocc_a,Nocc_b);

    fprintf('\n\nMolecular Orbital Basis\n')
    fprintf('------------------------\n')
    fprintf('	Nuclear Repulsion Energy = %4.12f Ha\n',Vnuc)
    fprintf('	Reference Energy = %4.12f Ha\n\n',hf_energy+Vnuc)
    fprintf('	MO Index     Energy (Ha)      Occupied\n')
    fprintf('	----------------------------------------\n')
    for i = 1:2:sys_cc.Norb
	    if i <= sys_cc.Nelec
	    	fprintf('	%d	    %4.8f         %s\n', (i+1)/2,sys_cc.FM(i,i),'*');
	    else
		fprintf('	%d	    %4.8f         %s\n', (i+1)/2,sys_cc.FM(i,i),'');
	    end
    end

    t = datetime;
    d = datetime(t,'Format','eeee, MMMM d, yyyy h:mm a');
    fprintf('\nCalculation started on %s (%s) \n',d, t);

    switch run.calc_type

        case 'CCSD'

            ccopts.tol = run.cc_tol;
            ccopts.diis_size = run.diis_size;
            ccopts.maxit = run.cc_maxit;
            [t1,t2,Ecorr_ccsd] = ccsd(sys_cc,ccopts);

            if run.left_cc == 1

                [HBar] = build_HBar_debug(t1,t2,sys_cc);

                lccopts.diis_size = run.diis_size;
                lccopts.maxit = run.lcc_maxit;
                lccopts.tol = run.cc_tol;
                [lambda1,lambda2,lcc_resid] = lccsd(t1,t2,HBar,sys_cc,lccopts);
            end

        case 'UCCSD'

            ccopts.tol = run.cc_tol;
            ccopts.diis_size = run.diis_size;
            ccopts.maxit = run.cc_maxit;
            [t1a,t1b,t2a,t2b,t2c,Ecorr_ucc] = uccsd(sys_ucc,ccopts);

        case 'CCSDT'

            ccopts.tol = run.cc_tol;
            ccopts.diis_size = run.diis_size;
            ccopts.maxit = run.cc_maxit;
            [t1,t2,t3,Ecorr_ccsdt] = ccsdt(sys_cc,ccopts);

        case 'CCSDT3'

            ccopts.tol = run.cc_tol;
            ccopts.diis_size = run.diis_size;
            ccopts.maxit = run.cc_maxit;
            [t1,t2,t3,Ecorr_ccsdt3] = ccsdt3(sys_cc,ccopts);

        case 'CRCC23'
    
            ccopts.tol = run.cc_tol;
            ccopts.diis_size = run.diis_size;
            ccopts.maxit = run.cc_maxit;
            [t1,t2,Ecorr_ccsd] = ccsd(sys_cc,ccopts);

            [HBar] = build_HBar_debug(t1,t2,sys_cc);

            lccopts.diis_size = run.diis_size;
            lccopts.maxit = run.lcc_maxit;
            lccopts.tol = run.cc_tol;
            [lambda1,lambda2,lcc_resid] = lccsd(t1,t2,HBar,sys_cc,lccopts);

            [Ecorr_crcc23A,Ecorr_crcc23B,Ecorr_crcc23C,Ecorr_crcc23D] = ...
                                        crcc23_opt(t1,t2,lambda1,lambda2,HBar,sys_cc);

        case 'EOMCCSD'

            ccopts.tol = run.cc_tol;
            ccopts.diis_size = run.diis_size;
            ccopts.maxit = run.cc_maxit;
            [t1,t2,Ecorr_ccsd] = ccsd(sys_cc,ccopts);

            [HBar] = build_HBar_debug(t1,t2,sys_cc);

            eomopts.nroot = run.nroot;
            eomopts.maxit = run.eom_maxit;
            eomopts.tol = run.eom_tol;
            eomopts.nvec_per_root = run.nvec_per_root;
            eomopts.max_nvec_per_root = run.max_nvec_per_root;
            eomopts.flag_verbose = 1;
            eomopts.init_guess = 'cis';
            eomopts.thresh_vec = 10*run.eom_tol;
	    
            [Rvec, omega, r0, eom_residual] = eomccsd(HBar,t1,t2,sys_cc,eomopts);

            if run.left_cc == 1

                lccopts.diis_size = run.diis_size;
                lccopts.maxit = run.lcc_maxit;
                lccopts.tol = run.cc_tol;
                [lambda1,lambda2,lcc_resid] = lccsd(t1,t2,HBar,sys_cc,lccopts);

                lccopts.maxit = run.eom_maxit;
                lccopts.diis_size = run.diis_size;
                lccopts.tol = run.eom_tol; 

                [Lvec,omega_lcc,eom_lcc_resid] = lefteomccsd(omega,Rvec,HBar,t1,t2,sys_cc,lccopts);

            end

        otherwise
                disp('Calculation type not defined!')
        end

   
	t = datetime;
	d = datetime(t,'Format','eeee, MMMM d, yyyy h:mm a');
	fprintf('\nCalculation ended on %s (%s) \n',d, t);

end

