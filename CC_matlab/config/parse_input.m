function [Nelec, NFZ_core, NFZ_vir, Nact_h, Nact_p, onebody_path, twobody_path, run] = parse_input(input_file)

    Nelec = []; NFZ_core = []; NFZ_vir = []; Nact_h = []; Nact_p = [];
    onebody_path = []; twobody_path = []; 

    run.nroot = []; run.calc_type = []; run.cc_tol = [];
    run.flag_RHF = []; run.diis_size = []; run.eom_tol = [];
    run.nvec_per_root = []; run.left_cc = [];
    run.cc_maxit = [];  run.eom_maxit = []; run.lcc_maxit = [];

    fid = fopen(input_file);


    tline = fgetl(fid);
    while ischar(tline)

%         if contains(tline, '%')
%             continue
%         end

        LINE = strip(split(tline,'='));

        inp = LINE{1};

        if strcmp(inp,'n') || strcmp(inp,'nelec')
            LL = strip(split(LINE{2})); LL = LL{1};
            Nelec = str2double(LL);
        elseif strcmp(inp,'nfzc')
            LL = strip(split(LINE{2})); LL = LL{1};
            NFZ_core = str2double(LL);
        elseif strcmp(inp,'nfzv')
            LL = strip(split(LINE{2})); LL = LL{1};
            NFZ_vir = str2double(LL);
        elseif strcmp(inp,'nact_h')
            LL = strip(split(LINE{2})); LL = LL{1};
            Nact_h = str2double(LL);
        elseif strcmp(inp,'nact_p')
            LL = strip(split(LINE{2})); LL = LL{1};
            Nact_p = str2double(LL);
        elseif strcmp(inp,'onebody')
            LL = strip(split(LINE{2})); LL = LL{1};
            onebody_path = LL;
        elseif strcmp(inp,'twobody')
            LL = strip(split(LINE{2})); LL = LL{1};
            twobody_path = LL;
        elseif strcmp(inp,'rhf') || strcmp(inp,'RHF') || strcmp(inp,'flag_RHF')
            LL = strip(split(LINE{2})); ifRHF = LL{1};
            if strcmp(ifRHF,'False') || strcmp(ifRHF,'false') || strcmp(ifRHF,'F') || strcmp(ifRHF,'f')
                run.flag_RHF = 0;
            else
                run.flag_RHF = 1;
            end
        elseif strcmp(inp,'diis_size') || strcmp(inp,'DIIS_size') 
            LL = strip(split(LINE{2})); LL = LL{1};
            run.diis_size = str2double(LL);
	elseif strcmp(inp,'eom_nvec_per_root') || strcmp(inp,'nvec_per_root')
        LL = strip(split(LINE{2})); LL = LL{1};
	    run.nvec_per_root = str2double(LL);
        elseif strcmp(inp,'eom_max_nvec') || strcmp(inp,'max_nvec_per_root')
            LL = strip(split(LINE{2})); LL = LL{1};
            run.max_nvec_per_root = str2double(LL);
        elseif strcmp(inp,'cc_conv_tol') || strcmp(inp,'cc_tol')
            LL = strip(split(LINE{2})); LL = LL{1};
            run.cc_tol = 10^(-str2double(LL));
        elseif strcmp(inp,'eom_conv_tol') || strcmp(inp,'eom_tol')
            LL = strip(split(LINE{2})); LL = LL{1};
            run.eom_tol = 10^(-str2double(LL));
        elseif strcmp(inp,'nroot')
            LL = strip(split(LINE{2})); LL = LL{1};
            run.nroot = str2double(LL);
        elseif strcmp(inp,'do_left')  || strcmp(inp,'left') || strcmp(inp,'DO_LEFT') || strcmp(inp,'LEFT')
            LL = strip(split(LINE{2})); LL = LL{1};
            ifLCC = LL;
            if strcmp(ifLCC,'True') || strcmp(ifLCC,'true') || strcmp(ifLCC,'T') || strcmp(ifLCC,'t')
                run.left_cc = 1;
            else
                run.left_cc = 0;
            end
        elseif strcmp(inp,'cc_maxit') 
            LL = strip(split(LINE{2})); LL = LL{1};
            run.cc_maxit = str2double(LL);
        elseif strcmp(inp,'eom_maxit')
            LL = strip(split(LINE{2})); LL = LL{1};
            run.eom_maxit = str2double(LL);
        elseif strcmp(inp,'lcc_maxit')  
            LL = strip(split(LINE{2})); LL = LL{1};
            run.lcc_maxit = str2double(LL);
        elseif strcmp(inp,'calc_type')
            LL = strip(split(LINE{2})); LL = LL{1};
            calc = LL;
            if strcmp(calc,'CCSD') || strcmp(calc,'ccsd')
                run.calc_type = 'CCSD';
            elseif strcmp(calc,'UCCSD') || strcmp(calc,'uccsd')
                run.calc_type = 'UCCSD';
            elseif strcmp(calc,'CRCC23') || strcmp(calc,'crcc23')
                run.calc_type = 'CRCC23';
            elseif strcmp(calc,'CCSDT') || strcmp(calc,'ccsdt')
                run.calc_type = 'CCSDT';
            elseif strcmp(calc,'CCSDT3') || strcmp(calc,'ccsdt3') || strcmp(calc,'CCSDT-3') || strcmp(calc,'ccsdt-3')
                run.calc_type = 'CCSDT3';
            elseif strcmp(calc,'EOMCCSD') || strcmp(calc,'eomccsd') || strcmp(calc,'eom-ccsd') || strcmp(calc,'EOM-CCSD')
                run.calc_type = 'EOMCCSD';
            else
                disp('Calculation type not recognized!');
            end
        end
        tline = fgetl(fid);
    end

    fclose(fid);
        
    [Nelec, NFZ_core, NFZ_vir, Nact_h, Nact_p, run] = ...
check_run_config(Nelec,NFZ_core,NFZ_vir,Nact_h,Nact_p,run);
    

end

