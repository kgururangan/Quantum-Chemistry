function [R, omega, res, cc_t] = eomuccsd(HBar_t,cc_t,sys,opts,B0)

    fprintf('\n==================================++Entering EOM-UCCSD Routine++=============================\n')
    
    nroot = opts.nroot;
    mult = opts.mult;

    if nargin < 4
        B0 = [];
    end

    % Matrix-vector product function
    HRmat = @(x) HR_ucc_matmat(x,HBar_t,cc_t,sys);

    if strcmp(opts.init_guess,'cis')
        % use CIS guess as initial basis vectors
        %[omega, C1] = cis_spinadapt(nroot,sys,mult);
        [cis_energy, C1A, C1B] = cis_spinint(floor(sys.singles_dim/2),sys,mult);

        B1A = zeros(size(C1A,1),nroot);
        B1B = zeros(size(C1B,1),nroot);
        omega_cis = zeros(1,nroot);

        ct = 1;

        if mult == 1
           for i = 1:length(cis_energy)
               chk = norm(C1A(:,i)-C1B(:,i));
               if abs(chk) < 1e-9
                   B1A(:,ct) = C1A(:,i);
                   B1B(:,ct) = C1B(:,i);
                   omega_cis(ct) = cis_energy(i);
                   if ct == nroot
                        break
                   end 
                  ct = ct + 1;
               end
            end
         end

         if mult == 3
             for i = 1:length(cis_energy)
                 chk = norm(C1A(:,i)-C1B(:,i))/norm(C1A(:,i));
                 if abs(chk-2) < 1e-6
                     B1A(:,ct) = C1A(:,i);
                     B1B(:,ct) = C1B(:,i);
                     omega_cis(ct) = cis_energy(i);
                     if ct == nroot
                         break
                     end
                     ct = ct + 1;
                 end
              end
          end

    %     if mult == 1
    %         % closed-shell singlet spin-adaption sets r1a = r1b
    %         B1A = C1A;
    %         B1B = C1B;
    %     end
    %     if mult == 3
    %         % closed-shell triplet spin-adaption sets r1a = -r1b
    %         B1A = C1A;
    %         B1B = C1B;
    %     end


        B2A = zeros(length(sys.posv{3}),size(B1A,2));
        B2B = zeros(length(sys.posv{4}),size(B1A,2));
        B2C = zeros(length(sys.posv{5}),size(B1A,2));

        B0 = cat(1,B1A,B1B,B2A,B2B,B2C);

        fprintf('Initial CIS energies (S = %d):\n',0.5*(mult-1))
        for i = 1:length(omega_cis)
            fprintf('      E%d = %4.8f\n',i,omega_cis(i))
        end
        fprintf('\n')
    end

    if opts.solver == 2
        [R, omega, res, ~] = davidson_solver2(HRmat, HBar_t, cc_t, nroot, B0, sys, opts);
    elseif opts.solver == 1
        [R, omega, res, ~] = davidson_solver1(HRmat, HBar_t, cc_t, nroot, B0, sys, opts);
    else
        [R, omega, res, ~] = davidson_solver3(HRmat, HBar_t, cc_t, nroot, B0, sys, opts);
    end
    omega = real(omega);
    
    % calculate r0 for each excited state
    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;

    for i = 1:nroot
        r1a = reshape(R(sys.posv{1},i),sys.size{1});
        r1b = reshape(R(sys.posv{2},i),sys.size{2});
        r2a = reshape(R(sys.posv{3},i),sys.size{3});
        r2b = reshape(R(sys.posv{4},i),sys.size{4});
        r2c = reshape(R(sys.posv{5},i),sys.size{5});

        HR0 = einsum_kg(H1A.ov,r1a,'me,em->')...
               + einsum_kg(H1B.ov,r1b,'me,em->')...
               + 0.25*einsum_kg(H2A.oovv,r2a,'mnef,efmn->')...
               + einsum_kg(H2B.oovv,r2b,'mnef,efmn->')...
               + 0.25*einsum_kg(H2C.oovv,r2c,'mnef,efmn->');

        r0 = HR0/omega(i);
        
        cc_t.r1a{i} = r1a;
        cc_t.r1b{i} = r1b;
        cc_t.r2a{i} = r2a;
        cc_t.r2b{i} = r2b;
        cc_t.r2c{i} = r2c;
        cc_t.r0(i) = real(r0);
    end


end