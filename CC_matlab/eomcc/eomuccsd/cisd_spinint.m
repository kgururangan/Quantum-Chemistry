function [omega,c1a,c1b] = cisd_spinint(nroot,sys,mult,nact_h,nact_p)

    Nocc_a = sys.Nocc_alpha; Nocc_b = sys.Nocc_beta;
    Nunocc_a = sys.Nvir_alpha; Nunocc_b = sys.Nvir_beta;

    H_1A_1A = zeros(Nunocc_a*Nocc_a,Nunocc_a*Nocc_a);
    ct1 = 1;
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            ct2 = 1;
            for j = 1:Nocc_a
                for b = 1:Nunocc_a
                    H_1A_1A(ct1,ct2) = sys.fa_vv(a,b)*(i==j)...
                                  -sys.fa_oo(j,i)*(a==b)...
                                  +sys.vA_voov(a,j,i,b);
                    ct2=ct2+1;
                end
            end
            ct1=ct1+1;
        end
    end

    H_1A_1B = zeros(Nunocc_a*Nocc_a,Nunocc_b*Nocc_b);
    ct1 = 1;
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            ct2 = 1;
            for j = 1:Nocc_b
                for b = 1:Nunocc_b
                    H_1A_1B(ct1,ct2) = sys.vB_voov(a,j,i,b);
                    ct2=ct2+1;
                end
            end
            ct1=ct1+1;
        end
    end

    H_1B_1B = zeros(Nunocc_b*Nocc_b,Nunocc_b*Nocc_b);
    ct1 = 1;
    for i = 1:Nocc_b
        for a = 1:Nunocc_b
            ct2 = 1;
            for j = 1:Nocc_b
                for b = 1:Nunocc_b
                    H_1B_1B(ct1,ct2) = sys.fb_vv(a,b)*(i==j)...
                                  -sys.fb_oo(j,i)*(a==b)...
                                  +sys.vC_voov(a,j,i,b);
                    ct2=ct2+1;
                end
            end
            ct1=ct1+1;
        end
    end

%     H_cis = [H1, H2;
%              H2, H4];
% 
%     [c1,omega] = eig(H_cis); omega = diag(omega);
%     [~,idx] = sort(omega,'ascend');
%     c1 = c1(:,idx(1:nroot));
%     omega = omega(idx(1:nroot));
% 
%     c1a = c1(1:Nunocc_a*Nocc_a,:);
%     c1b = c1(Nunocc_a*Nocc_a+1:end,:);