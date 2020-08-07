function [omega,c1] = cis_spinadapt(nroot,sys,mult)

    Nocc = sys.Nocc_alpha;
    Nunocc = sys.Nvir_alpha;
    
    if Nocc ~= sys.Nocc_beta || Nunocc ~= sys.Nvir_beta
        disp('Error - Singlet spin-adapted CIS routine requires RHF reference!')
    end

    H_cis = zeros(Nunocc*Nocc,Nunocc*Nocc);

    % closed-shell singlet spin adaptation
    if mult == 1
        ct1 = 1;
        for a = 1:Nunocc
            for i = 1:Nocc
                ct2 = 1;
                for b = 1:Nunocc
                    for j = 1:Nocc
                        tmp = sys.fa_vv(a,b)*deltafcn(i,j)...
                              -sys.fa_oo(j,i)*deltafcn(a,b)...
                              +2*sys.e2int_voov(a,j,i,b)-sys.e2int_vovo(a,j,b,i);
                        H_cis(ct1,ct2) = tmp;
                        ct2 = ct2 + 1;
                    end
                end
                ct1 = ct1 + 1;
            end
        end
    end

    % closed-shell triplet spin adaptation
    if mult == 3
        ct1 = 1;
        for a = 1:Nunocc
            for i = 1:Nocc
                ct2 = 1;
                for b = 1:Nunocc
                    for j = 1:Nocc
                        tmp = sys.fa_vv(a,b)*deltafcn(i,j)...
                              -sys.fa_oo(j,i)*deltafcn(a,b)...
                              -sys.e2int_vovo(a,j,b,i);
                        H_cis(ct1,ct2) = tmp;
                        ct2 = ct2 + 1;
                    end
                end
                ct1 = ct1 + 1;
            end
        end
    end

    [c1,omega] = eig(H_cis); omega = diag(omega);
    [~,idx] = sort(omega,'ascend');
    c1 = c1(:,idx(1:nroot));
    omega = omega(idx(1:nroot));
    
end

function [xout] = deltafcn(a,b)
    if a == b
        xout = 1;
    else
        xout = 0;
    end
end