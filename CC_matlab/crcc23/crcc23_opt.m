function [Ecorr_crcc23A,Ecorr_crcc23B,Ecorr_crcc23C,Ecorr_crcc23D] = crcc23_opt(t1,t2,L1,L2,HBar,sys)

    fprintf('\n==================================++Entering CR-CC(2,3) Routine++=============================\n')

    fprintf('\n')

    fprintf('MM Correction calculation...')
    tic

    Nocc = sys.Nocc;
    Nunocc = sys.Nunocc;

    L1_ia = permute(L1,[2,1]);
    L2_ijab = permute(L2,[3,4,1,2]);
    I_abie = HBar{2}{2,2,1,2}-einsum_kg(HBar{1}{1,2},t2,'me,abim->abie');
    I_bcek = -permute(I_abie,[1,2,4,3]);
    h_amij = HBar{2}{2,1,1,1};

    deltaA = 0.0; % using MP denominator -(f_aa - f_ii + f_bb - f_jj + f_cc - f_kk)
    deltaB = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) | phi_{ijk}^{abc}>
    deltaC = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) | phi_{ijk}^{abc}>
    deltaD = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) + H_3(CCSD) | phi_{ijk}^{abc}>

    for i = 1:Nocc
        for j = i+1:Nocc
            for k = j+1:Nocc

                % calculate < 0 | L H(CCSD) | phi_{ijkabc} >
                temp_ijk = build_partial_L3(i,j,k,L1_ia,L2_ijab,HBar);
                temp_kji = build_partial_L3(k,j,i,L1_ia,L2_ijab,HBar);
                temp_ikj = build_partial_L3(i,k,j,L1_ia,L2_ijab,HBar);
                NUM = temp_ijk - temp_kji - temp_ikj;

                % calculate < phi_{ijkabc} | H(CCSD) | 0 >
                temp_ijk = build_partial_MM23(i,j,k,t1,t2,h_amij,I_bcek);
                temp_kji = build_partial_MM23(k,j,i,t1,t2,h_amij,I_bcek);
                temp_ikj = build_partial_MM23(i,k,j,t1,t2,h_amij,I_bcek);
                BETA = temp_ijk - temp_kji - temp_ikj;

                for a = 1:Nunocc
                    for b = a+1:Nunocc
                        for c = b+1:Nunocc

                            MM23 = BETA(a,b,c) - BETA(b,a,c) - BETA(c,b,a);

                            L3 = NUM(a,b,c) - NUM(b,a,c) - NUM(c,b,a) ...
                                   - NUM(a,c,b) + NUM(b,c,a) + NUM(c,a,b);

                            TMP = MM23*L3;

                            DMP = sys.fvv(a,a) + sys.fvv(b,b) + sys.fvv(c,c) ...
                                  -sys.foo(i,i) - sys.foo(j,j) - sys.foo(k,k);

                            [D1,D2,D3] = HBar_CCSD_T_diagonal(a,b,c,i,j,k,HBar);

                            D_A = -DMP;
                            D_B = -D1;
                            D_C = -(D1+D2);
                            D_D = -(D1+D2+D3);

                            deltaA = deltaA + TMP/D_A;
                            deltaB = deltaB + TMP/D_B;
                            deltaC = deltaC + TMP/D_C;
                            deltaD = deltaD + TMP/D_D;

                        end
                    end
                end

            end
        end
    end


    Ecorr_ccsd = cc_energy(t1,t2,sys);

    Ecorr_crcc23A = Ecorr_ccsd + deltaA;
    Ecorr_crcc23B = Ecorr_ccsd + deltaB;
    Ecorr_crcc23C = Ecorr_ccsd + deltaC;
    Ecorr_crcc23D = Ecorr_ccsd + deltaD;

    E_crcc23A = sys.Escf + Ecorr_ccsd + deltaA;
    E_crcc23B = sys.Escf + Ecorr_ccsd + deltaB;
    E_crcc23C = sys.Escf + Ecorr_ccsd + deltaC;
    E_crcc23D = sys.Escf + Ecorr_ccsd + deltaD;

    fprintf(' finished in %4.2f s\n',toc)

    fprintf('\n')
    fprintf('CR-CC(2,3)_A = %4.12f Ha     Ecorr = %4.12f Ha     Delta_A = %4.12f Ha\n',E_crcc23A,Ecorr_crcc23A,deltaA)
    fprintf('CR-CC(2,3)_B = %4.12f Ha     Ecorr = %4.12f Ha     Delta_B = %4.12f Ha\n',E_crcc23B,Ecorr_crcc23B,deltaB)
    fprintf('CR-CC(2,3)_C = %4.12f Ha     Ecorr = %4.12f Ha     Delta_C = %4.12f Ha\n',E_crcc23C,Ecorr_crcc23C,deltaC)
    fprintf('CR-CC(2,3)_D = %4.12f Ha     Ecorr = %4.12f Ha     Delta_D = %4.12f Ha\n',E_crcc23D,Ecorr_crcc23D,deltaD)
                   
end


