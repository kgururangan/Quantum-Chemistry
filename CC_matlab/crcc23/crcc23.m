function [Ecorr_crcc23A,Ecorr_crcc23B,Ecorr_crcc23C,Ecorr_crcc23D] = crcc23(t1,t2,L1,L2,HBar,sys)

    fprintf('\n==================================++Entering CR-CC(2,3) Routine++=============================\n')

    fprintf('\n')

    MM23 = build_MM23(t1,t2,HBar);
    L3 = build_L3_approx(L1,L2,HBar);

    deltaA = 0.0; % using MP denominator -(f_aa - f_ii + f_bb - f_jj + f_cc - f_kk)
    deltaB = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) | phi_{ijk}^{abc}>
    deltaC = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) | phi_{ijk}^{abc}>
    deltaD = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) + H_3(CCSD) | phi_{ijk}^{abc}>
    
    fprintf('\nCalculating MM correction...\n')
    tic

    for a = 1:sys.Nunocc
        for b = a+1:sys.Nunocc
            for c = b+1:sys.Nunocc
                for i = 1:sys.Nocc
                    for j = i+1:sys.Nocc
                        for k = j+1:sys.Nocc

                            temp = L3(i,j,k,a,b,c)*MM23(a,b,c,i,j,k);

                            DMP =    sys.fvv(a,a) + sys.fvv(b,b) + sys.fvv(c,c) ...
                                    -sys.foo(i,i) - sys.foo(j,j) - sys.foo(k,k);

                            [D1,D2,D3] = HBar_CCSD_T_diagonal(a,b,c,i,j,k,HBar);

                            D_A = -DMP;
                            D_B = -D1;
                            D_C = -(D1+D2);
                            D_D = -(D1+D2+D3);

                            deltaA = deltaA + temp/D_A;
                            deltaB = deltaB + temp/D_B;
                            deltaC = deltaC + temp/D_C;
                            deltaD = deltaD + temp/D_D;
                            
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

    fprintf('MM correction completed in %4.2f s\n',toc)

    fprintf('\n')
    fprintf('CR-CC(2,3)_A = %4.12f Ha     Delta_A = %4.12f Ha\n',Ecorr_crcc23A,deltaA)
    fprintf('CR-CC(2,3)_B = %4.12f Ha     Delta_B = %4.12f Ha\n',Ecorr_crcc23B,deltaB)
    fprintf('CR-CC(2,3)_C = %4.12f Ha     Delta_C = %4.12f Ha\n',Ecorr_crcc23C,deltaC)
    fprintf('CR-CC(2,3)_D = %4.12f Ha     Delta_D = %4.12f Ha \n',Ecorr_crcc23D,deltaD)
                



end

