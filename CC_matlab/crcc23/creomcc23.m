function [omega_crcc23A,omega_crcc23B,omega_crcc23C,omega_crcc23D,LM3] = creomcc23(t1,t2,R1,R2,L1,L2,r0,omega,HBar,sys)

    fprintf('\n==================================++Entering CR-EOMCC(2,3) Routine++=============================\n')

    fprintf('\n')

    %EOMMM23 = build_EOMMM23(t1,t2,R1,R2,HBar);
    MM23 = build_MM23(t1,t2,HBar);
    L3 = build_L3_approx(L1,L2,HBar);

    
    EOMMM23 = build_EOMMM23_v2(t1,t2,R1,R2,HBar,sys);
    

    deltaA = 0.0; % using MP denominator omega-(f_aa - f_ii + f_bb - f_jj + f_cc - f_kk)
    deltaB = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) | phi_{ijk}^{abc}>
    deltaC = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) | phi_{ijk}^{abc}>
    deltaD = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) + H_3(CCSD) | phi_{ijk}^{abc}>
    
    fprintf('MM Correction calculation...')
    tic
    
    LR = 0.0;
    
    LM3 = zeros(sys.Nocc,sys.Nocc,sys.Nocc,sys.Nunocc,sys.Nunocc,sys.Nunocc);

    for a = 1:sys.Nunocc
        for b = a+1:sys.Nunocc
            for c = b+1:sys.Nunocc
                for i = 1:sys.Nocc
                    for j = i+1:sys.Nocc
                        for k = j+1:sys.Nocc
                            
                            LR = LR + L3(i,j,k,a,b,c)*EOMMM23(a,b,c,i,j,k);

                            temp = L3(i,j,k,a,b,c)*(EOMMM23(a,b,c,i,j,k) + r0*MM23(a,b,c,i,j,k));
                            
                            LM3(i,j,k,a,b,c) = temp;

                            DMP =    sys.fvv(a,a) + sys.fvv(b,b) + sys.fvv(c,c) ...
                                    -sys.foo(i,i) - sys.foo(j,j) - sys.foo(k,k);

                            [D1,D2,D3] = HBar_CCSD_T_diagonal(a,b,c,i,j,k,HBar);

                            D_A = omega-DMP;
                            D_B = omega-D1;
                            D_C = omega-(D1+D2);
                            D_D = omega-(D1+D2+D3);

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

    omega_crcc23A = omega + deltaA;
    omega_crcc23B = omega + deltaB;
    omega_crcc23C = omega + deltaC;
    omega_crcc23D = omega + deltaD;

    Etot_ccsd = sys.Escf + Ecorr_ccsd + omega;
    Etot_crcc23A = sys.Escf + Ecorr_ccsd + omega_crcc23A;
    Etot_crcc23B = sys.Escf + Ecorr_ccsd + omega_crcc23B;
    Etot_crcc23C = sys.Escf + Ecorr_ccsd + omega_crcc23C;
    Etot_crcc23D = sys.Escf + Ecorr_ccsd + omega_crcc23D;

    fprintf(' finished in %4.2f s\n',toc)
    
    fprintf('<L3(2)*MM3(2)> = %4.12f\n',LR)
    
    HatoeV = 27.2113957;

    fprintf('\n')
    fprintf('EOMCCSD = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n\n',Etot_ccsd,omega,omega*HatoeV)
    fprintf('CR-EOMCC(2,3)_A = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n',Etot_crcc23A,omega_crcc23A,omega_crcc23A*HatoeV)
    fprintf('CR-EOMCC(2,3)_B = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n',Etot_crcc23B,omega_crcc23B,omega_crcc23B*HatoeV)
    fprintf('CR-EOMCC(2,3)_C = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n',Etot_crcc23C,omega_crcc23C,omega_crcc23C*HatoeV)
    fprintf('CR-EOMCC(2,3)_D = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n',Etot_crcc23D,omega_crcc23D,omega_crcc23D*HatoeV)
                



end
