function [deltaA,deltaB,deltaC,deltaD] = creomucc23(cc_t,omega,HBar_t,sys,iroot)

    tic_Start = tic;
    
    % get system dimensions
    Nocc_a = sys.Nocc_alpha; Nocc_b = sys.Nocc_beta;
    Nunocc_a = sys.Nvir_alpha; Nunocc_b = sys.Nvir_beta;
    
    % get 2-body HBar components
    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    
    % CCSD correlation energy
%     [Ecorr_ccsd] = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);
    
    % get the 3-body HBar triples diagonal
    [D3A_V,D3A_O,D3B_V,D3B_O,D3C_V,D3C_O,D3D_V,D3D_O] = get_triples_diagonal(cc_t,sys);

    % MM correction containers
    deltaA = 0.0; % using MP denominator -(f_aa - f_ii + f_bb - f_jj + f_cc - f_kk)
    deltaB = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) | phi_{ijk}^{abc}>
    deltaC = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) | phi_{ijk}^{abc}>
    deltaD = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) + H_3(CCSD) | phi_{ijk}^{abc}>
        
    r0 = cc_t.r0(iroot);

    % MM23A correction
    
    [EOMMM23A] = build_EOMMM23A(cc_t,HBar_t,iroot);
    [MM23A] = build_MM23A(cc_t,HBar_t,sys);
    [L3A] = build_L3A_approx(cc_t,HBar_t,sys,iroot);
    
    LM3A = zeros(Nunocc_a,Nunocc_a,Nunocc_a,Nocc_a,Nocc_a,Nocc_a);
    
    ticA = tic;
    for a = 1:Nunocc_a
        for b = a+1:Nunocc_a
            for c = b+1:Nunocc_a
                for i = 1:Nocc_a
                    for j = i+1:Nocc_a
                        for k = j+1:Nocc_a

                            temp = L3A(i,j,k,a,b,c)*(EOMMM23A(a,b,c,i,j,k) + r0*MM23A(a,b,c,i,j,k));
                            
                            LM3A(a,b,c,i,j,k) = temp;                        

                            DMP = sys.fa_vv(a,a) + sys.fa_vv(b,b) + sys.fa_vv(c,c) ...
                                  -sys.fa_oo(i,i) - sys.fa_oo(j,j) - sys.fa_oo(k,k); 
                                
                            D1 = H1A.vv(a,a) + H1A.vv(b,b) + H1A.vv(c,c) ...
                                -H1A.oo(i,i) - H1A.oo(j,j) - H1A.oo(k,k);

                            D2 = -H2A.voov(a,i,i,a)-H2A.voov(b,i,i,b)-H2A.voov(c,i,i,c)...
                                 -H2A.voov(a,j,j,a)-H2A.voov(b,j,j,b)-H2A.voov(c,j,j,c)...
                                 -H2A.voov(a,k,k,a)-H2A.voov(b,k,k,b)-H2A.voov(c,k,k,c)...
                                 -H2A.oooo(j,i,j,i)-H2A.oooo(k,i,k,i)-H2A.oooo(k,j,k,j)...
                                 -H2A.vvvv(b,a,b,a)-H2A.vvvv(c,a,c,a)-H2A.vvvv(c,b,c,b);
                            D2 = -D2;
                            
                            D3 = -D3A_O(a,i,j)-D3A_O(a,i,k)-D3A_O(a,j,k)...
                                 -D3A_O(b,i,j)-D3A_O(b,i,k)-D3A_O(b,j,k)...
                                 -D3A_O(c,i,j)-D3A_O(c,i,k)-D3A_O(c,j,k)...
                                 +D3A_V(a,i,b)+D3A_V(a,i,c)+D3A_V(b,i,c)...
                                 +D3A_V(a,j,b)+D3A_V(a,j,c)+D3A_V(b,j,c)...
                                 +D3A_V(a,k,b)+D3A_V(a,k,c)+D3A_V(b,k,c);

                            D_A = omega(iroot)-DMP;
                            D_B = omega(iroot)-D1;
                            D_C = omega(iroot)-(D1+D2);
                            D_D = omega(iroot)-(D1+D2+D3);

%                             if abs(temp) > 1e-6
%                                 fprintf('\n%d  %d  %d->  %d %d %d\n',i,j,k,a+Nocc_a,b+Nocc_a,c+Nocc_a)
%                                 fprintf('M3A*L3A = %4.12f\n',temp)
%                                 fprintf('M3A = %4.12f\n',EOMMM23A(a,b,c,i,j,k) + r0*MM23A(a,b,c,i,j,k))
%                                 fprintf('L3A = %4.12f\n',L3A(i,j,k,a,b,c))
%                                 fprintf('D_A = %4.12f\n',D_A)
%                                 fprintf('D_B = %4.12f\n',D_B)
%                                 fprintf('D_C = %4.12f\n',D_C)
%                                 fprintf('D_D = %4.12f\n',D_D)
%                             end
        

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
    %clear EOMMM23A MM23A L3A
    fprintf('EOM-MM23A correction completed in %4.4fs\n\n',toc(ticA))

    % MM23B correction

    [EOMMM23B] = build_EOMMM23B(cc_t,HBar_t,iroot);
    [MM23B] = build_MM23B(cc_t,HBar_t,sys);
    [L3B] = build_L3B_approx(cc_t,HBar_t,sys,iroot);
    
    LM3B = zeros(Nunocc_a,Nunocc_a,Nunocc_b,Nocc_a,Nocc_a,Nocc_b);

    ticB = tic;
    for a = 1:Nunocc_a
        for b = a+1:Nunocc_a
            for c = 1:Nunocc_b
                for i = 1:Nocc_a
                    for j = i+1:Nocc_a
                        for k = 1:Nocc_b

                            temp = L3B(i,j,k,a,b,c)*(EOMMM23B(a,b,c,i,j,k) + r0*MM23B(a,b,c,i,j,k));
                            
                            LM3B(a,b,c,i,j,k) = temp;                           

                            DMP = sys.fa_vv(a,a) + sys.fa_vv(b,b) + sys.fb_vv(c,c) ...
                                  -sys.fa_oo(i,i) - sys.fa_oo(j,j) - sys.fb_oo(k,k);
                                
                            D1 = H1A.vv(a,a) + H1A.vv(b,b) + H1B.vv(c,c) ...
                                -H1A.oo(i,i) - H1A.oo(j,j) - H1B.oo(k,k);
                                
                            D2 = -H2A.voov(a,i,i,a)-H2A.voov(b,i,i,b)+H2B.vovo(c,i,c,i)...
                                 -H2A.voov(a,j,j,a)-H2A.voov(b,j,j,b)+H2B.vovo(c,j,c,j)...
                                 +H2B.ovov(k,a,k,a)+H2B.ovov(k,b,k,b)-H2C.voov(c,k,k,c)...
                                 -H2A.oooo(j,i,j,i)-H2B.oooo(k,i,k,i)-H2B.oooo(k,j,k,j)...
                                 -H2A.vvvv(b,a,b,a)-H2B.vvvv(c,a,c,a)-H2B.vvvv(c,b,c,b);
                            D2 = -D2;
                            
                            D3 = -D3A_O(a,i,j)-D3B_O(a,i,k)-D3B_O(a,j,k)...
                                 -D3A_O(b,i,j)-D3B_O(b,i,k)-D3B_O(b,j,k)...
                                 -D3C_O(c,i,k)-D3C_O(c,j,k)...
                                 +D3A_V(a,i,b)+D3B_V(a,i,c)+D3B_V(b,i,c)...
                                 +D3A_V(a,j,b)+D3B_V(a,j,c)+D3B_V(b,j,c)...
                                 +D3C_V(a,k,c)+D3C_V(b,k,c);

                            D_A = omega(iroot)-DMP;
                            D_B = omega(iroot)-D1;
                            D_C = omega(iroot)-(D1+D2);
                            D_D = omega(iroot)-(D1+D2+D3);

%                             if abs(temp) > 1e-6
%                                 fprintf('\n%d  %d  %d->  %d %d %d\n',i,j,k,a+Nocc_a,b+Nocc_a,c+Nocc_b)
%                                 fprintf('M3B*L3B = %4.12f\n',temp)
%                                 fprintf('M3B = %4.12f\n',EOMMM23B(a,b,c,i,j,k) + r0*MM23B(a,b,c,i,j,k))
%                                 fprintf('L3B = %4.12f\n',L3B(i,j,k,a,b,c))
%                                 fprintf('D_A = %4.12f\n',D_A)
%                                 fprintf('D_B = %4.12f\n',D_B)
%                                 fprintf('D_C = %4.12f\n',D_C)
%                                 fprintf('D_D = %4.12f\n',D_D)
%                             end

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
    %clear MM23B L3B
    fprintf('EOM-MM23B correction completed in %4.4fs\n\n',toc(ticB))
    
    % MM23C correction

    [EOMMM23C] = build_EOMMM23C(cc_t,HBar_t,iroot);
    [MM23C] = build_MM23C(cc_t,HBar_t,sys);
    [L3C] = build_L3C_approx(cc_t,HBar_t,sys,iroot);
    
    LM3C = zeros(Nunocc_a,Nunocc_b,Nunocc_b,Nocc_a,Nocc_b,Nocc_b);

    ticC = tic;
    for a = 1:Nunocc_a
        for b = 1:Nunocc_b
            for c = b+1:Nunocc_b
                for i = 1:Nocc_a
                    for j = 1:Nocc_b
                        for k = j+1:Nocc_b

                            temp = L3C(i,j,k,a,b,c)*(EOMMM23C(a,b,c,i,j,k) + r0*MM23C(a,b,c,i,j,k));
                            
                            LM3C(a,b,c,i,j,k) = temp;

                            DMP = sys.fa_vv(a,a) + sys.fb_vv(b,b) + sys.fb_vv(c,c) ...
                                  -sys.fa_oo(i,i) - sys.fb_oo(j,j) - sys.fb_oo(k,k);
                                
                            D1 = H1A.vv(a,a) + H1B.vv(b,b) + H1B.vv(c,c) ...
                                -H1A.oo(i,i) - H1B.oo(j,j) - H1B.oo(k,k);

                            D2 = -H2A.voov(a,i,i,a)+H2B.vovo(b,i,b,i)+H2B.vovo(c,i,c,i)...
                                 +H2B.ovov(j,a,j,a)-H2C.voov(b,j,j,b)-H2C.voov(c,j,j,c)...
                                 +H2B.ovov(k,a,k,a)-H2C.voov(b,k,k,b)-H2C.voov(c,k,k,c)...
                                 -H2B.oooo(j,i,j,i)-H2B.oooo(k,i,k,i)-H2C.oooo(k,j,k,j)...
                                 -H2B.vvvv(b,a,b,a)-H2B.vvvv(c,a,c,a)-H2C.vvvv(c,b,c,b);
                            D2 = -D2;
                            
                            D3 = -D3B_O(a,i,j)-D3B_O(a,i,k)...
                                 -D3C_O(b,i,j)-D3C_O(b,i,k)-D3D_O(b,j,k)...
                                 -D3C_O(c,i,j)-D3C_O(c,i,k)-D3D_O(c,j,k)...
                                 +D3B_V(a,i,b)+D3B_V(a,i,c)...
                                 +D3C_V(a,j,b)+D3C_V(a,j,c)+D3D_V(b,j,c)...
                                 +D3C_V(a,k,b)+D3C_V(a,k,c)+D3D_V(b,k,c);

                            D_A = omega(iroot)-DMP;
                            D_B = omega(iroot)-D1;
                            D_C = omega(iroot)-(D1+D2);
                            D_D = omega(iroot)-(D1+D2+D3);

%                             if abs(temp) > 1e-6
%                                 fprintf('\n%d  %d  %d->  %d %d %d\n',i,j,k,a+Nocc_a,b+Nocc_b,c+Nocc_b)
%                                 fprintf('M3C*L3C = %4.12f\n',temp)
%                                 fprintf('M3C = %4.12f\n',EOMMM23C(a,b,c,i,j,k) + r0*MM23C(a,b,c,i,j,k))
%                                 fprintf('L3C = %4.12f\n',L3C(i,j,k,a,b,c))
%                                 fprintf('D_A = %4.12f\n',D_A)
%                                 fprintf('D_B = %4.12f\n',D_B)
%                                 fprintf('D_C = %4.12f\n',D_C)
%                                 fprintf('D_D = %4.12f\n',D_D)
%                             end

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
    %clear MM23C L3C
    fprintf('MM23C correction completed in %4.4fs\n\n',toc(ticC))
    
    % MM23D correction

    [EOMMM23D] = build_EOMMM23D(cc_t,HBar_t,iroot);
    [MM23D] = build_MM23D(cc_t,HBar_t,sys);
    [L3D] = build_L3D_approx(cc_t,HBar_t,sys,iroot);
    
    LM3D = zeros(Nunocc_b,Nunocc_b,Nunocc_b,Nocc_b,Nocc_b,Nocc_b);

    ticD = tic;
    for a = 1:Nunocc_b
        for b = a+1:Nunocc_b
            for c = b+1:Nunocc_b
                for i = 1:Nocc_b
                    for j = i+1:Nocc_b
                        for k = j+1:Nocc_b

                            temp = L3D(i,j,k,a,b,c)*(EOMMM23D(a,b,c,i,j,k) + r0*MM23D(a,b,c,i,j,k));
                            
                            LM3D(a,b,c,i,j,k) = temp;

                            DMP = sys.fb_vv(a,a) + sys.fb_vv(b,b) + sys.fb_vv(c,c) ...
                                  -sys.fb_oo(i,i) - sys.fb_oo(j,j) - sys.fb_oo(k,k);
                                
                            D1 = H1B.vv(a,a) + H1B.vv(b,b) + H1B.vv(c,c) ...
                                -H1B.oo(i,i) - H1B.oo(j,j) - H1B.oo(k,k);

                            D2 = -H2C.voov(a,i,i,a)-H2C.voov(b,i,i,b)-H2C.voov(c,i,i,c)...
                                 -H2C.voov(a,j,j,a)-H2C.voov(b,j,j,b)-H2C.voov(c,j,j,c)...
                                 -H2C.voov(a,k,k,a)-H2C.voov(b,k,k,b)-H2C.voov(c,k,k,c)...
                                 -H2C.oooo(j,i,j,i)-H2C.oooo(k,i,k,i)-H2C.oooo(k,j,k,j)...
                                 -H2C.vvvv(b,a,b,a)-H2C.vvvv(c,a,c,a)-H2C.vvvv(c,b,c,b);
                            D2 = -D2;
                            
                            D3 = -D3D_O(a,i,j)-D3D_O(a,i,k)-D3D_O(a,j,k)...
                                 -D3D_O(b,i,j)-D3D_O(b,i,k)-D3D_O(b,j,k)...
                                 -D3D_O(c,i,j)-D3D_O(c,i,k)-D3D_O(c,j,k)...
                                 +D3D_V(a,i,b)+D3D_V(a,i,c)+D3D_V(b,i,c)...
                                 +D3D_V(a,j,b)+D3D_V(a,j,c)+D3D_V(b,j,c)...
                                 +D3D_V(a,k,b)+D3D_V(a,k,c)+D3D_V(b,k,c);

                            D_A = omega(iroot)-DMP;
                            D_B = omega(iroot)-D1;
                            D_C = omega(iroot)-(D1+D2);
      % 
%                             if abs(temp) > 1e-6
%                                 fprintf('\n%d  %d  %d->  %d %d %d\n',i,j,k,a+Nocc_b,b+Nocc_b,c+Nocc_b)
%                                 fprintf('M3D*L3D = %4.12f\n',temp)
%                                 fprintf('M3D = %4.12f\n',EOMMM23D(a,b,c,i,j,k) + r0*MM23D(a,b,c,i,j,k))
%                                 fprintf('L3D = %4.12f\n',L3D(i,j,k,a,b,c))
%                                 fprintf('D_A = %4.12f\n',D_A)
%                                 fprintf('D_B = %4.12f\n',D_B)
%                                 fprintf('D_C = %4.12f\n',D_C)
%                                 fprintf('D_D = %4.12f\n',D_D)
%                             end
                      D_D = omega(iroot)-(D1+D2+D3);

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
    %clear MM23D L3D
    fprintf('MM23D correction completed in %4.4fs\n\n',toc(ticD))

%     omega_crcc23A = omega(iroot) + deltaA;
%     omega_crcc23B = omega(iroot) + deltaB;
%     omega_crcc23C = omega(iroot) + deltaC;
%     omega_crcc23D = omega(iroot) + deltaD;

%     Etot_ccsd = sys.Escf + Ecorr_ccsd + omega(iroot);
%     Etot_crcc23A = sys.Escf + Ecorr_ccsd + omega_crcc23A;
%     Etot_crcc23B = sys.Escf + Ecorr_ccsd + omega_crcc23B;
%     Etot_crcc23C = sys.Escf + Ecorr_ccsd + omega_crcc23C;
%     Etot_crcc23D = sys.Escf + Ecorr_ccsd + omega_crcc23D;
%     
%     HatoeV = 27.2113957;
%
%     fprintf('\n')
%     fprintf('EOMCCSD = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n\n',Etot_ccsd,omega(iroot),omega(iroot)*HatoeV)
%     fprintf('CR-EOMCC(2,3)_A = %4.12f Eh     VEE (rel. CCSD) = %4.12f Eh   (%4.8f eV)\n',Etot_crcc23A,omega_crcc23A,omega_crcc23A*HatoeV)
%     fprintf('CR-EOMCC(2,3)_B = %4.12f Eh     VEE (rel. CCSD) = %4.12f Eh   (%4.8f eV)\n',Etot_crcc23B,omega_crcc23B,omega_crcc23B*HatoeV)
%     fprintf('CR-EOMCC(2,3)_C = %4.12f Eh     VEE (rel. CCSD) = %4.12f Eh   (%4.8f eV)\n',Etot_crcc23C,omega_crcc23C,omega_crcc23C*HatoeV)
%     fprintf('CR-EOMCC(2,3)_D = %4.12f Eh     VEE (rel. CCSD) = %4.12f Eh   (%4.8f eV)\n',Etot_crcc23D,omega_crcc23D,omega_crcc23D*HatoeV)
                

    fprintf('CR-EOMCC(2,3) calculation completed in %4.4f s\n',toc(tic_Start));

end
