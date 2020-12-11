function [deltaA,deltaB,deltaC,deltaD] = crucc23(cc_t,HBar_t,sys)

%     fprintf('\n==================================++Entering CR-UCC(2,3) Routine++=============================\n')
    
    tic_Start = tic;
    
    % get system dimensions
    Nocc_a = sys.Nocc_alpha; Nocc_b = sys.Nocc_beta;
    Nunocc_a = sys.Nvir_alpha; Nunocc_b = sys.Nvir_beta;
    
    % get T amplitudes
    t1a = cc_t.t1a; t1b = cc_t.t1b;
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
    
    % get 2-body HBar components
    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    
    % get the 3-body HBar triples diagonal
    [D3A_V,D3A_O,D3B_V,D3B_O,D3C_V,D3C_O,D3D_V,D3D_O] = get_triples_diagonal(cc_t,sys);
%     [d3a1,d3a2,d3b1,d3b2,d3c1,d3c2,d3d1,d3d2] = get_triples_diagonal_v2(cc_t,sys);

    % MM correction containers
    deltaA = 0.0; % using MP denominator -(f_aa - f_ii + f_bb - f_jj + f_cc - f_kk)
    deltaB = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) | phi_{ijk}^{abc}>
    deltaC = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) | phi_{ijk}^{abc}>
    deltaD = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) + H_3(CCSD) | phi_{ijk}^{abc}>
    
    LM3A = zeros(Nunocc_a,Nunocc_a,Nunocc_a,Nocc_a,Nocc_a,Nocc_a);
    LM3B = zeros(Nunocc_a,Nunocc_a,Nunocc_b,Nocc_a,Nocc_a,Nocc_b);
    LM3C = zeros(Nunocc_a,Nunocc_b,Nunocc_b,Nocc_a,Nocc_b,Nocc_b);
    LM3D = zeros(Nunocc_b,Nunocc_b,Nunocc_b,Nocc_b,Nocc_b,Nocc_b);
        

    % MM23A correction
    [MM23A] = build_MM23A(cc_t,HBar_t,sys);
    [L3A] = build_L3A_approx(cc_t,HBar_t,sys,0);
    ticA = tic;
    for a = 1:Nunocc_a
        for b = a+1:Nunocc_a
            for c = b+1:Nunocc_a
                for i = 1:Nocc_a
                    for j = i+1:Nocc_a
                        for k = j+1:Nocc_a

                            temp = L3A(i,j,k,a,b,c)*MM23A(a,b,c,i,j,k); 
                            
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

%                             D3 =  d3a1(a,j,i) + d3a1(a,k,i) + d3a1(a,k,j)...
%                                  +d3a1(b,j,i) + d3a1(b,k,i) + d3a1(b,k,j)...
%                                  +d3a1(c,j,i) + d3a1(c,k,i) + d3a1(c,k,j)...
%                                  -d3a2(b,a,i) - d3a2(b,a,j) - d3a2(b,a,k)...
%                                  -d3a2(c,a,i) - d3a2(c,a,j) - d3a2(c,a,k)...
%                                  -d3a2(c,b,i) - d3a2(c,b,j) - d3a2(c,b,k);
%                             D3 = -D3;
            

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
    %clear MM23A L3A
    fprintf('\nMM23A correction completed in %4.4fs\n\n',toc(ticA))
    
    % MM23B correction
    [MM23B] = build_MM23B(cc_t,HBar_t,sys);
    [L3B] = build_L3B_approx(cc_t,HBar_t,sys,0);
    ticB = tic;
    for a = 1:Nunocc_a
        for b = a+1:Nunocc_a
            for c = 1:Nunocc_b
                for i = 1:Nocc_a
                    for j = i+1:Nocc_a
                        for k = 1:Nocc_b

                            temp = L3B(i,j,k,a,b,c)*MM23B(a,b,c,i,j,k);
                            
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

%                             D3 =  d3a1(a,j,i) + d3b1(a,k,i) + d3b1(a,k,j)...
%                                  +d3a1(b,j,i) + d3b1(b,k,i) + d3b1(b,k,j)...
%                                  +d3c1(c,k,i) + d3c1(c,k,j)...
%                                  -d3a2(b,a,i) - d3a2(b,a,j)...
%                                  -d3b2(c,a,i) - d3b2(c,a,j) - d3c2(c,a,k)...
%                                  -d3b2(c,b,i) - d3b2(c,b,j) - d3c2(c,b,k);
%                             D3 = -D3;

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
    %clear MM23B L3B
    fprintf('MM23B correction completed in %4.4fs\n\n',toc(ticB))
    
    % closed shell
%     deltaA = 2*deltaA;
%     deltaB = 2*deltaB;
%     deltaC = 2*deltaC;
%     deltaD = 2*deltaD;
    
    % MM23C correction
    [MM23C] = build_MM23C(cc_t,HBar_t,sys);
    [L3C] = build_L3C_approx(cc_t,HBar_t,sys,0);
    ticC = tic;
    for a = 1:Nunocc_a
        for b = 1:Nunocc_b
            for c = b+1:Nunocc_b
                for i = 1:Nocc_a
                    for j = 1:Nocc_b
                        for k = j+1:Nocc_b

                            temp = L3C(i,j,k,a,b,c)*MM23C(a,b,c,i,j,k);
                            
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
    %clear MM23C L3C
    fprintf('MM23C correction completed in %4.4fs\n\n',toc(ticC))
    
    % MM23D correction
    [MM23D] = build_MM23D(cc_t,HBar_t,sys);
    [L3D] = build_L3D_approx(cc_t,HBar_t,sys,0);
    ticD = tic;
    for a = 1:Nunocc_b
        for b = a+1:Nunocc_b
            for c = b+1:Nunocc_b
                for i = 1:Nocc_b
                    for j = i+1:Nocc_b
                        for k = j+1:Nocc_b

                            temp = L3D(i,j,k,a,b,c)*MM23D(a,b,c,i,j,k);
                            
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
    %clear MM23D L3D
    fprintf('MM23D correction completed in %4.4fs\n\n',toc(ticD))
    
%     Ecorr_ccsd = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);

%     Ecorr_crcc23A = Ecorr_ccsd + deltaA;
%     Ecorr_crcc23B = Ecorr_ccsd + deltaB;
%     Ecorr_crcc23C = Ecorr_ccsd + deltaC;
%     Ecorr_crcc23D = Ecorr_ccsd + deltaD;
% 
%     E_crcc23A = sys.Escf + Ecorr_ccsd + deltaA;
%     E_crcc23B = sys.Escf + Ecorr_ccsd + deltaB;
%     E_crcc23C = sys.Escf + Ecorr_ccsd + deltaC;
%     E_crcc23D = sys.Escf + Ecorr_ccsd + deltaD;
% 
%     fprintf('\n')
%     fprintf('CR-UCC(2,3)_A = %4.12f Eh     Ecorr_A = %4.12f Eh     Delta_A = %4.12f Eh\n',E_crcc23A,Ecorr_crcc23A,deltaA)
%     fprintf('CR-UCC(2,3)_B = %4.12f Eh     Ecorr_B = %4.12f Eh     Delta_B = %4.12f Eh\n',E_crcc23B,Ecorr_crcc23B,deltaB)
%     fprintf('CR-UCC(2,3)_C = %4.12f Eh     Ecorr_C = %4.12f Eh     Delta_C = %4.12f Eh\n',E_crcc23C,Ecorr_crcc23C,deltaC)
%     fprintf('CR-UCC(2,3)_D = %4.12f Eh     Ecorr_D = %4.12f Eh     Delta_D = %4.12f Eh\n',E_crcc23D,Ecorr_crcc23D,deltaD)
    
    fprintf('CR-UCC(2,3) program completed in %4.4f s\n',toc(tic_Start));

end