function [deltaA,deltaB,deltaC,deltaD] = cct3_2(cc_t,HBar_t,Nact_h,Nact_p,sys)

    fprintf('\n==================================++Entering CR-UCC(t,3) Routine++=============================\n')
    
    tic_Start = tic;
    
    flag_debug_diagonal = false;
    
    % get size of unique spin-integrated triples space
    n1 = get_number_unique('3A',sys);
    n2 = get_number_unique('3B',sys);
    n3 = get_number_unique('3C',sys);
    n4 = get_number_unique('3D',sys);
    ntot = n1 + n2 + n3 + n4;
    
    % get system dimensions
    Nocc_a = sys.Nocc_alpha; Nocc_b = sys.Nocc_beta;
    Nunocc_a = sys.Nvir_alpha; Nunocc_b = sys.Nvir_beta;
    
    % active orbital ranges
    act_h_alpha_range = [sys.Nocc_alpha - Nact_h + 1 : sys.Nocc_alpha];
    act_h_beta_range = [sys.Nocc_beta - Nact_h + 1 : sys.Nocc_beta];
    act_p_alpha_range = [1 : Nact_p];
    act_p_beta_range = [1 : Nact_p];
    
    % get 2-body HBar components
    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    
    % get the 3-body HBar triples diagonal
    [D3A_V,D3A_O,D3B_V,D3B_O,D3C_V,D3C_O,D3D_V,D3D_O] = get_triples_diagonal(cc_t,sys);
    %[d3a1,d3a2,d3b1,d3b2,d3c1,d3c2,d3d1,d3d2] = get_triples_diagonal_v2(cc_t,sys);
    
    if flag_debug_diagonal
        [D3A_V,D3A_O,D3B_V,D3B_O,D3C_V,D3C_O,D3D_V,D3D_O] = get_triples_diagonal(cc_t,sys);

        fprintf('Error in D3A_O = %4.12f\n',get_error(d3a1,permute(D3A_O,[1,3,2])))
        fprintf('Error in D3A_V = %4.12f\n',get_error(d3a2,permute(D3A_V,[3,1,2])))
        fprintf('Error in D3B_O = %4.12f\n',get_error(d3b1,permute(D3B_O,[1,3,2])))
        fprintf('Error in D3B_V = %4.12f\n',get_error(d3b2,permute(D3B_V,[3,1,2])))
        fprintf('Error in D3C_O = %4.12f\n',get_error(d3c1,permute(D3C_O,[1,3,2])))
        fprintf('Error in D3C_V = %4.12f\n',get_error(d3c2,permute(D3C_V,[3,1,2])))
        fprintf('Error in D3D_O = %4.12f\n',get_error(d3d1,permute(D3D_O,[1,3,2])))
        fprintf('Error in D3D_V = %4.12f\n',get_error(d3d2,permute(D3D_V,[3,1,2])))
    end

    % MM correction containers
    deltaA = 0.0; % using MP denominator -(f_aa - f_ii + f_bb - f_jj + f_cc - f_kk)
    deltaB = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) | phi_{ijk}^{abc}>
    deltaC = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) | phi_{ijk}^{abc}>
    deltaD = 0.0; % using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) + H_3(CCSD) | phi_{ijk}^{abc}>
    
    corr_a = 0.0;
    corr_b = 0.0;
    corr_c = 0.0;
    corr_d = 0.0;
    
    LM3A = zeros(Nunocc_a,Nunocc_a,Nunocc_a,Nocc_a,Nocc_a,Nocc_a);
    LM3B = zeros(Nunocc_a,Nunocc_a,Nunocc_b,Nocc_a,Nocc_a,Nocc_b);
    LM3C = zeros(Nunocc_a,Nunocc_b,Nunocc_b,Nocc_a,Nocc_b,Nocc_b);
    LM3D = zeros(Nunocc_b,Nunocc_b,Nunocc_b,Nocc_b,Nocc_b,Nocc_b);
        
    ctP = 0;
    ctQ = 0;

    % MM23A correction
    [MM23A] = build_MM23A(cc_t,HBar_t,sys);
    [L3A] = build_L3A_approx(cc_t,HBar_t,sys,0);
    
%     MM23A = zero_MM_P_space(MM23A,p_space);
%     L3A = zero_L_P_space(L3A,p_space);
    
    ticA = tic;
    for a = 1:Nunocc_a
        for b = a+1:Nunocc_a
            for c = b+1:Nunocc_a
                for i = 1:Nocc_a
                    for j = i+1:Nocc_a
                        for k = j+1:Nocc_a
                                                                                  
                            num_act_h = num_active([i,j,k],act_h_alpha_range);
                            num_act_p = num_active([a,b,c],act_p_alpha_range);

                            if num_act_h >=2 && num_act_p >= 2

                                %fprintf('%d  %d  %d  %d  %d  %d\n',a+Nocc_a,b+Nocc_a,c+Nocc_a,i,j,k)
                                ctP = ctP + 1;
                                continue
                            else
                                
                                ctQ = ctQ + 1;

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
                                 
%                                 D3 =  d3a1(a,j,i) + d3a1(a,k,i) + d3a1(a,k,j)...
%                                      +d3a1(b,j,i) + d3a1(b,k,i) + d3a1(b,k,j)...
%                                      +d3a1(c,j,i) + d3a1(c,k,i) + d3a1(c,k,j)...
%                                      -d3a2(b,a,i) - d3a2(b,a,j) - d3a2(b,a,k)...
%                                      -d3a2(c,a,i) - d3a2(c,a,j) - d3a2(c,a,k)...
%                                      -d3a2(c,b,i) - d3a2(c,b,j) - d3a2(c,b,k);
%                                 D3 = -D3;

                                D_A = -DMP;
                                D_B = -D1;
                                D_C = -(D1+D2);
                                D_D = -(D1+D2+D3);

                                deltaA = deltaA + temp/D_A;
                                deltaB = deltaB + temp/D_B;
                                deltaC = deltaC + temp/D_C;
                                deltaD = deltaD + temp/D_D;
                                
                                corr_a = corr_a + temp/D_D;
                            end
                            
                        end
                    end
                end
            end
        end
    end
    %clear MM23A L3A
    fprintf('\nMM23A correction completed in %4.4fs\n\n',toc(ticA))
    
    %corr_a
    
    % MM23B correction
    [MM23B] = build_MM23B(cc_t,HBar_t,sys);
    [L3B] = build_L3B_approx(cc_t,HBar_t,sys,0);
    
%     MM23B = zero_MM_P_space(MM23B,p_space);
%     L3B = zero_L_P_space(L3B,p_space);
    
    ticB = tic;
    for a = 1:Nunocc_a
        for b = a+1:Nunocc_a
            for c = 1:Nunocc_b
                for i = 1:Nocc_a
                    for j = i+1:Nocc_a
                        for k = 1:Nocc_b
                            
                            num_act_h = num_active([i,j],act_h_alpha_range) + num_active([k],act_h_beta_range);
                            num_act_p = num_active([a,b],act_p_alpha_range) + num_active([c],act_p_beta_range);

                            if num_act_h >=2 && num_act_p >= 2
                                %fprintf('%d  %d  %d  %d  %d  %d\n',a+Nocc_a,b+Nocc_a,c+Nocc_b,i,j,k)
                                ctP = ctP + 1;
                                continue
                            else
                                
                                ctQ = ctQ + 1;

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

%                                 D3 =  d3a1(a,j,i) + d3b1(a,k,i) + d3b1(a,k,j)...
%                                      +d3a1(b,j,i) + d3b1(b,k,i) + d3b1(b,k,j)...
%                                      +d3c1(c,k,i) + d3c1(c,k,j)...
%                                      -d3a2(b,a,i) - d3a2(b,a,j)...
%                                      -d3b2(c,a,i) - d3b2(c,a,j) - d3c2(c,a,k)...
%                                      -d3b2(c,b,i) - d3b2(c,b,j) - d3c2(c,b,k);
%                                 D3 = -D3;

                                D_A = -DMP;
                                D_B = -D1;
                                D_C = -(D1+D2);
                                D_D = -(D1+D2+D3);

                                deltaA = deltaA + temp/D_A;
                                deltaB = deltaB + temp/D_B;
                                deltaC = deltaC + temp/D_C;
                                deltaD = deltaD + temp/D_D;
                                
                                corr_b = corr_b + temp/D_D;
                                
                            end
                            
                        end
                    end
                end
            end
        end
    end
    %clear MM23B L3B
    fprintf('MM23B correction completed in %4.4fs\n\n',toc(ticB))
    
    %corr_b
    
    % closed shell
%     deltaA = 2*deltaA;
%     deltaB = 2*deltaB;
%     deltaC = 2*deltaC;
%     deltaD = 2*deltaD;
    
    % MM23C correction
    [MM23C] = build_MM23C(cc_t,HBar_t,sys);
    [L3C] = build_L3C_approx(cc_t,HBar_t,sys,0);
    
    %MM23C = zero_MM_P_space(MM23C,p_space);
    %L3C = zero_L_P_space(L3C,p_space);
    
    ticC = tic;
    for a = 1:Nunocc_a
        for b = 1:Nunocc_b
            for c = b+1:Nunocc_b
                for i = 1:Nocc_a
                    for j = 1:Nocc_b
                        for k = j+1:Nocc_b

                            num_act_h = num_active([i],act_h_alpha_range) + num_active([j,k],act_h_beta_range);
                            num_act_p = num_active([a],act_p_alpha_range) + num_active([b,c],act_p_beta_range);

                            if num_act_h >=2 && num_act_p >= 2
                                %fprintf('%d  %d  %d  %d  %d  %d\n',a+Nocc_a,b+Nocc_b,c+Nocc_b,i,j,k)
                                ctP = ctP + 1;
                                continue
                            else
                                
                                ctQ = ctQ + 1;
                                
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
                
%                                 D3 =  d3b1(a,j,i) + d3b1(a,k,i) ...
%                                      +d3c1(b,j,i) + d3c1(b,k,i) + d3d1(b,k,j)...
%                                      +d3c1(c,j,i) + d3c1(c,k,i) + d3d1(c,k,j)...
%                                      -d3b2(b,a,i) - d3c2(b,a,j) - d3c2(b,a,k)...
%                                      -d3b2(c,a,i) - d3c2(c,a,j) - d3c2(c,a,k)...
%                                      -d3d2(c,b,j) - d3d2(c,b,k);
%                                 D3 = -D3;

                                D_A = -DMP;
                                D_B = -D1;
                                D_C = -(D1+D2);
                                D_D = -(D1+D2+D3);

                                deltaA = deltaA + temp/D_A;
                                deltaB = deltaB + temp/D_B;
                                deltaC = deltaC + temp/D_C;
                                deltaD = deltaD + temp/D_D;
                                
                                corr_c = corr_c + temp/D_D;
                                
                            end
                            
                        end
                    end
                end
            end
        end
    end
    %clear MM23C L3C
    fprintf('MM23C correction completed in %4.4fs\n\n',toc(ticC))
    
    %corr_c
    
    % MM23D correction
    [MM23D] = build_MM23D(cc_t,HBar_t,sys);
    [L3D] = build_L3D_approx(cc_t,HBar_t,sys,0);
    
    %MM23D = zero_MM_P_space(MM23D,p_space);
    %L3D = zero_L_P_space(L3D,p_space);
    
    
    ticD = tic;
    for a = 1:Nunocc_b
        for b = a+1:Nunocc_b
            for c = b+1:Nunocc_b
                for i = 1:Nocc_b
                    for j = i+1:Nocc_b
                        for k = j+1:Nocc_b

                            num_act_h = num_active([i,j,k],act_h_beta_range);
                            num_act_p = num_active([a,b,c],act_p_beta_range);

                            if num_act_h >=2 && num_act_p >= 2
                                %fprintf('%d  %d  %d  %d  %d  %d\n',a+Nocc_b,b+Nocc_b,c+Nocc_b,i,j,k)
                                ctP = ctP + 1;
                                continue
                            else
                                
                                ctQ = ctQ + 1;
                                
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

%                                 D3 =  d3d1(a,j,i) + d3d1(a,k,i) + d3d1(a,k,j)...
%                                      +d3d1(b,j,i) + d3d1(b,k,i) + d3d1(b,k,j)...
%                                      +d3d1(c,j,i) + d3d1(c,k,i) + d3d1(c,k,j)...
%                                      -d3d2(b,a,i) - d3d2(b,a,j) - d3d2(b,a,k)...
%                                      -d3d2(c,a,i) - d3d2(c,a,j) - d3d2(c,a,k)...
%                                      -d3d2(c,b,i) - d3d2(c,b,j) - d3d2(c,b,k);
%                                 D3 = -D3;

                                D_A = -DMP;
                                D_B = -D1;
                                D_C = -(D1+D2);
                                D_D = -(D1+D2+D3);

                                deltaA = deltaA + temp/D_A;
                                deltaB = deltaB + temp/D_B;
                                deltaC = deltaC + temp/D_C;
                                deltaD = deltaD + temp/D_D;
                                
                                corr_d = corr_d + temp/D_D;
                                
                            end
                            
                        end
                    end
                end
            end
        end
    end
    %clear MM23D L3D
    fprintf('MM23D correction completed in %4.4fs\n\n',toc(ticD))
    
    %corr_d
    
    fprintf('CC(t;3) correction completed in %4.4f s\n',toc(tic_Start));
    fprintf('\n%d P space triples\n%d Q space triples\nSum = %d\n# Triples w/o symmetry = %d\n',ctP,ctQ,ctP+ctQ,ntot)

end

function [flag] = is_P(x,act_range)

    nact = length(find(x >= act_range(1) & x <= act_range(end)));
    if nact >= 1
        flag = true;
    else
        flag = false;
    end

end


function [MM3] = zero_MM_P_space(MM3,p_space)

    for a = 1:size(MM3,1)
        for b = 1:size(MM3,2)
            for c = 1:size(MM3,3)
                for i = 1:size(MM3,4)
                    for j = size(MM3,5)
                        for k = size(MM3,6)
                            if p_space(a,b,c,i,j,k) == 1
                                MM3(a,b,c,i,j,k) = 0;
                            end
                        end
                    end
                end
            end
        end
    end
    
end

function [L] = zero_L_P_space(L,p_space)

    for i = 1:size(L,1)
        for j = 1:size(L,2)
            for k = 1:size(L,3)
                for a = 1:size(L,4)
                    for b = size(L,5)
                        for c = size(L,6)
                            if p_space(a,b,c,i,j,k) == 1
                                L(i,j,k,a,b,c) = 0;
                            end
                        end
                    end
                end
            end
        end
    end  
end

function [val] = num_active(x,act_range)
    val = length(find(x >= act_range(1) & x <= act_range(end)));
end