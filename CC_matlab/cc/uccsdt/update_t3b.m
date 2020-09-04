function [t3b] = update_t3b(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, HBar_t, VT3_t, sys, shift)

    H1A = HBar_t.H1A;
    H1B = HBar_t.H1B;
    H2A = HBar_t.H2A;
    H2B = HBar_t.H2B;
    H2C = HBar_t.H2C;
    VTA = VT3_t.A;
    VTB = VT3_t.B;  
    
    % MM23B
    
    I2B_ovoo = H2B.ovoo + einsum_kg(H1A.ov,t2b,'me,ecjk->mcjk');
    I2B_vooo = H2B.vooo + einsum_kg(H1B.ov,t2b,'me,aeik->amik');
    I2A_vooo = H2A.vooo + einsum_kg(H1A.ov,t2a,'me,aeij->amij');    
   
    M23_D1 = einsum_kg(H2B.vvvo + VTB.vvvo,t2a,'bcek,aeij->abcijk');

    M23_D2 = -einsum_kg(I2B_ovoo + VTB.ovoo,t2a,'mcjk,abim->abcijk');
    
    M23_D3 = +einsum_kg(H2B.vvov + VTB.vvov,t2b,'acie,bejk->abcijk');

    M23_D4 = -einsum_kg(I2B_vooo + VTB.vooo,t2b,'amik,bcjm->abcijk');

    M23_D5 = +einsum_kg(H2A.vvov + VTA.vvov,t2b,'abie,ecjk->abcijk');

    M23_D6 = -einsum_kg(I2A_vooo + VTA.vooo,t2b,'amij,bcmk->abcijk');

%     M23_D1 = einsum_kg(H2B.vvvo,t2a,'bcek,aeij->abcijk');
% 
%     M23_D2 = -einsum_kg(I2B_ovoo,t2a,'mcjk,abim->abcijk');
%     
%     M23_D3 = +einsum_kg(H2B.vvov,t2b,'acie,bejk->abcijk');
% 
%     M23_D4 = -einsum_kg(I2B_vooo,t2b,'amik,bcjm->abcijk');
% 
%     M23_D5 = +einsum_kg(H2A.vvov,t2b,'abie,ecjk->abcijk');
% 
%     M23_D6 = -einsum_kg(I2A_vooo,t2b,'amij,bcmk->abcijk');
    
    M23_34 = M23_D3 + M23_D4;
    M23_34 = M23_34 - permute(M23_34,[2,1,3,4,5,6]) - permute(M23_34,[1,2,3,5,4,6]) + permute(M23_34,[2,1,3,5,4,6]);
        
    M23_25 = M23_D2 + M23_D5;
    M23_25 = M23_25 - permute(M23_25,[1,2,3,5,4,6]);

    M23_16 = M23_D1 + M23_D6;
    M23_16 = M23_16 - permute(M23_16,[2,1,3,4,5,6]);
        
    MM23B = M23_16 + M23_25 + M23_34;
%     
%     tmp = -0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,abfmjn->abej') - einsum_kg(sys.vB_oovv,t3b,'mnef,abfmjn->abej');
%     VT_1 = einsum_kg(tmp,t2b,'abej,ecik->abcijk');
%     VT_1 = VT_1 - permute(VT_1,[1,2,3,5,4,6]);
%     
%     tmp = -0.5*einsum_kg(sys.vA_oovv,t3b,'mnef,afcmnk->acek') - einsum_kg(sys.vB_oovv,t3c,'mnef,afcmnk->acek');
%     VT_2 = einsum_kg(tmp,t2a,'acek,ebij->abcijk');
%     VT_2 = VT_2 - permute(VT_2,[2,1,3,4,5,6]);
%     
%     tmp = -einsum_kg(sys.vB_oovv,t3b,'nmfe,bfcjnm->bcje') - 0.5*einsum_kg(sys.vC_oovv,t3c,'nmfe,bfcjnm->bcje');
%     VT_3 = einsum_kg(tmp,t2b,'bcje,aeik->abcijk');
%     VT_3 = VT_3 - permute(VT_3,[2,1,3,4,5,6]) - permute(VT_3,[1,2,3,5,4,6]) + permute(VT_3,[2,1,3,5,4,6]);
%         
%     tmp = 0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,ebfijn->mbij') + einsum_kg(sys.vB_oovv,t3b,'mnef,ebfijn->mbij');
%     VT_4 = -einsum_kg(tmp,t2b,'mbij,acmk->abcijk');
%     VT_4 = VT_4 - permute(VT_4,[2,1,3,4,5,6]);
%     
%     tmp = 0.5*einsum_kg(sys.vA_oovv,t3b,'mnef,efcink->mcik') + einsum_kg(sys.vB_oovv,t3c,'mnef,efcink->mcik');
%     VT_5 = -einsum_kg(tmp,t2a,'mcik,abmj->abcijk');
%     VT_5 = VT_5 - permute(VT_5,[1,2,3,5,4,6]);
%     
%     tmp = einsum_kg(sys.vB_oovv,t3b,'nmfe,afeink->amik') + 0.5*einsum_kg(sys.vC_oovv,t3c,'nmfe,afeink->amik');
%     VT_6 = -einsum_kg(tmp,t2b,'amik,bcjm->abcijk');
%     VT_6 = VT_6 - permute(VT_6,[2,1,3,4,5,6]) - permute(VT_6,[1,2,3,5,4,6]) + permute(VT_6,[2,1,3,5,4,6]);
%     
%     MM23B = MM23B + VT_1 + VT_2 + VT_3 + VT_4 + VT_5 + VT_6;

    % (HBar*T3)_C
    
    D1 = -einsum_kg(H1A.oo-diag(diag(sys.fa_oo)),t3b,'mi,abcmjk->abcijk');

    D2 = -einsum_kg(H1B.oo-diag(diag(sys.fb_oo)),t3b,'mk,abcijm->abcijk');

    D3 = einsum_kg(H1A.vv-diag(diag(sys.fa_vv)),t3b,'ae,ebcijk->abcijk');

    D4 = einsum_kg(H1B.vv-diag(diag(sys.fb_vv)),t3b,'ce,abeijk->abcijk');
    
    D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
    
    D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
    
    D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
    
    D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk');
    
    D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');    
    D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk');    
    D_9_10 = D9 + D10;
    D_9_10 = D_9_10 - permute(D_9_10,[2,1,3,4,5,6]) - permute(D_9_10,[1,2,3,5,4,6]) + permute(D_9_10,[2,1,3,5,4,6]);
    
    D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
    
    D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
    
    D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
    D_3_8_13 = D3 + D8 + D13;
    D_3_8_13 = D_3_8_13 - permute(D_3_8_13,[2,1,3,4,5,6]);
    
    D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
    D_1_6_14 = D1 + D6 + D14;
    D_1_6_14 = D_1_6_14 - permute(D_1_6_14,[1,2,3,5,4,6]);
     
   
    X3B_abcijk = MM23B + D2 + D4 + D5 + D7 + D11 + D12 + D_3_8_13 + D_1_6_14 + D_9_10;

%     D1 = -einsum_kg(H1A.oo-diag(diag(sys.fa_oo)),t3b,'mi,abcmjk->abcijk');
%     D1 = D1 - permute(D1,[1,2,3,5,4,6]);
%     
%     D2 = -einsum_kg(H1B.oo-diag(diag(sys.fb_oo)),t3b,'mk,abcijm->abcijk');
%     
%     D3 = einsum_kg(H1A.vv-diag(diag(sys.fa_vv)),t3b,'ae,ebcijk->abcijk');
%     D3 = D3 - permute(D3,[2,1,3,4,5,6]);
%     
%     D4 = einsum_kg(H1B.vv-diag(diag(sys.fb_vv)),t3b,'ce,abeijk->abcijk');
%     
%     D5 = 0.5*einsum_kg(H2A.oooo,t3b,'mnij,abcmnk->abcijk');
%     
%     D6 = einsum_kg(H2B.oooo,t3b,'mnjk,abcimn->abcijk');
%     D6 = D6 - permute(D6,[1,2,3,5,4,6]);
%     
%     D7 = 0.5*einsum_kg(H2A.vvvv,t3b,'abef,efcijk->abcijk');
%     
%     D8 = einsum_kg(H2B.vvvv,t3b,'bcef,aefijk->abcijk');
%     D8 = D8 - permute(D8,[2,1,3,4,5,6]);
%     
%     D9 = einsum_kg(H2A.voov,t3b,'amie,ebcmjk->abcijk');
%     D9 = D9 - permute(D9,[1,2,3,5,4,6]) - permute(D9,[2,1,3,4,5,6]) + permute(D9,[2,1,3,5,4,6]);
%     
%     D10 = einsum_kg(H2B.voov,t3c,'amie,becjmk->abcijk');
%     D10 = D10 - permute(D10,[1,2,3,5,4,6]) - permute(D10,[2,1,3,4,5,6]) + permute(D10,[2,1,3,5,4,6]);
%     
%     D11 = einsum_kg(H2B.ovvo,t3a,'mcek,abeijm->abcijk');
%     
%     D12 = einsum_kg(H2C.voov,t3b,'cmke,abeijm->abcijk');
%     
%     D13 = -einsum_kg(H2B.vovo,t3b,'amek,ebcijm->abcijk');
%     D13 = D13 - permute(D13,[2,1,3,4,5,6]);
%     
%     D14 = -einsum_kg(H2B.ovov,t3b,'mcie,abemjk->abcijk');
%     D14 = D14 - permute(D14,[1,2,3,5,4,6]);
%     
%     X3B_abcijk = MM23B + D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 ...
%                        + D11 + D12 + D13 + D14;

     for a = 1:sys.Nvir_alpha
        for b = a+1:sys.Nvir_alpha
            for c = 1:sys.Nvir_beta
                for i = 1:sys.Nocc_alpha
                    for j = i+1:sys.Nocc_alpha
                        for k = 1:sys.Nocc_beta

                            t3b(a,b,c,i,j,k) = X3B_abcijk(a,b,c,i,j,k)/...
                                (sys.fa_oo(i,i)+sys.fa_oo(j,j)+sys.fb_oo(k,k)-sys.fa_vv(a,a)-sys.fa_vv(b,b)-sys.fb_vv(c,c)-shift);
                            t3b(b,a,c,i,j,k) = -t3b(a,b,c,i,j,k);
                            t3b(a,b,c,j,i,k) = -t3b(a,b,c,i,j,k);
                            t3b(b,a,c,j,i,k) = t3b(a,b,c,i,j,k);

                        end
                    end
                end
            end
        end
    end


end

