function [X3B_abcijk] = build_t3b_uccsdt(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, sys)

    [HBar_t, VT3_t] = build_ucc_hbar_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys);

    H1A = HBar_t.H1A;
    H1B = HBar_t.H1B;
    H2A = HBar_t.H2A;
    H2B = HBar_t.H2B;
    H2C = HBar_t.H2C;
    VTA = VT3_t.A;
    VTB = VT3_t.B;
    VTC = VT3_t.C;    
    
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

    M23_34 = M23_D3 + M23_D4;
    M23_34 = M23_34 - permute(M23_34,[2,1,3,4,5,6]) - permute(M23_34,[1,2,3,5,4,6]) + permute(M23_34,[2,1,3,5,4,6]);
        
    M23_25 = M23_D2 + M23_D5;
    M23_25 = M23_25 - permute(M23_25,[1,2,3,5,4,6]);

    M23_16 = M23_D1 + M23_D6;
    M23_16 = M23_16 - permute(M23_16,[2,1,3,4,5,6]);
        
    MM23B = M23_16 + M23_25 + M23_34;

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
     
   
    X3B_abcijk = MM23B + D2 + D4 + D5  + D7 + D11 + D12 + D_3_8_13 + D_1_6_14 + D_9_10;

end
