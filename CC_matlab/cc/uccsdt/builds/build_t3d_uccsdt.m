function [X3D_abcijk] = build_t3d_uccsdt(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, sys)

    [HBar_t, VT3_t] = build_ucc_hbar_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys);

    H1B = HBar_t.H1B;
    H2B = HBar_t.H2B;
    H2C = HBar_t.H2C;
    VTC = VT3_t.C;

    I2C_vvov = H2C.vvov-einsum_kg(H1B.ov,t2c,'me,abim->abie');
    
    % MM23A
    M23_D1 =   -einsum_kg(H2C.vooo + VTC.vooo,t2c,'amij,bcmk->abcijk'); 
    M23_D1 = M23_D1 - permute(M23_D1,[1,2,3,6,5,4]) - permute(M23_D1,[1,2,3,4,6,5]) - permute(M23_D1,[3,2,1,4,5,6])...
                    - permute(M23_D1,[2,1,3,4,5,6]) + permute(M23_D1,[2,1,3,6,5,4]) + permute(M23_D1,[3,2,1,6,5,4]) ...
                    + permute(M23_D1,[2,1,3,4,6,5]) + permute(M23_D1,[3,2,1,4,6,5]);

         
    M23_D2 =   +einsum_kg(I2C_vvov + VTC.vvov,t2c,'abie,ecjk->abcijk');
    M23_D2 = M23_D2 - permute(M23_D2,[1,2,3,5,4,6]) - permute(M23_D2,[1,2,3,6,5,4]) ...
                    - permute(M23_D2,[3,2,1,4,5,6]) - permute(M23_D2,[1,3,2,4,5,6]) + permute(M23_D2,[3,2,1,5,4,6]) ...
                    + permute(M23_D2,[1,3,2,5,4,6]) + permute(M23_D2,[3,2,1,6,5,4]) + permute(M23_D2,[1,3,2,6,5,4]);
        
    MM23D = M23_D1 + M23_D2;

    % (HBar*T3)_C

    D1 = -einsum_kg(H1B.oo-diag(diag(sys.fb_oo)),t3d,'mk,abcijm->abcijk');

    D2 = einsum_kg(H1B.vv-diag(diag(sys.fb_vv)),t3d,'ce,abeijk->abcijk');

    D3 = 0.5*einsum_kg(H2C.oooo,t3d,'mnij,abcmnk->abcijk');

    D4 = 0.5*einsum_kg(H2C.vvvv,t3d,'abef,efcijk->abcijk');

    D5 = einsum_kg(H2B.ovvo,t3c,'maei,ebcmjk->abcijk');

    D6 = einsum_kg(H2C.voov,t3d,'amie,ebcmjk->abcijk');

    D13 = D1 + D3;
    D13 = D13 - permute(D13,[1,2,3,6,4,5]) - permute(D13,[1,2,3,4,6,5]);

    D24 = D2 + D4;
    D24 = D24 - permute(D24,[3,2,1,4,5,6]) - permute(D24,[1,3,2,4,5,6]);

    D56 = D5 + D6;
    D56 = D56 - permute(D56,[1,2,3,6,4,5]) - permute(D56,[1,2,3,5,4,6])...
              -permute(D56,[2,1,3,4,5,6]) - permute(D56,[3,2,1,4,5,6]) ...
              +permute(D56,[2,1,3,6,4,5]) + permute(D56,[2,1,3,5,4,6]) ...
              +permute(D56,[3,2,1,6,4,5]) + permute(D56,[3,2,1,5,4,6]); 
     
     X3D_abcijk = MM23D + D13 + D24 + D56;
     

end