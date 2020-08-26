function [X3A_abcijk] = build_t3a_uccsdt(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, sys)

    [HBar_t, VT3_t] = build_ucc_hbar_intermediates(t1a,t1b,t2a,t2b,t2c,t3a,t3b,t3c,t3d,sys);

    H1A = HBar_t.H1A;
    H1B = HBar_t.H1B;
    H2A = HBar_t.H2A;
    H2B = HBar_t.H2B;
    H2C = HBar_t.H2C;
    VTA = VT3_t.A;
    VTB = VT3_t.B;
    VTC = VT3_t.C;    

    H1A = HBar_t.H1A;
    H2A = HBar_t.H2A;
    H2B = HBar_t.H2B;
    VTA = VT3_t.A;

    I12 = H2A.vvov-einsum_kg(H1A.ov,t2a,'me,abim->abie');
    
    % MM23A
    M23_D1 =   -einsum_kg(H2A.vooo + VTA.vooo,t2a,'amij,bcmk->abcijk'); 
    M23_D1 = M23_D1 - permute(M23_D1,[1,2,3,6,5,4]) - permute(M23_D1,[1,2,3,4,6,5]) - permute(M23_D1,[3,2,1,4,5,6])...
                    - permute(M23_D1,[2,1,3,4,5,6]) + permute(M23_D1,[2,1,3,6,5,4]) + permute(M23_D1,[3,2,1,6,5,4]) ...
                    + permute(M23_D1,[2,1,3,4,6,5]) + permute(M23_D1,[3,2,1,4,6,5]);

         
    M23_D2 =   +einsum_kg(I12 + VTA.vvov,t2a,'abie,ecjk->abcijk');
    M23_D2 = M23_D2 - permute(M23_D2,[1,2,3,5,4,6]) - permute(M23_D2,[1,2,3,6,5,4]) ...
                    - permute(M23_D2,[3,2,1,4,5,6]) - permute(M23_D2,[1,3,2,4,5,6]) + permute(M23_D2,[3,2,1,5,4,6]) ...
                    + permute(M23_D2,[1,3,2,5,4,6]) + permute(M23_D2,[3,2,1,6,5,4]) + permute(M23_D2,[1,3,2,6,5,4]);
        
    MM23A = M23_D1 + M23_D2;

    % (HBar*T3)_C
    
    D1 = -einsum_kg(H1A.oo-diag(diag(sys.fa_oo)),t3a,'mk,abcijm->abcijk');
    
    D2 = einsum_kg(H1A.vv-diag(diag(sys.fa_vv)),t3a,'ce,abeijk->abcijk');
    
    D3 = 0.5*einsum_kg(H2A.oooo,t3a,'mnij,abcmnk->abcijk');
    D13 = D1 + D3;
    D13 = D13 - permute(D13,[1,2,3,6,5,4]) - permute(D13,[1,2,3,4,6,5]);
    
    D4 = 0.5*einsum_kg(H2A.vvvv,t3a,'abef,efcijk->abcijk');
    D24 = D2 + D4;
    D24 = D24 - permute(D24,[3,2,1,4,5,6]) - permute(D24,[1,3,2,4,5,6]);
    
    D5 = einsum_kg(H2A.voov,t3a,'cmke,abeijm->abcijk');
    D6 = einsum_kg(H2B.voov,t3b,'cmke,abeijm->abcijk');
    D56 = D5 + D6;

    D56 = D56 - permute(D56,[1,2,3,6,5,4]) - permute(D56,[1,2,3,4,6,5])...
              - permute(D56,[3,2,1,4,5,6]) - permute(D56,[1,3,2,4,5,6]) ...
              + permute(D56,[3,2,1,6,5,4]) + permute(D56,[3,2,1,4,6,5]) ...
              + permute(D56,[1,3,2,6,5,4]) + permute(D56,[1,3,2,4,6,5]);    
     
     X3A_abcijk = MM23A + D13 + D24 + D56;
end