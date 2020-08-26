function [t3c] = update_t3c(t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d, HBar_t, VT3_t, sys, shift)

    H1A = HBar_t.H1A;
    H1B = HBar_t.H1B;
    H2A = HBar_t.H2A;
    H2B = HBar_t.H2B;
    H2C = HBar_t.H2C;
    VTB = VT3_t.B;  
    VTC = VT3_t.C;
    
    % MM23C
    I2B_vooo = H2B.vooo + einsum_kg(H1B.ov,t2b,'me,aeij->amij');
    I2C_ovoo = H2C.ovoo + einsum_kg(H1B.ov,t2c,'me,ecjk->mcjk');
    I2B_ovoo = H2B.ovoo + einsum_kg(H1A.ov,t2b,'me,ebij->mbij');
    
    M23_D1 = +einsum_kg(H2B.vvov + VTB.vvov,t2c,'abie,ecjk->abcijk');

    M23_D2 = -einsum_kg(I2B_vooo + VTB.vooo,t2c,'amij,bcmk->abcijk');
        
    M23_D3 = +einsum_kg(H2C.vvvo + VTC.vvvo,t2b,'bcek,aeij->abcijk');

    M23_D4 = -einsum_kg(I2C_ovoo + VTC.ovoo,t2b,'mcjk,abim->abcijk');

    M23_D5 = +einsum_kg(H2B.vvvo + VTB.vvvo,t2b,'abej,ecik->abcijk');

    M23_D6 = -einsum_kg(I2B_ovoo + VTB.ovoo,t2b,'mbij,acmk->abcijk');

    M23_14 = M23_D1 + M23_D4;
    M23_14 = M23_14 - permute(M23_14,[1,3,2,4,5,6]);

    M23_23 = M23_D2 + M23_D3;
    M23_23 = M23_23 - permute(M23_23,[1,2,3,4,6,5]);

    M23_56 = M23_D5 + M23_D6;
    M23_56 = M23_56 - permute(M23_56,[1,3,2,4,5,6]) - permute(M23_56,[1,2,3,4,6,5]) + permute(M23_56,[1,3,2,4,6,5]);
    
    MM23C = M23_14 + M23_23 + M23_56;

    % (HBar*T3)_C
    
    D1 = -einsum_kg(H1A.oo-diag(diag(sys.fa_oo)),t3c,'mi,abcmjk->abcijk');

    D2 = -einsum_kg(H1B.oo-diag(diag(sys.fb_oo)),t3c,'mj,abcimk->abcijk');

    D3 = +einsum_kg(H1A.vv-diag(diag(sys.fa_vv)),t3c,'ae,ebcijk->abcijk');
   
    D4 = +einsum_kg(H1B.vv-diag(diag(sys.fb_vv)),t3c,'be,aecijk->abcijk');

    D5 = 0.5*einsum_kg(H2C.oooo,t3c,'mnjk,abcimn->abcijk');

    D6 = einsum_kg(H2B.oooo,t3c,'mnij,abcmnk->abcijk');

    D7 = 0.5*einsum_kg(H2C.vvvv,t3c,'bcef,aefijk->abcijk');

    D8 = einsum_kg(H2B.vvvv,t3c,'abef,efcijk->abcijk');

    D9 = einsum_kg(H2A.voov,t3c,'amie,ebcmjk->abcijk');

    D10 = einsum_kg(H2B.voov,t3d,'amie,ebcmjk->abcijk');

    D11 = einsum_kg(H2B.ovvo,t3b,'mbej,aecimk->abcijk');

    D12 = einsum_kg(H2C.voov,t3c,'bmje,aecimk->abcijk');

    D13 = -einsum_kg(H2B.ovov,t3c,'mbie,aecmjk->abcijk');

    D14 = -einsum_kg(H2B.vovo,t3c,'amej,ebcimk->abcijk');

    D_2_6_14 = D2 + D6 + D14;
    D_2_6_14 = D_2_6_14 - permute(D_2_6_14,[1,2,3,4,6,5]);

    D_4_8_13 = D4 + D8 + D13;
    D_4_8_13 = D_4_8_13 - permute(D_4_8_13,[1,3,2,4,5,6]);

    D_11_12 = D11 + D12;
    D_11_12 = D_11_12 - permute(D_11_12,[1,3,2,4,5,6]) - permute(D_11_12,[1,2,3,4,6,5]) + permute(D_11_12,[1,3,2,4,6,5]);

    X3C_abcijk = MM23C + D1 + D_2_6_14 + D3 + D_4_8_13 + D5 + D7 + D9 + D10 + D_11_12;

     t3c = zeros(size(t3c));
     for a = 1:sys.Nvir_alpha
        for b = 1:sys.Nvir_beta
            for c = b+1:sys.Nvir_beta
                for i = 1:sys.Nocc_alpha
                    for j = 1:sys.Nocc_beta
                        for k = j+1:sys.Nocc_beta

                            t3c(a,b,c,i,j,k) = X3C_abcijk(a,b,c,i,j,k)/...
                                (sys.fa_oo(i,i)+sys.fb_oo(j,j)+sys.fb_oo(k,k)-sys.fa_vv(a,a)-sys.fb_vv(b,b)-sys.fb_vv(c,c)-shift);
                            t3c(a,c,b,i,j,k) = -t3c(a,b,c,i,j,k);
                            t3c(a,b,c,i,k,j) = -t3c(a,b,c,i,j,k);
                            t3c(a,c,b,i,k,j) = t3c(a,b,c,i,j,k);

                        end
                    end
                end
            end
        end
    end


end

