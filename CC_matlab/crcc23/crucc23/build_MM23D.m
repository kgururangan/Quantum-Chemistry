function [MM23D] = build_MM23D(cc_t,HBar_t,sys)

    fprintf('MMCC(2,3)D construction... ')
    
    tic

    h_me = HBar_t.H1B.ov;
    h_amij = HBar_t.H2C.vooo;
    h_abie = HBar_t.H2C.vvov;
    t2c = cc_t.t2c;

    I12 = h_abie-einsum_kg(h_me,t2c,'me,abim->abie');
    
    % < phi_{ijkabc} | H_{CCSD} | 0 >
    % = -A(k/ij)A(a/bc) h(amij)*t(bcmk) + A(i/jk)A(c/ab)(h(abie)-h(me)*t(abim))*t(ecjk)
    MM23D =   -einsum_kg(h_amij,t2c,'amij,bcmk->abcijk') ... % (1)
             +einsum_kg(h_amij,t2c,'amkj,bcmi->abcijk') ... % (ik)
             +einsum_kg(h_amij,t2c,'amik,bcmj->abcijk') ... % (jk)
             +einsum_kg(h_amij,t2c,'cmij,bamk->abcijk') ... % (ac)
             +einsum_kg(h_amij,t2c,'bmij,acmk->abcijk') ... % (ab)
             -einsum_kg(h_amij,t2c,'bmkj,acmi->abcijk') ... % (ab)(ik)
             -einsum_kg(h_amij,t2c,'cmkj,bami->abcijk') ... % (ac)(ik)
             -einsum_kg(h_amij,t2c,'bmik,acmj->abcijk') ... % (ab)(jk)
             -einsum_kg(h_amij,t2c,'cmik,bamj->abcijk');    % (ac)(jk)
         
    MM23D =  MM23D +einsum_kg(I12,t2c,'abie,ecjk->abcijk') ... % (1)
                  -einsum_kg(I12,t2c,'abje,ecik->abcijk') ... % (ij)
                  -einsum_kg(I12,t2c,'abke,ecji->abcijk') ... % (ik)
                  -einsum_kg(I12,t2c,'cbie,eajk->abcijk') ... % (ac)
                  -einsum_kg(I12,t2c,'acie,ebjk->abcijk') ... % (bc)
                  +einsum_kg(I12,t2c,'cbje,eaik->abcijk') ... % (ac)(ij)
                  +einsum_kg(I12,t2c,'acje,ebik->abcijk') ... % (bc)(ij)
                  +einsum_kg(I12,t2c,'cbke,eaji->abcijk') ... % (ac)(ik)
                  +einsum_kg(I12,t2c,'acke,ebji->abcijk');    % (bc)(ik)

    fprintf(' finished in %4.2f s\n',toc)

end