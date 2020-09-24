function [MM23A] = build_MM23A(cc_t,HBar_t,sys)

    fprintf('MMCC(2,3)A construction... ')
    
    tic

    h_me = HBar_t.H1A.ov;
    h_amij = HBar_t.H2A.vooo;
    h_abie = HBar_t.H2A.vvov;
    t2a = cc_t.t2a;

    I12 = h_abie+einsum_kg(h_me,t2a,'me,abim->abie');
    
    % < phi_{ijkabc} | H_{CCSD} | 0 >
    % = -A(k/ij)A(a/bc) h(amij)*t(bcmk) + A(i/jk)A(c/ab)(h(abie)-h(me)*t(abim))*t(ecjk)
    MM23A =   -einsum_kg(h_amij,t2a,'amij,bcmk->abcijk') ... % (1)
              +einsum_kg(h_amij,t2a,'amkj,bcmi->abcijk') ... % (ik)
              +einsum_kg(h_amij,t2a,'amik,bcmj->abcijk') ... % (jk)
              +einsum_kg(h_amij,t2a,'cmij,bamk->abcijk') ... % (ac)
              +einsum_kg(h_amij,t2a,'bmij,acmk->abcijk') ... % (ab)
              -einsum_kg(h_amij,t2a,'bmkj,acmi->abcijk') ... % (ab)(ik)
              -einsum_kg(h_amij,t2a,'cmkj,bami->abcijk') ... % (ac)(ik)
              -einsum_kg(h_amij,t2a,'bmik,acmj->abcijk') ... % (ab)(jk)
              -einsum_kg(h_amij,t2a,'cmik,bamj->abcijk');    % (ac)(jk)
         
    MM23A =  MM23A +einsum_kg(I12,t2a,'abie,ecjk->abcijk') ... % (1)
                  -einsum_kg(I12,t2a,'abje,ecik->abcijk') ... % (ij)
                  -einsum_kg(I12,t2a,'abke,ecji->abcijk') ... % (ik)
                  -einsum_kg(I12,t2a,'cbie,eajk->abcijk') ... % (ac)
                  -einsum_kg(I12,t2a,'acie,ebjk->abcijk') ... % (bc)
                  +einsum_kg(I12,t2a,'cbje,eaik->abcijk') ... % (ac)(ij)
                  +einsum_kg(I12,t2a,'acje,ebik->abcijk') ... % (bc)(ij)
                  +einsum_kg(I12,t2a,'cbke,eaji->abcijk') ... % (ac)(ik)
                  +einsum_kg(I12,t2a,'acke,ebji->abcijk');    % (bc)(ik)

    fprintf('finished in %4.2f s\n',toc)

end

