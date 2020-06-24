function [MM23] = build_MM23(t1,t2,HBar)
    
    tic

    h_me = HBar{1}{1,2};
    h_amij = HBar{2}{2,1,1,1};
    h_abie = HBar{2}{2,2,1,2};
    
    I12 = h_abie-einsum_kg(h_me,t2,'me,abim->abie');
    
    % < phi_{ijkabc} | H_{CCSD} | 0 >
    % = -A(k/ij)A(a/bc) h(amij)*t(bcmk) + A(i/jk)A(c/ab)(h(abie)-h(me)*t(abim))*t(ecjk)
    MM23 =   -einsum_kg(h_amij,t2,'amij,bcmk->abcijk') ... % (1)
             +einsum_kg(h_amij,t2,'amkj,bcmi->abcijk') ... % (ik)
             +einsum_kg(h_amij,t2,'amik,bcmj->abcijk') ... % (jk)
             +einsum_kg(h_amij,t2,'cmij,bamk->abcijk') ... % (ac)
             +einsum_kg(h_amij,t2,'bmij,acmk->abcijk') ... % (ab)
             -einsum_kg(h_amij,t2,'bmkj,acmi->abcijk') ... % (ab)(ik)
             -einsum_kg(h_amij,t2,'cmkj,bami->abcijk') ... % (ac)(ik)
             -einsum_kg(h_amij,t2,'bmik,acmj->abcijk') ... % (ab)(jk)
             -einsum_kg(h_amij,t2,'cmik,bamj->abcijk');    % (ac)(jk)
         
    MM23 =  MM23 +einsum_kg(I12,t2,'abie,ecjk->abcijk') ... % (1)
                 -einsum_kg(I12,t2,'abje,ecik->abcijk') ... % (ij)
                 -einsum_kg(I12,t2,'abke,ecji->abcijk') ... % (ik)
                 -einsum_kg(I12,t2,'cbie,eajk->abcijk') ... % (ac)
                 -einsum_kg(I12,t2,'acie,ebjk->abcijk') ... % (bc)
                 +einsum_kg(I12,t2,'cbje,eaik->abcijk') ... % (ac)(ij)
                 +einsum_kg(I12,t2,'acje,ebik->abcijk') ... % (bc)(ij)
                 +einsum_kg(I12,t2,'cbke,eaji->abcijk') ... % (ac)(ik)
                 +einsum_kg(I12,t2,'acke,ebji->abcijk');    % (bc)(ik)

    fprintf('MM(2,3) constructed in %4.2f s\n',toc)
end

