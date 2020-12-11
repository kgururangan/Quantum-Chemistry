function [MM23B] = build_MM23B(cc_t,HBar_t,sys)

    fprintf('MMCC(2,3)B construction... ')
    
    tic
    
    h_bcek = HBar_t.H2B.vvvo; 
    h_mcjk = HBar_t.H2B.ovoo;
    h_acie = HBar_t.H2B.vvov;
    h_amik = HBar_t.H2B.vooo;
    h_amij = HBar_t.H2A.vooo;
    h_abie = HBar_t.H2A.vvov;
    
    h1A_me = HBar_t.H1A.ov;
    h1B_me = HBar_t.H1B.ov;
    
    t2a = cc_t.t2a; t2b = cc_t.t2b;
    
    I1 = h_mcjk - einsum_kg(h1A_me,t2b,'me,ecjk->mcjk');
    I2 = h_amik - einsum_kg(h1B_me,t2b,'me,aeik->amik');
    I3 = h_amij - einsum_kg(h1A_me,t2a,'me,aeij->amij');

%     I1 = h_mcjk + einsum_kg(h1A_me,t2b,'me,ecjk->mcjk');
%     I2 = h_amik + einsum_kg(h1B_me,t2b,'me,aeik->amik');
%     I3 = h_amij + einsum_kg(h1A_me,t2a,'me,aeij->amij');    
   
    MM23B = einsum_kg(h_bcek,t2a,'bcek,aeij->abcijk')...
            - einsum_kg(h_bcek,t2a,'acek,beij->abcijk')...
            - einsum_kg(I1,t2a,'mcjk,abim->abcijk')...
            + einsum_kg(I1,t2a,'mcik,abjm->abcijk');
        
    MM23B = MM23B ...
            + einsum_kg(h_acie,t2b,'acie,bejk->abcijk')...
            - einsum_kg(h_acie,t2b,'bcie,aejk->abcijk')...
            - einsum_kg(h_acie,t2b,'acje,beik->abcijk')...
            + einsum_kg(h_acie,t2b,'bcje,aeik->abcijk')...
            - einsum_kg(I2,t2b,'amik,bcjm->abcijk')...
            + einsum_kg(I2,t2b,'bmik,acjm->abcijk')...
            + einsum_kg(I2,t2b,'amjk,bcim->abcijk')...
            - einsum_kg(I2,t2b,'bmjk,acim->abcijk');
        
    MM23B = MM23B ...
            + einsum_kg(h_abie,t2b,'abie,ecjk->abcijk')...
            - einsum_kg(h_abie,t2b,'abje,ecik->abcijk')...
            - einsum_kg(I3,t2b,'amij,bcmk->abcijk')...
            + einsum_kg(I3,t2b,'bmij,acmk->abcijk');
        
    fprintf('finished in %4.2f s\n',toc)

end

