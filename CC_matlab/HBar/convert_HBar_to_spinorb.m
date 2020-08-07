function [HBar] = convert_HBar_to_spinorb(HBar_t)

    % The only spin-cases for 1-body are
    % A - (a,a)
    % B - (b,b)
    
    % The only spin-cases for 2-body are
    % A - (a,a,a,a)
    % B - (a,b,a,b), (a,b,b,a), (b,a,a,b), (b,a,b,a)
    % C - (b,b,b,b)
    
    % The only spin-cases for 3-body are
    % A - ((a,a,a,a,a,a), 
    % B - (a,a,b,a,a,b), (a,b,a,a,a,b), (b,a,a,a,a,b)
    %     (a,a,b,a,b,a), (a,b,a,a,b,a), (b,a,a,a,b,a)
    %     (a,a,b,a,a,b), (a,b,a,b,a,a), (b,a,a,b,a,a)
    % C - (a,b,b,a,b,b), (b,a,b,a,b,b), (b,b,a,a,b,b)
    %     (a,b,b,b,a,b), (b,a,b,b,a,b), (b,b,a,b,a,b)
    %     (a,b,b,b,b,a), (b,a,b,b,b,a), (b,b,a,b,b,a)
    % D - (b,b,b,b,b,b)
    
    % Note: we are not allowed to have spin cases that are like H(a,a,b,b)
    % or H(b,b,a,a).
    
    % E.g.
    % If you need to get Hvooo(ub,oa,ob,oa), you need to permute Hovoo since
    % Hvooo can't produce ub. The permutation is permute(H2B.ovoo,[2,1,4,3]) 
    % which is 2 swaps so a sign of +1. The sign must be added because 
    % H2B.ovoo is not intrisincally antisymmetric (but A and C are).
    
    fprintf('\n======================================================\n')
    fprintf('Converting spin-integrated HBar into spinorbital HBar')
    fprintf('\n======================================================\n')

    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
    H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3C = HBar_t.H3C; H3D = HBar_t.H3D;
    
    H1 = cell(2,2);
    H2 = cell(2,2,2,2);
    H3 = cell(2,2,2,2,2,2);
    
    [Nocc_a, Nunocc_a] = size(H1A.ov); [Nocc_b, Nunocc_b] = size(H1B.ov);
    Nocc = Nocc_a + Nocc_b; 
    Nunocc = Nunocc_a + Nunocc_b;
    
    Norb = Nunocc + Nocc;
    occ = 1:Nocc; unocc = [Nocc+1:Norb] - Nocc;
    oa = occ(1:2:end); ob = occ(2:2:end);
    ua = unocc(1:2:end); ub = unocc(2:2:end);


    % H(mi)
    % spinorbital container
    Hchk_oo = zeros(Nocc,Nocc);
    Hchk_oo(oa,oa) = H1A.oo;
    Hchk_oo(ob,ob) = H1B.oo;
    H1{1,1} = Hchk_oo;
    
    % H(ae)
    % spinorbital container
    Hchk_vv = zeros(Nunocc,Nunocc);
    Hchk_vv(ua,ua) = H1A.vv;
    Hchk_vv(ub,ub) = H1B.vv;
    H1{2,2} = Hchk_vv;

    % H(me)
    % spinorbital container
    Hchk_ov = zeros(Nocc,Nunocc);
    Hchk_ov(oa,ua) = H1A.ov;
    Hchk_ov(ob,ub) = H1B.ov;
    H1{1,2} = Hchk_ov;

    % H(amij)
    % spinorbital container
    Hchk_vooo = zeros(Nunocc,Nocc,Nocc,Nocc);
    Hchk_vooo(ua,oa,oa,oa) = H2A.vooo;
    Hchk_vooo(ub,ob,ob,ob) = H2C.vooo;
    Hchk_vooo(ua,ob,oa,ob) = H2B.vooo;
    Hchk_vooo(ub,oa,ob,oa) = permute(H2B.ovoo,[2,1,4,3]);
    Hchk_vooo(ub,oa,oa,ob) = -permute(H2B.ovoo,[2,1,3,4]);
    Hchk_vooo(ua,ob,ob,oa) = -permute(H2B.vooo,[1,2,4,3]);
    H2{2,1,1,1} = Hchk_vooo;
    H2{1,2,1,1} = -permute(Hchk_vooo,[2,1,3,4]);

    % H(abie)
    % spinorbital container
    % (ua,ua,oa,ua), (ub,ub,ob,ub), (ua,ub,oa,ub),
    % (ub,ua,oa,ub), (ua,ub,ob,ua), (ub,ua,ob,ua)
    Hchk_vvov = zeros(Nunocc,Nunocc,Nocc,Nunocc);
    Hchk_vvov(ua,ua,oa,ua) = H2A.vvov;
    Hchk_vvov(ub,ub,ob,ub) = H2C.vvov;
    Hchk_vvov(ua,ub,oa,ub) = H2B.vvov;
    Hchk_vvov(ub,ua,oa,ub) = -permute(H2B.vvov,[2,1,3,4]);
    Hchk_vvov(ua,ub,ob,ua) = -permute(H2B.vvvo,[1,2,4,3]);
    Hchk_vvov(ub,ua,ob,ua) = permute(H2B.vvvo,[2,1,4,3]);
    H2{2,2,1,2} = Hchk_vvov;
    H2{2,2,2,1} = -permute(Hchk_vvov,[1,2,4,3]);


    % H(mnij)
    % spinorbital container
    Hchk_oooo = zeros(Nocc,Nocc,Nocc,Nocc);
    Hchk_oooo(oa,oa,oa,oa) = H2A.oooo;
    Hchk_oooo(oa,ob,oa,ob) = H2B.oooo;
    Hchk_oooo(ob,oa,oa,ob) = -permute(H2B.oooo,[2,1,3,4]);
    Hchk_oooo(oa,ob,ob,oa) = -permute(H2B.oooo,[1,2,4,3]);
    Hchk_oooo(ob,oa,ob,oa) = permute(H2B.oooo,[2,1,4,3]);
    Hchk_oooo(ob,ob,ob,ob) = H2C.oooo;
    H2{1,1,1,1} = Hchk_oooo;


    % H(abef)
    % spinorbital container
    Hchk_vvvv = zeros(Nunocc,Nunocc,Nunocc,Nunocc);
    Hchk_vvvv(ua,ua,ua,ua) = H2A.vvvv;
    Hchk_vvvv(ua,ub,ua,ub) = H2B.vvvv;
    Hchk_vvvv(ub,ua,ua,ub) = -permute(H2B.vvvv,[2,1,3,4]);
    Hchk_vvvv(ua,ub,ub,ua) = -permute(H2B.vvvv,[1,2,4,3]);
    Hchk_vvvv(ub,ua,ub,ua) = permute(H2B.vvvv,[2,1,4,3]);
    Hchk_vvvv(ub,ub,ub,ub) = H2C.vvvv;
    H2{2,2,2,2} = Hchk_vvvv;

    % H(amie)
    % spinorbital container
    % (ua,oa,oa,ua), (ub,ob,ob,ub), (ua,ob,oa,ub),
    % (ub,oa,ob,ua), (ua,ob,ob,ua), (ub,oa,oa,ub),
    Hchk_voov = zeros(Nunocc,Nocc,Nocc,Nunocc);
    Hchk_voov(ua,oa,oa,ua) = H2A.voov;
    Hchk_voov(ua,ob,oa,ub) = H2B.voov;
    Hchk_voov(ub,oa,ob,ua) = permute(H2B.ovvo,[2,1,4,3]);
    Hchk_voov(ua,ob,ob,ua) = -permute(H2B.vovo,[1,2,4,3]);
    Hchk_voov(ub,oa,oa,ub) = -permute(H2B.ovov,[2,1,3,4]);
    Hchk_voov(ub,ob,ob,ub) = H2C.voov;
    H2{2,1,1,2} = Hchk_voov;
    H2{1,2,1,2} = -permute(Hchk_voov,[2,1,3,4]);
    H2{2,1,2,1} = -permute(Hchk_voov,[1,2,4,3]);
    H2{1,2,2,1} = permute(Hchk_voov,[2,1,4,3]);


    % H(amef)
    % (ua,oa,ua,ua), (ub,ob,ub,ub), (ua,ob,ua,ub),
    % (ua,ob,ub,ua), (ub,oa,ua,ub), (ub,oa,ub,ua)
    % spinorbital container
    Hchk_vovv = zeros(Nunocc,Nocc,Nunocc,Nunocc);
    Hchk_vovv(ua,oa,ua,ua) = H2A.vovv;
    Hchk_vovv(ub,ob,ub,ub) = H2C.vovv;
    Hchk_vovv(ua,ob,ua,ub) = H2B.vovv;
    Hchk_vovv(ua,ob,ub,ua) = -permute(H2B.vovv,[1,2,4,3]);
    Hchk_vovv(ub,oa,ua,ub) = -permute(H2B.ovvv,[2,1,3,4]);
    Hchk_vovv(ub,oa,ub,ua) = permute(H2B.ovvv,[2,1,4,3]);
    H2{2,1,2,2} = Hchk_vovv;
    H2{1,2,2,2} = -permute(Hchk_vovv,[2,1,3,4]);

    % H(mnie)
    % spinorbital container
    % (oa,oa,oa,ua), (ob,ob,ob,ub), (oa,ob,oa,ub),
    % (ob,oa,oa,ub), (oa,ob,ob,ua), (ob,oa,ob,ua)
    Hchk_ooov = zeros(Nocc,Nocc,Nocc,Nunocc);
    Hchk_ooov(oa,oa,oa,ua) = H2A.ooov;
    Hchk_ooov(ob,ob,ob,ub) = H2C.ooov;
    Hchk_ooov(oa,ob,oa,ub) = H2B.ooov;
    Hchk_ooov(ob,oa,oa,ub) = -permute(H2B.ooov,[2,1,3,4]);
    Hchk_ooov(oa,ob,ob,ua) = -permute(H2B.oovo,[1,2,4,3]);
    Hchk_ooov(ob,oa,ob,ua) = permute(H2B.oovo,[2,1,4,3]);
    H2{1,1,1,2} = Hchk_ooov;
    H2{1,1,2,1} = -permute(Hchk_ooov,[1,2,4,3]);

    % H(mnef)
    % spinorbital container
    Hchk_oovv = zeros(Nocc,Nocc,Nunocc,Nunocc);
    Hchk_oovv(oa,oa,ua,ua) = H2A.oovv;
    Hchk_oovv(ob,ob,ub,ub) = H2C.oovv;
    Hchk_oovv(oa,ob,ua,ub) = H2B.oovv;
    Hchk_oovv(ob,oa,ua,ub) = -permute(H2B.oovv,[2,1,3,4]);
    Hchk_oovv(oa,ob,ub,ua) = -permute(H2B.oovv,[1,2,4,3]);
    Hchk_oovv(ob,oa,ub,ua) = permute(H2B.oovv,[2,1,4,3]);
    H2{1,1,2,2} = Hchk_oovv;

    % H(abmije)
    % H3A: (ua,ua,oa,oa,oa,ua)
    
    % H3D: (ub,ub,ob,ob,ob,ub)
    
    % H3B: (ua,ua,ob,oa,oa,ub), (ub,ua,oa,oa,oa,ub), (ua,ub,oa,oa,oa,ub)
    %      (ua,ua,ob,oa,ob,ua), (ub,ua,oa,oa,ob,ua), (ua,ub,oa,oa,ob,ua)
    %      (ua,ua,ob,ob,oa,ua), (ub,ua,oa,ob,oa,ua), (ua,ub,oa,ob,oa,ua)
    
    % H3C: (ua,ub,ob,oa,ob,ub), (ub,ua,ob,oa,ob,ub), (ub,ub,oa,oa,ob,ub)
    %      (ua,ub,ob,ob,oa,ub), (ub,ua,ob,ob,oa,ub), (ub,ub,oa,ob,oa,ub)
    %      (ua,ub,ob,ob,ob,ua), (ub,ua,ob,ob,ob,ua), (ub,ub,oa,ob,ob,ua)

    % spinorbital container
    Hchk_vvooov = zeros(Nunocc,Nunocc,Nocc,Nocc,Nocc,Nunocc);

    Hchk_vvooov(ua,ua,oa,oa,oa,ua) = H3A.vvooov;
    
    Hchk_vvooov(ua,ua,ob,oa,oa,ub) = H3B.vvooov; 
    Hchk_vvooov(ua,ub,oa,oa,oa,ub) = -permute(H3B.vovoov,[1,3,2,4,5,6]);
    Hchk_vvooov(ub,ua,oa,oa,oa,ub) = permute(H3B.vovoov,[3,1,2,4,5,6]);
    
    Hchk_vvooov(ua,ua,ob,oa,ob,ua) = -permute(H3B.vvoovo,[1,2,3,4,6,5]);
    Hchk_vvooov(ua,ub,oa,oa,ob,ua) = permute(H3B.vovovo,[1,3,2,4,6,5]);
    Hchk_vvooov(ub,ua,oa,oa,ob,ua) = -permute(H3B.vovovo,[3,1,2,4,6,5]);
    
    Hchk_vvooov(ua,ua,ob,ob,oa,ua) = permute(H3B.vvoovo,[1,2,3,6,4,5]);
    Hchk_vvooov(ua,ub,oa,ob,oa,ua) = -permute(H3B.vovovo,[1,3,2,6,4,5]);
    Hchk_vvooov(ub,ua,oa,ob,oa,ua) = permute(H3B.vovovo,[3,1,2,6,4,5]);
    
    Hchk_vvooov(ua,ub,ob,oa,ob,ub) = H3C.vvooov;
    Hchk_vvooov(ub,ua,ob,oa,ob,ub) = -permute(H3C.vvooov,[2,1,3,4,5,6]);
    Hchk_vvooov(ub,ub,oa,oa,ob,ub) = -permute(H3C.ovvovo,[2,3,1,4,6,5]);
    
    Hchk_vvooov(ua,ub,ob,ob,oa,ub) = -permute(H3C.vvooov,[1,2,3,5,4,6]);
    Hchk_vvooov(ub,ua,ob,ob,oa,ub) = permute(H3C.vovovo,[3,1,2,6,4,5]);
    Hchk_vvooov(ub,ub,oa,ob,oa,ub) = permute(H3C.ovvovo,[2,3,1,6,4,5]);
    
    Hchk_vvooov(ua,ub,ob,ob,ob,ua) = permute(H3C.vvovoo,[1,2,3,5,6,4]);
    Hchk_vvooov(ub,ua,ob,ob,ob,ua) = -permute(H3C.vvovoo,[2,1,3,5,6,4]);
    Hchk_vvooov(ub,ub,oa,ob,ob,ua) = permute(H3C.ovvvoo,[2,3,1,5,6,4]);

    Hchk_vvooov(ub,ub,ob,ob,ob,ub) = H3D.vvooov;
    
    H3{2,2,1,1,1,2} = Hchk_vvooov;

 
    % H(abmief)
    % spinorbital container
    Hchk_vvoovv = zeros(Nunocc,Nunocc,Nocc,Nocc,Nunocc,Nunocc);
    
    Hchk_vvoovv(ua,ua,oa,oa,ua,ua) = H3A.vvoovv;
    
    Hchk_vvoovv(ua,ua,ob,oa,ua,ub) = H3B.vvoovv;
    Hchk_vvoovv(ua,ub,oa,oa,ua,ub) = -permute(H3B.vovovv,[1,3,2,4,5,6]);
    Hchk_vvoovv(ub,ua,oa,oa,ua,ub) = permute(H3B.vovovv,[3,1,2,4,5,6]);
    
    Hchk_vvoovv(ua,ua,ob,oa,ub,ua) = -permute(H3B.vvoovv,[1,2,3,4,6,5]);
    Hchk_vvoovv(ua,ub,oa,oa,ub,ua) = permute(H3B.vovovv,[1,3,2,4,6,5]);
    Hchk_vvoovv(ub,ua,oa,oa,ub,ua) = -permute(H3B.vovovv,[3,1,2,4,6,5]);
     
    Hchk_vvoovv(ua,ua,ob,ob,ua,ua) = permute(H3B.vvovvo,[1,2,3,6,4,5]);
    Hchk_vvoovv(ua,ub,oa,ob,ua,ua) = -permute(H3B.vovvvo,[1,3,2,6,4,5]);
    Hchk_vvoovv(ub,ua,oa,ob,ua,ua) = permute(H3B.vovvvo,[3,1,2,6,4,5]);
    
    Hchk_vvoovv(ua,ub,ob,oa,ub,ub) = H3C.vvoovv;
    Hchk_vvoovv(ub,ua,ob,oa,ub,ub) = -permute(H3C.vvoovv,[2,1,3,4,5,6]);
    Hchk_vvoovv(ub,ub,oa,oa,ub,ub) = permute(H3C.ovvovv,[2,3,1,4,5,6]);
   
    Hchk_vvoovv(ua,ub,ob,ob,ua,ub) = -permute(H3C.vvovov,[1,2,3,5,4,6]);
    Hchk_vvoovv(ub,ua,ob,ob,ua,ub) = permute(H3C.vvovov,[2,1,3,5,4,6]);
    Hchk_vvoovv(ub,ub,oa,ob,ua,ub) = -permute(H3C.ovvvov,[2,3,1,5,4,6]);
    
    Hchk_vvoovv(ua,ub,ob,ob,ub,ua) = permute(H3C.vvovov,[1,2,3,5,6,4]);
    Hchk_vvoovv(ub,ua,ob,ob,ub,ua) = -permute(H3C.vvovov,[2,1,3,5,6,4]);
    Hchk_vvoovv(ub,ub,oa,ob,ub,ua) = permute(H3C.ovvvov,[2,3,1,5,6,4]);
     
    Hchk_vvoovv(ub,ub,ob,ob,ub,ub) = H3D.vvoovv;
     
    H3{2,2,1,1,2,2} = Hchk_vvoovv;
    H3{2,1,2,1,2,2} = -permute(Hchk_vvoovv,[1,3,2,4,5,6]);
    H3{2,2,1,2,1,2} = -permute(Hchk_vvoovv,[1,2,3,5,4,6]);
    H3{2,1,2,2,2,1} = -permute(Hchk_vvoovv,[1,3,2,5,6,4]);
     

%     % H(amnije)
%     
%     Hchk_voooov(ua,oa,oa,oa,oa,ua) = H3A.voooov;
%     
%     Hchk_voooov(ua,oa,ob,oa,oa,ub) = H3B.voooov;
%     Hchk_voooov(ua,ob,oa,oa,oa,ub) = -permute(H3B.voooov,[1,3,2,4,5,6]);
%     Hchk_voooov(ub,oa,oa,oa,oa,ub) = permute(H3B.oovoov,[3,1,2,6,4,5]);
%     
%     Hchk_voooov(ua,oa,ob,oa,ob,ua) = permute(H3B.ovoovo
%     Hchk_voooov(ua,ob,oa,oa,ob,ua)
%     Hchk_voooov(ub,oa,oa,oa,ob,ua)
%     
%     Hchk_voooov(ua,oa,ob,ob,oa,ua)
%     Hchk_voooov(ua,ob,oa,ob,oa,ua)
%     Hchk_voooov(ub,oa,oa,ob,oa,ua)
%     



    HBar = {H1, H2, H3};

end