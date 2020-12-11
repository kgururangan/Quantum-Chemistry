function [T3A,T3B,T3C,T3D] = make_t3_act_struct(cc_t,sys)

%     t1a = cc_t.t1a; t1b = cc_t.t1b;
%     t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
    t3a = cc_t.t3a; t3b = cc_t.t3b; t3c = cc_t.t3c; t3d = cc_t.t3d;

    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

    % t3A
    T3A.PPPHHH = t3a(PA,PA,PA,HA,HA,HA);
    T3A.PPPhHH = t3a(PA,PA,PA,hA,HA,HA);
    T3A.PPPhhH = t3a(PA,PA,PA,hA,hA,HA);
    T3A.PPPhhh = t3a(PA,PA,PA,hA,hA,hA);
    T3A.PPpHHH = t3a(PA,PA,pA,HA,HA,HA);
    T3A.PPphHH = t3a(PA,PA,pA,hA,HA,HA);
    T3A.PPphhH = t3a(PA,PA,pA,hA,hA,HA);
    T3A.PPphhh = t3a(PA,PA,pA,hA,hA,hA);
    T3A.PppHHH = t3a(PA,pA,pA,HA,HA,HA);
    T3A.PpphHH = t3a(PA,pA,pA,hA,HA,HA);
    T3A.PpphhH = t3a(PA,pA,pA,hA,hA,HA);
    T3A.Ppphhh = t3a(PA,pA,pA,hA,hA,hA);
    T3A.pppHHH = t3a(PA,pA,pA,HA,HA,HA);
    T3A.ppphHH = t3a(PA,pA,pA,hA,HA,HA);
    T3A.ppphhH = t3a(PA,pA,pA,hA,hA,HA);
    T3A.ppphhh = t3a(PA,pA,pA,hA,hA,hA);


    % t3B
    T3B.PPPHHH = t3b(PA,PA,PB,HA,HA,HB);
    T3B.PPPhHH = t3b(PA,PA,PB,hA,HA,HB);
    T3B.PPPHHh = t3b(PA,PA,PB,HA,HA,hB);
    T3B.PPPhHh = t3b(PA,PA,PB,hA,HA,hB);
    T3B.PPPhhH = t3b(PA,PA,PB,hA,hA,HB);
    T3B.PPPhhh = t3b(PA,PA,PB,hA,hA,hB);
    T3B.PpPHHH = t3b(PA,pA,PB,HA,HA,HB);
    T3B.PpPhHH = t3b(PA,pA,PB,hA,HA,HB);
    T3B.PpPHHh = t3b(PA,pA,PB,HA,HA,hB);
    T3B.PpPhHh = t3b(PA,pA,PB,hA,HA,hB);
    T3B.PpPhhH = t3b(PA,pA,PB,hA,hA,HB);
    T3B.PpPhhh = t3b(PA,pA,PB,hA,hA,hB);
    T3B.PpPHHh = t3b(PA,pA,PB,HA,HA,hB);
    T3B.PPpHHH = t3b(PA,PA,pB,HA,HA,HB);
    T3B.PPphHH = t3b(PA,PA,pB,hA,HA,HB);
    T3B.PPpHHh = t3b(PA,PA,pB,HA,HA,hB);
    T3B.PPphHh = t3b(PA,PA,pB,hA,HA,hB);
    T3B.PPphhH = t3b(PA,PA,pB,hA,hA,HB);
    T3B.PPphhh = t3b(PA,PA,pB,hA,hA,hB);
    T3B.PPpHHh = t3b(PA,PA,pB,HA,HA,hB);
    T3B.PppHHH = t3b(PA,pA,pB,HA,HA,HB);
    T3B.PpphHH = t3b(PA,pA,pB,hA,HA,HB);
    T3B.PppHHh = t3b(PA,pA,pB,HA,HA,hB);
    T3B.PpphHh = t3b(PA,pA,pB,hA,HA,hB);
    T3B.PpphhH = t3b(PA,pA,pB,hA,hA,HB);
    T3B.Ppphhh = t3b(PA,pA,pB,hA,hA,hB);
    T3B.PppHHh = t3b(PA,pA,pB,HA,HA,hB);
    T3B.ppPHHH = t3b(pA,pA,PB,HA,HA,HB);
    T3B.ppPhHH = t3b(pA,pA,PB,hA,HA,HB);
    T3B.ppPHHh = t3b(pA,pA,PB,HA,HA,hB);
    T3B.ppPhHh = t3b(pA,pA,PB,hA,HA,hB);
    T3B.ppPhhH = t3b(pA,pA,PB,hA,hA,HB);
    T3B.ppPhhh = t3b(pA,pA,PB,hA,hA,hB);
    T3B.ppPHHh = t3b(pA,pA,PB,HA,HA,hB);   
    T3B.pppHHH = t3b(pA,pA,pB,HA,HA,HB);
    T3B.ppphHH = t3b(pA,pA,pB,hA,HA,HB);
    T3B.pppHHh = t3b(pA,pA,pB,HA,HA,hB);
    T3B.ppphHh = t3b(pA,pA,pB,hA,HA,hB);
    T3B.ppphhH = t3b(pA,pA,pB,hA,hA,HB);
    T3B.ppphhh = t3b(pA,pA,pB,hA,hA,hB);
    T3B.pppHHh = t3b(pA,pA,pB,HA,HA,hB);

    % t3C
    T3C.PPPHHH = t3c(PA,PB,PB,HA,HB,HB);
    T3C.PPPhHH = t3c(PA,PB,PB,hA,HB,HB);
    T3C.PPPHhH = t3c(PA,PB,PB,HA,hB,HB);
    T3C.PPPhhH = t3c(PA,PB,PB,hA,hB,HB);
    T3C.PPPHhh = t3c(PA,PB,PB,HA,hB,hB);
    T3C.PPPhhh = t3c(PA,PB,PB,hA,hB,hB); 
    T3C.pPPHHH = t3c(pA,PB,PB,HA,HB,HB);
    T3C.pPPhHH = t3c(pA,PB,PB,hA,HB,HB);
    T3C.pPPHhH = t3c(pA,PB,PB,HA,hB,HB);
    T3C.pPPhhH = t3c(pA,PB,PB,hA,hB,HB);
    T3C.pPPHhh = t3c(pA,PB,PB,HA,hB,hB);
    T3C.pPPhhh = t3c(pA,PB,PB,hA,hB,hB);
    T3C.PPpHHH = t3c(PA,PB,pB,HA,HB,HB);
    T3C.PPphHH = t3c(PA,PB,pB,hA,HB,HB);
    T3C.PPpHhH = t3c(PA,PB,pB,HA,hB,HB);
    T3C.PPphhH = t3c(PA,PB,pB,hA,hB,HB);
    T3C.PPpHhh = t3c(PA,PB,pB,HA,hB,hB);
    T3C.PPphhh = t3c(PA,PB,pB,hA,hB,hB);
    T3C.pPpHHH = t3c(pA,PB,pB,HA,HB,HB);
    T3C.pPphHH = t3c(pA,PB,pB,hA,HB,HB);
    T3C.pPpHhH = t3c(pA,PB,pB,HA,hB,HB);
    T3C.pPphhH = t3c(pA,PB,pB,hA,hB,HB);
    T3C.pPpHhh = t3c(pA,PB,pB,HA,hB,hB);
    T3C.pPphhh = t3c(pA,PB,pB,hA,hB,hB);
    T3C.PppHHH = t3c(PA,pB,pB,HA,HB,HB);
    T3C.PpphHH = t3c(PA,pB,pB,hA,HB,HB);
    T3C.PppHhH = t3c(PA,pB,pB,HA,hB,HB);
    T3C.PpphhH = t3c(PA,pB,pB,hA,hB,HB);
    T3C.PppHhh = t3c(PA,pB,pB,HA,hB,hB);
    T3C.Ppphhh = t3c(PA,pB,pB,hA,hB,hB);
    T3C.pppHHH = t3c(pA,pB,pB,HA,HB,HB);
    T3C.ppphHH = t3c(pA,pB,pB,hA,HB,HB);
    T3C.pppHhH = t3c(pA,pB,pB,HA,hB,HB);
    T3C.ppphhH = t3c(pA,pB,pB,hA,hB,HB);
    T3C.pppHhh = t3c(pA,pB,pB,HA,hB,hB);
    T3C.ppphhh = t3c(pA,pB,pB,hA,hB,hB);
    

    % t3D
    T3D.PPPHHH = t3d(PB,PB,PB,HB,HB,HB);
    T3D.PPPhHH = t3d(PB,PB,PB,hB,HB,HB);
    T3D.PPPhhH = t3d(PB,PB,PB,hB,hB,HB);
    T3D.PPPhhh = t3d(PB,PB,PB,hB,hB,hB);
    T3D.PPpHHH = t3d(PB,PB,pB,HB,HB,HB);
    T3D.PPphHH = t3d(PB,PB,pB,hB,HB,HB);
    T3D.PPphhH = t3d(PB,PB,pB,hB,hB,HB);
    T3D.PPphhh = t3d(PB,PB,pB,hB,hB,hB);
    T3D.PppHHH = t3d(PB,pB,pB,HB,HB,HB);
    T3D.PpphHH = t3d(PB,pB,pB,hB,HB,HB);
    T3D.PpphhH = t3d(PB,pB,pB,hB,hB,HB);
    T3D.Ppphhh = t3d(PB,pB,pB,hB,hB,hB);
    T3D.pppHHH = t3d(PB,pB,pB,HB,HB,HB);
    T3D.ppphHH = t3d(PB,pB,pB,hB,HB,HB);
    T3D.ppphhH = t3d(PB,pB,pB,hB,hB,HB);
    T3D.ppphhh = t3d(PB,pB,pB,hB,hB,hB);
    
end