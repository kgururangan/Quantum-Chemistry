function [t3] = make_act_struct_into_t3(T_struct,spin_case,sys)

    % Problem:
    % I think the error lies in how we are making the active t3 tensor out
    % of the unique active space components. Since we want all tensor
    % components of t3 such that at least 2 occupied/unocciped
    % indices lie in the active space, we would have, for instance, the
    % element t3a(A,B,C,I,J,k) in the full t3a tensor. This is related to 
    % the unique entry t3A(A,B,C,k,I,J) by symmetry (sign +1), however, by 
    % naïvely setting t3a(PA,PA,PA,hA,HA,HA) = T3A.PPPhHH, we are neglecting
    % elements such as this that are related by permutational symmetry of
    % t3. 
    %
    % How to fix:
    % (1) Change the t1 and t2 updates to the active-space variants such
    % that only the unique active-defined T3A, T3B, T3C, and T3D entries
    % are used to construct the full t1/t2 updates
    % (2) use a full loop. E.g.
    % for a = 1:Nact_p_alpha
    %   for b = a+1:Nact_p_alpha
    %       for c = b+1:Nact_p_alpha
    %           for i = 1:Nunact_h_alpha
    %               for j = i+1:Nact_h_alpha
    %                   for k = j+1:Nact_h_alpha
    %                       t3a(PA(a),PA(b),PA(c),hA(i),HA(j),HA(k)) = T3A.PPPhHH(a,b,c,i,j,k);
    %                       t3a(PA(a),PA(b),PA(c),HA(j),hA(i),HA(k)) = -T3A.PPPhHH(a,b,c,i,j,k);
    %                       ... and so on

    Noa = sys.Nocc_alpha; Nob = sys.Nocc_beta;
    Nua = sys.Nvir_alpha; Nub = sys.Nvir_beta;

    switch spin_case

            case {'A','a'}

                    t3 = zeros(Nua,Nua,Nua,Noa,Noa,Noa);
%                     t3 = t3_in;

                    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;

                    t3(PA,PA,PA,HA,HA,HA) = T_struct.PPPHHH;

                    t3(PA,PA,PA,hA,HA,HA) = T_struct.PPPhHH;
                    t3(PA,PA,PA,HA,hA,HA) = -permute(T_struct.PPPhHH,[1,2,3,5,4,6]);
                    t3(PA,PA,PA,HA,HA,hA) = permute(T_struct.PPPhHH,[1,2,3,5,6,4]);
    
                    %t3(PA,PA,PA,hA,hA,HA) = T_struct.PPPhhH;
                    %t3(PA,PA,PA,hA,hA,hA) = T_struct.PPPhhh;

                    t3(PA,PA,pA,HA,HA,HA) = T_struct.PPpHHH;
                    t3(PA,pA,PA,HA,HA,HA) = -permute(T_struct.PPpHHH,[1,3,2,4,5,6]);
                    t3(pA,PA,PA,HA,HA,HA) = permute(T_struct.PPpHHH,[3,1,2,4,5,6]);

                    t3(PA,PA,pA,hA,HA,HA) = T_struct.PPphHH;
                    t3(PA,PA,pA,HA,hA,HA) = -permute(T_struct.PPphHH,[1,2,3,5,4,6]);
                    t3(PA,PA,pA,HA,HA,hA) = permute(T_struct.PPphHH,[1,2,3,5,6,4]);

                    t3(PA,pA,PA,hA,HA,HA) = -permute(T_struct.PPphHH,[1,3,2,4,5,6]);
                    t3(PA,pA,PA,HA,hA,HA) = permute(T_struct.PPphHH,[1,3,2,5,4,6]);
                    t3(PA,pA,PA,HA,HA,hA) = -permute(T_struct.PPphHH,[1,3,2,5,6,4]);

                    t3(pA,PA,PA,hA,HA,HA) = permute(T_struct.PPphHH,[3,1,2,4,5,6]);
                    t3(pA,PA,PA,HA,hA,HA) = -permute(T_struct.PPphHH,[3,1,2,5,4,6]);
                    t3(pA,PA,PA,HA,HA,hA) = permute(T_struct.PPphHH,[3,1,2,5,6,4]);




                    %t3(PA,PA,pA,hA,hA,HA) = T_struct.PPphhH;
                    %t3(PA,PA,pA,hA,hA,hA) = T_struct.PPphhh;

%                     t3(PA,pA,pA,HA,HA,HA) = T_struct.PppHHH;
%                     t3(PA,pA,pA,hA,HA,HA) = T_struct.PpphHH;
%                     t3(PA,pA,pA,hA,hA,HA) = T_struct.PpphhH;
%                     t3(PA,pA,pA,hA,hA,hA) = T_struct.Ppphhh;
% 
%                     t3(pA,pA,pA,HA,HA,HA) = T_struct.pppHHH;
%                     t3(pA,pA,pA,hA,HA,HA) = T_struct.ppphHH;
%                     t3(pA,pA,pA,hA,hA,HA) = T_struct.ppphhH;
%                     t3(pA,pA,pA,hA,hA,hA) = T_struct.ppphhh;

                
            case {'B','b'}

                    t3 = zeros(Nua,Nua,Nub,Noa,Noa,Nob);
%                     t3 = t3_in;

                    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
                    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

                    t3(PA,PA,PB,HA,HA,HB) = T_struct.PPPHHH;

                    t3(PA,PA,PB,hA,HA,HB) = T_struct.PPPhHH;
                    t3(PA,PA,PB,HA,hA,HB) = -permute(T_struct.PPPhHH,[1,2,3,5,4,6]);

                    t3(PA,PA,PB,HA,HA,hB) = T_struct.PPPHHh;

%                     t3(PA,PA,PB,hA,HA,hB) = T_struct.PPPhHh;
%                     t3(PA,PA,PB,hA,hA,HB) = T_struct.PPPhhH;
%                     t3(PA,PA,PB,hA,hA,hB) = T_struct.PPPhhh;

                    t3(PA,pA,PB,HA,HA,HB) = T_struct.PpPHHH;
                    t3(pA,PA,PB,HA,HA,HB) = -permute(T_struct.PpPHHH,[2,1,3,4,5,6]);

                    t3(PA,pA,PB,hA,HA,HB) = T_struct.PpPhHH;
                    t3(PA,pA,PB,HA,hA,HB) = -permute(T_struct.PpPhHH,[1,2,3,5,4,6]);
                    t3(pA,PA,PB,hA,HA,HB) = -permute(T_struct.PpPhHH,[2,1,3,4,5,6]);
                    t3(pA,PA,PB,HA,hA,HB) = permute(T_struct.PpPhHH,[2,1,3,5,4,6]);

                    t3(PA,pA,PB,HA,HA,hB) = T_struct.PpPHHh;
                    t3(pA,PA,PB,HA,HA,hB) = -permute(T_struct.PpPHHh,[2,1,3,4,5,6]);

%                     t3(PA,pA,PB,hA,HA,hB) = T_struct.PpPhHh;
%                     t3(PA,pA,PB,hA,hA,HB) = T_struct.PpPhhH;
%                     t3(PA,pA,PB,hA,hA,hB) = T_struct.PpPhhh;

                    t3(PA,PA,pB,HA,HA,HB) = T_struct.PPpHHH;

                    t3(PA,PA,pB,hA,HA,HB) = T_struct.PPphHH;
                    t3(PA,PA,pB,HA,hA,HB) = -permute(T_struct.PPphHH,[1,2,3,5,4,6]);

                    t3(PA,PA,pB,HA,HA,hB) = T_struct.PPpHHh;

%                     t3(PA,PA,pB,hA,HA,hB) = T_struct.PPphHh;
%                     t3(PA,PA,pB,hA,hA,HB) = T_struct.PPphhH;
%                     t3(PA,PA,pB,hA,hA,hB) = T_struct.PPphhh;

%                     t3(PA,pA,pB,HA,HA,HB) = T_struct.PppHHH;
%                     t3(PA,pA,pB,hA,HA,HB) = T_struct.PpphHH;
%                     t3(PA,pA,pB,HA,HA,hB) = T_struct.PppHHh;
%                     t3(PA,pA,pB,hA,HA,hB) = T_struct.PpphHh;
%                     t3(PA,pA,pB,hA,hA,HB) = T_struct.PpphhH;
%                     t3(PA,pA,pB,hA,hA,hB) = T_struct.Ppphhh;
% 
%                     t3(pA,pA,PB,HA,HA,HB) = T_struct.ppPHHH;
%                     t3(pA,pA,PB,hA,HA,HB) = T_struct.ppPhHH;
%                     t3(pA,pA,PB,HA,HA,hB) = T_struct.ppPHHh;
%                     t3(pA,pA,PB,hA,HA,hB) = T_struct.ppPhHh;
%                     t3(pA,pA,PB,hA,hA,HB) = T_struct.ppPhhH;
%                     t3(pA,pA,PB,hA,hA,hB) = T_struct.ppPhhh;
% 
%                     t3(pA,pA,pB,HA,HA,HB) = T_struct.pppHHH;
%                     t3(pA,pA,pB,hA,HA,HB) = T_struct.ppphHH;
%                     t3(pA,pA,pB,HA,HA,hB) = T_struct.pppHHh;
%                     t3(pA,pA,pB,hA,HA,hB) = T_struct.ppphHh;
%                     t3(pA,pA,pB,hA,hA,HB) = T_struct.ppphhH;
%                     t3(pA,pA,pB,hA,hA,hB) = T_struct.ppphhh;

                
            case {'C','c'}
    
                    t3 = zeros(Nua,Nub,Nub,Noa,Nob,Nob);
%                     t3 = t3_in;

                    PA = sys.PA; pA = sys.pA; HA = sys.HA; hA = sys.hA;
                    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

                    t3(PA,PB,PB,HA,HB,HB) = T_struct.PPPHHH;

                    t3(PA,PB,PB,hA,HB,HB) = T_struct.PPPhHH;

                    t3(PA,PB,PB,HA,hB,HB) = T_struct.PPPHhH;        
                    t3(PA,PB,PB,HA,HB,hB) = -permute(T_struct.PPPHhH,[1,2,3,4,6,5]);

%                     t3(PA,PB,PB,hA,hB,HB) = T_struct.PPPhhH;
%                     t3(PA,PB,PB,HA,hB,hB) = T_struct.PPPHhh;
%                     t3(PA,PB,PB,hA,hB,hB) = T_struct.PPPhhh;

                    t3(pA,PB,PB,HA,HB,HB) = T_struct.pPPHHH;

                    t3(pA,PB,PB,hA,HB,HB) = T_struct.pPPhHH;

                    t3(pA,PB,PB,HA,hB,HB) = T_struct.pPPHhH;
                    t3(pA,PB,PB,HA,HB,hB) = -permute(T_struct.pPPHhH,[1,2,3,4,6,5]);

%                     t3(pA,PB,PB,hA,hB,HB) = T_struct.pPPhhH;
%                     t3(pA,PB,PB,HA,hB,hB) = T_struct.pPPHhh;
%                     t3(pA,PB,PB,hA,hB,hB) = T_struct.pPPhhh;

                    t3(PA,PB,pB,HA,HB,HB) = T_struct.PPpHHH;        
                    t3(PA,pB,PB,HA,HB,HB) = -permute(T_struct.PPpHHH,[1,3,2,4,5,6]);

                    t3(PA,PB,pB,hA,HB,HB) = T_struct.PPphHH;
                    t3(PA,pB,PB,hA,HB,HB) = -permute(T_struct.PPphHH,[1,3,2,4,5,6]);

                    t3(PA,PB,pB,HA,hB,HB) = T_struct.PPpHhH;
                    t3(PA,PB,pB,HA,HB,hB) = -permute(T_struct.PPpHhH,[1,2,3,4,6,5]);
                    t3(PA,pB,PB,HA,hB,HB) = -permute(T_struct.PPpHhH,[1,3,2,4,5,6]);
                    t3(PA,pB,PB,HA,HB,hB) = permute(T_struct.PPpHhH,[1,3,2,4,6,5]);

%                     t3(PA,PB,pB,hA,hB,HB) = T_struct.PPphhH;
%                     t3(PA,PB,pB,HA,hB,hB) = T_struct.PPpHhh;
%                     t3(PA,PB,pB,hA,hB,hB) = T_struct.PPphhh;

%                     t3(pA,PB,pB,HA,HB,HB) = T_struct.pPpHHH;
%                     t3(pA,PB,pB,hA,HB,HB) = T_struct.pPphHH;
%                     t3(pA,PB,pB,HA,hB,HB) = T_struct.pPpHhH;
%                     t3(pA,PB,pB,hA,hB,HB) = T_struct.pPphhH;
%                     t3(pA,PB,pB,HA,hB,hB) = T_struct.pPpHhh;
%                     t3(pA,PB,pB,hA,hB,hB) = T_struct.pPphhh;
% 
%                     t3(PA,pB,pB,HA,HB,HB) = T_struct.PppHHH;
%                     t3(PA,pB,pB,hA,HB,HB) = T_struct.PpphHH;
%                     t3(PA,pB,pB,HA,hB,HB) = T_struct.PppHhH;
%                     t3(PA,pB,pB,hA,hB,HB) = T_struct.PpphhH;
%                     t3(PA,pB,pB,HA,hB,hB) = T_struct.PppHhh;
%                     t3(PA,pB,pB,hA,hB,hB) = T_struct.Ppphhh;
% 
%                     t3(pA,pB,pB,HA,HB,HB) = T_struct.pppHHH;
%                     t3(pA,pB,pB,hA,HB,HB) = T_struct.ppphHH;
%                     t3(pA,pB,pB,HA,hB,HB) = T_struct.pppHhH;
%                     t3(pA,pB,pB,hA,hB,HB) = T_struct.ppphhH;
%                     t3(pA,pB,pB,HA,hB,hB) = T_struct.pppHhh;
%                     t3(pA,pB,pB,hA,hB,hB) = T_struct.ppphhh;


            case {'D','d'}

                    t3 = zeros(Nub,Nub,Nub,Nob,Nob,Nob);
%                     t3 = t3_in;

                    PB = sys.PB; pB = sys.pB; HB = sys.HB; hB = sys.hB;

                    t3(PB,PB,PB,HB,HB,HB) = T_struct.PPPHHH;

                    t3(PB,PB,PB,hB,HB,HB) = T_struct.PPPhHH;
                    t3(PB,PB,PB,HB,hB,HB) = -permute(T_struct.PPPhHH,[1,2,3,5,4,6]);
                    t3(PB,PB,PB,HB,HB,hB) = permute(T_struct.PPPhHH,[1,2,3,5,6,4]);
%     
%                     t3(PB,PB,PB,hB,hB,HB) = T_struct.PPPhhH;
%                     t3(PB,PB,PB,hB,hB,hB) = T_struct.PPPhhh;

                    t3(PB,PB,pB,HB,HB,HB) = T_struct.PPpHHH;
                    t3(PB,pB,PB,HB,HB,HB) = -permute(T_struct.PPpHHH,[1,3,2,4,5,6]);
                    t3(pB,PB,PB,HB,HB,HB) = permute(T_struct.PPpHHH,[3,1,2,4,5,6]);

                    t3(PB,PB,pB,hB,HB,HB) = T_struct.PPphHH;
                    t3(PB,PB,pB,HB,hB,HB) = -permute(T_struct.PPphHH,[1,2,3,5,4,6]);
                    t3(PB,PB,pB,HB,HB,hB) = permute(T_struct.PPphHH,[1,2,3,5,6,4]);

                    t3(PB,pB,PB,hB,HB,HB) = -permute(T_struct.PPphHH,[1,3,2,4,5,6]);
                    t3(PB,pB,PB,HB,hB,HB) = permute(T_struct.PPphHH,[1,3,2,5,4,6]);
                    t3(PB,pB,PB,HB,HB,hB) = -permute(T_struct.PPphHH,[1,3,2,5,6,4]);

                    t3(pB,PB,PB,hB,HB,HB) = permute(T_struct.PPphHH,[3,1,2,4,5,6]);
                    t3(pB,PB,PB,HB,hB,HB) = -permute(T_struct.PPphHH,[3,1,2,5,4,6]);
                    t3(pB,PB,PB,HB,HB,hB) = permute(T_struct.PPphHH,[3,1,2,5,6,4]);

%                     t3(PB,PB,pB,hB,hB,HB) = T_struct.PPphhH;
%                     t3(PB,PB,pB,hB,hB,hB) = T_struct.PPphhh;
 
%                     t3(PB,pB,pB,HB,HB,HB) = T_struct.PppHHH;
%                     t3(PB,pB,pB,hB,HB,HB) = T_struct.PpphHH;
%                     t3(PB,pB,pB,hB,hB,HB) = T_struct.PpphhH;
%                     t3(PB,pB,pB,hB,hB,hB) = T_struct.Ppphhh;
% 
%                     t3(pB,pB,pB,HB,HB,HB) = T_struct.pppHHH;
%                     t3(pB,pB,pB,hB,HB,HB) = T_struct.ppphHH;
%                     t3(pB,pB,pB,hB,hB,HB) = T_struct.ppphhH;
%                     t3(pB,pB,pB,hB,hB,hB) = T_struct.ppphhh;


    end

end

