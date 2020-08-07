function [D, Dia, Dijab] = HBar_CCSD_SD_diagonal(HBar, t1, t2, sys)

%     occ = sys.occ;
%     unocc = sys.unocc;
%     
%     VM = sys.VM;
%     
%     Nocc = length(occ);
%     Nunocc = length(unocc);
    
    Nocc = sys.Nocc;
    Nunocc = sys.Nunocc;

    Dia = zeros(Nunocc,Nocc);
    Dijab = zeros(Nunocc,Nunocc,Nocc,Nocc);
    
    for a = 1:Nunocc
        for i = 1:Nocc
            
            Dia(a,i) = Dia(a,i) + HBar{1}{2,2}(a,a) - HBar{1}{1,1}(i,i);
            
            % WHY IS THIS A MINUS???
            Dia(a,i) = Dia(a,i) - HBar{2}{1,2,2,1}(i,a,a,i);
            
            for b = 1:Nunocc
                for j = 1:Nocc
                    
                    % these terms come in paris of A(ab) and A(ij) but the sign is +1 for each antisymmetrizer
                    % because determinants on both the left and right flip sign
                    Dijab(a,b,i,j) = Dijab(a,b,i,j) + ...
                                     HBar{1}{2,2}(a,a) + HBar{1}{2,2}(b,b) - ...
                                     HBar{1}{1,1}(i,i) - HBar{1}{1,1}(j,j);


                                 
                   %  note these 4 terms come from A(ij)A(ab) but the signs are all +1 because both determinants in 
                   %  <ijab|Hbar|ijab> get flipped so the antisymmetrizer sign is always +1
                   %  WHY IS THIS A MINUS???
                   %  Must be to do with the inherent sign of W_mbej... there is no funny business here!
                   
                    Dijab(a,b,i,j) = Dijab(a,b,i,j) + ...
                                    HBar{2}{1,2,2,1}(i,b,b,i) - ...
                                    HBar{2}{1,2,2,1}(j,b,b,j) - ...
                                    HBar{2}{1,2,2,1}(i,a,a,i) - ...
                                    HBar{2}{1,2,2,1}(j,a,a,j); 
                                
                                
                    Dijab(a,b,i,j) = Dijab(a,b,i,j) + HBar{2}{2,2,2,2}(a,b,a,b) + HBar{2}{1,1,1,1}(i,j,i,j);
                    
%                     for l = 1:Nocc
%                         Dijab(a,b,i,j) = Dijab(a,b,i,j) - ...
%                                          0.5*(VM(occ(j),occ(l),unocc(a),unocc(b))*t2(a,b,j,l) + ...
%                                               VM(occ(i),occ(l),unocc(a),unocc(b))*t2(a,b,i,l));
%                     end
%                     
%                     for d = 1:Nunocc
%                         Dijab(a,b,i,j) = Dijab(a,b,i,j) - ...
%                                          0.5*(VM(occ(i),occ(j),unocc(b),unocc(d))*t2(b,d,i,j) + ...
%                                               VM(occ(i),occ(j),unocc(a),unocc(d))*t2(a,d,i,j));
%                     end
                    
                end
                
            end
            
        end
        
    end
    
    D = cat(1,Dia(:),Dijab(:));

end
% 
%     for a in range(Nunocc):
%         for i in range(Nocc):
%             
%             Dia[a,i] += F_ae[a,a] - F_mi[i,i]
%             
%             # WHY IS THIS A MINUS???
%             Dia[a,i] -= W_mbej[i,a,a,i]
%             
%             for b in range(Nunocc):
%                 for j in range(Nocc):
%                     
%                     # these terms come in paris of A(ab) and A(ij) but the sign is +1 for each antisymmetrizer
%                     # because determinants on both the left and right flip sign
%                     Dijab[a,b,i,j] += F_ae[a,a] + F_ae[b,b] - F_mi[i,i] - F_mi[j,j]
%                     
%                     # note these 4 terms come from A(ij)A(ab) but the signs are all +1 because both determinants in 
%                     # <ijab|Hbar|ijab> get flipped so the antisymmetrizer sign is always +1
%                     # WHY IS THIS A MINUS???
%                     # Must be to do with the inherent sign of W_mbej... there is no funny business here!
%                     Dijab[a,b,i,j] -= W_mbej[i,b,b,i] + W_mbej[j,b,b,j] + W_mbej[i,a,a,i] + W_mbej[j,a,a,j] 
%                     
%                     Dijab[a,b,i,j] += W_abef[a,b,a,b] + W_mnij[i,j,i,j]
%                     
%                     for l in range(Nocc):
%                         Dijab[a,b,i,j] -= 0.5*(VM[occ[j],occ[l],unocc[a],unocc[b]]*t2[a,b,j,l] + \
%                                                VM[occ[i],occ[l],unocc[a],unocc[b]]*t2[a,b,i,l])
%                     for d in range(Nunocc):
%                         Dijab[a,b,i,j] -= 0.5*(VM[occ[i],occ[j],unocc[b],unocc[d]]*t2[b,d,i,j] + \
%                                                VM[occ[i],occ[j],unocc[a],unocc[d]]*t2[a,d,i,j])
%     
%                     
%     D = np.hstack(( Dia.flatten(), Dijab.flatten()))
%     idx = np.argsort(D)
%     return D
