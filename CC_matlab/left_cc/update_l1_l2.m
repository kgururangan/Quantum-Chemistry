function [L1,L2] = update_l1_l2(X_ia,X_ijab,omega,HBar,Dia,Dijab,FM,VM,occ,unocc)

    Nocc = length(occ); Nunocc = length(unocc);
    
    L1 = zeros(Nocc,Nunocc); L2 = zeros(Nocc,Nocc,Nunocc,Nunocc);
    
    weight = 1;
    
    for i = 1:Nocc
        for a = 1:Nunocc
            
            MP = HBar{1}{2,2}(a,a) - HBar{1}{1,1}(i,i);
%               MP = FM(unocc(a),unocc(a))-FM(occ(i),occ(i));
            
            L1(i,a) = (1-weight)*L1(i,a)-weight*X_ia(i,a)./(MP - omega);
            %L1(i,a) = L1(i,a)-X_ia(i,a)./(MP - omega);


            %for j = i+1:Nocc
            for j = 1:Nocc
                %for b = a+1:Nunocc
                for b = 1:Nunocc
                    
                    MP =   HBar{1}{2,2}(a,a) - HBar{1}{1,1}(i,i)...
                          +HBar{1}{2,2}(b,b) - HBar{1}{1,1}(j,j);
%                       MP = FM(unocc(a),unocc(a))-FM(occ(i),occ(i))...
%                            +FM(unocc(b),unocc(b))-FM(occ(j),occ(j));
                
                    L2(i,j,a,b) = (1-weight)*L2(i,j,a,b)-weight*X_ijab(i,j,a,b)./(MP - omega);
                    %L2(i,j,a,b) = L2(i,j,a,b)-X_ijab(i,j,a,b)./(MP - omega);
                    
                    %L2(j,i,a,b) = -L2(i,j,a,b);
                    %L2(i,j,b,a) = -L2(i,j,a,b);
                    %L2(j,i,b,a) = L2(i,j,a,b);
                end
            end
        end
    end

end
