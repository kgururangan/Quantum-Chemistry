function [t1, t2] = update_t1_t2(X_ai,X_abij,VM,FM,occ,unocc)
    
    Nocc = length(occ);
    Nunocc = length(unocc);
    
    foo = FM(occ,occ);
    fvv = FM(unocc,unocc);
    
    t1 = zeros(Nunocc,Nocc);
    t2 = zeros(Nunocc,Nunocc,Nocc,Nocc);
    
    for a = 1:Nunocc
        for i = 1:Nocc
            
            t1(a,i) = X_ai(a,i)/(foo(i,i)-fvv(a,a));
            
            for b = a+1:Nunocc
                for j = i+1:Nocc
                    t2(a,b,i,j) = X_abij(a,b,i,j)/(foo(i,i)+foo(j,j)-fvv(a,a)-fvv(b,b));
                    t2(b,a,i,j) = -t2(a,b,i,j);
                    t2(a,b,j,i) = -t2(a,b,i,j);
                    t2(b,a,j,i) = t2(a,b,i,j);
                end
            end
            
            
        end
    end
    
end
