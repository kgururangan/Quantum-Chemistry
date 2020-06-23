function [t1, t2, t3] = update_t1_t2_t3(X_ai,X_abij,X_abcijk,VM,FM,occ,unocc)
    
    Nocc = length(occ);
    Nunocc = length(unocc);
    
    foo = FM(occ,occ);
    fvv = FM(unocc,unocc);
    
    t1 = zeros(Nunocc,Nocc);
    t2 = zeros(Nunocc,Nunocc,Nocc,Nocc);
    t3 = zeros(Nunocc,Nunocc,Nunocc,Nocc,Nocc,Nocc);
    
    for a = 1:Nunocc
        for i = 1:Nocc
            
            t1(a,i) = X_ai(a,i)/(foo(i,i)-fvv(a,a));
            
            for b = a+1:Nunocc
                for j = i+1:Nocc
                    t2(a,b,i,j) = X_abij(a,b,i,j)/(foo(i,i)+foo(j,j)-fvv(a,a)-fvv(b,b));
                    t2(b,a,i,j) = -t2(a,b,i,j);
                    t2(a,b,j,i) = -t2(a,b,i,j);
                    t2(b,a,j,i) = t2(a,b,i,j);
                    
                    for c = b+1:Nunocc
                        for k = j+1:Nocc
                            
                            % (1)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(a,b,c,i,j,k) = X_abcijk(a,b,c,i,j,k)/...
                                (foo(i,i)+foo(j,j)+foo(k,k)-fvv(a,a)-fvv(b,b)-fvv(c,c));
                            t3(a,b,c,k,i,j) = t3(a,b,c,i,j,k);
                            t3(a,b,c,j,k,i) = t3(a,b,c,i,j,k);
                            t3(a,b,c,i,k,j) = -t3(a,b,c,i,j,k);
                            t3(a,b,c,j,i,k) = -t3(a,b,c,i,j,k);
                            t3(a,b,c,k,j,i) = -t3(a,b,c,i,j,k);
                            % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(b,a,c,i,j,k) = -t3(a,b,c,i,j,k);
                            t3(b,a,c,k,i,j) = -t3(a,b,c,i,j,k);
                            t3(b,a,c,j,k,i) = -t3(a,b,c,i,j,k);
                            t3(b,a,c,i,k,j) = t3(a,b,c,i,j,k);
                            t3(b,a,c,j,i,k) = t3(a,b,c,i,j,k);
                            t3(b,a,c,k,j,i) = t3(a,b,c,i,j,k);
                            % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(a,c,b,i,j,k) = -t3(a,b,c,i,j,k);
                            t3(a,c,b,k,i,j) = -t3(a,b,c,i,j,k);
                            t3(a,c,b,j,k,i) = -t3(a,b,c,i,j,k);
                            t3(a,c,b,i,k,j) = t3(a,b,c,i,j,k);
                            t3(a,c,b,j,i,k) = t3(a,b,c,i,j,k);
                            t3(a,c,b,k,j,i) = t3(a,b,c,i,j,k);
                            % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(c,b,a,i,j,k) = -t3(a,b,c,i,j,k);
                            t3(c,b,a,k,i,j) = -t3(a,b,c,i,j,k);
                            t3(c,b,a,j,k,i) = -t3(a,b,c,i,j,k);
                            t3(c,b,a,i,k,j) = t3(a,b,c,i,j,k);
                            t3(c,b,a,j,i,k) = t3(a,b,c,i,j,k);
                            t3(c,b,a,k,j,i) = t3(a,b,c,i,j,k);
                            % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(b,c,a,i,j,k) = t3(a,b,c,i,j,k);
                            t3(b,c,a,k,i,j) = t3(a,b,c,i,j,k);
                            t3(b,c,a,j,k,i) = t3(a,b,c,i,j,k);
                            t3(b,c,a,i,k,j) = -t3(a,b,c,i,j,k);
                            t3(b,c,a,j,i,k) = -t3(a,b,c,i,j,k);
                            t3(b,c,a,k,j,i) = -t3(a,b,c,i,j,k);
                            % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                            t3(c,a,b,i,j,k) = t3(a,b,c,i,j,k);
                            t3(c,a,b,k,i,j) = t3(a,b,c,i,j,k);
                            t3(c,a,b,j,k,i) = t3(a,b,c,i,j,k);
                            t3(c,a,b,i,k,j) = -t3(a,b,c,i,j,k);
                            t3(c,a,b,j,i,k) = -t3(a,b,c,i,j,k);
                            t3(c,a,b,k,j,i) = -t3(a,b,c,i,j,k);
                        end
                    end
                end
            end
        end
    end 
end
