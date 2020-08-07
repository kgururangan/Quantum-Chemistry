function [l1a_new,l1b_new,l2a_new,l2b_new,l2c_new] = update_L(X1A,X1B,X2A,X2B,X2C,HBar_t,sys,omega,shift)

    Nocc_a = sys.Nocc_alpha; Nocc_b = sys.Nocc_beta; 
    Nunocc_a = sys.Nvir_alpha; Nunocc_b = sys.Nvir_beta;
    
    H1A = HBar_t.H1A; H1B = HBar_t.H1B;

    l1a_new = zeros(Nunocc_a, Nocc_a);
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            denom = H1A.vv(a,a) - H1A.oo(i,i);
            l1a_new(a,i) = l1a_new(a,i)-X1A(a,i)/(denom - omega + shift); 
        end
    end
    
    l1b_new = zeros(Nunocc_b,Nocc_b);
    for i = 1:Nocc_b
        for a = 1:Nunocc_b
            denom = H1B.vv(a,a) - H1B.oo(i,i);
            l1b_new(a,i) = l1b_new(a,i)-X1B(a,i)/(denom - omega + shift); 
        end
    end
    
    l2a_new = zeros(Nunocc_a, Nunocc_a, Nocc_a, Nocc_a);
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            for j = i+1:Nocc_a
                for b = a+1:Nunocc_a
                    
                    denom = H1A.vv(a,a) - H1A.oo(i,i) + H1A.vv(b,b) - H1A.oo(j,j);
                      
                    l2a_new(a,b,i,j) = l2a_new(a,b,i,j)-X2A(a,b,i,j)/(denom - omega + shift);                   
                    l2a_new(b,a,i,j) = -l2a_new(a,b,i,j);
                    l2a_new(a,b,j,i) = -l2a_new(a,b,i,j);
                    l2a_new(b,a,j,i) = l2a_new(a,b,i,j);
                end
            end
        end
    end
    
   l2b_new = zeros(Nunocc_a, Nunocc_b, Nocc_a, Nocc_b);
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            for j = 1:Nocc_b
                for b = 1:Nunocc_b                    
                    denom = H1A.vv(a,a) - H1A.oo(i,i) + H1B.vv(b,b) - H1B.oo(j,j);                      
                    l2b_new(a,b,i,j) = l2b_new(a,b,i,j)-X2B(a,b,i,j)/(denom - omega + shift);                   
                end
            end
        end
    end
    
    l2c_new = zeros(Nunocc_b, Nunocc_b, Nocc_b, Nocc_b);
    for i = 1:Nocc_b
        for a = 1:Nunocc_b
            for j = i+1:Nocc_b
                for b = a+1:Nunocc_b
                    
                    denom = H1B.vv(a,a) - H1B.oo(i,i) + H1B.vv(b,b) - H1B.oo(j,j);
                      
                    l2c_new(a,b,i,j) = l2c_new(a,b,i,j)-X2C(a,b,i,j)/(denom - omega + shift);                   
                    l2c_new(b,a,i,j) = -l2c_new(a,b,i,j);
                    l2c_new(a,b,j,i) = -l2c_new(a,b,i,j);
                    l2c_new(b,a,j,i) = l2c_new(a,b,i,j);
                end
            end
        end
    end


end

