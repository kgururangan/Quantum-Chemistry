function [] = print_R(R,sys,thresh)

    if nargin < 3
        thresh = 0.1;
    end
    
    Nocc_a = sys.Nocc_alpha;
    Nocc_b = sys.Nocc_beta;

    r1a = reshape(R(sys.posv{1}),sys.size{1});
    r1b = reshape(R(sys.posv{2}),sys.size{2});
    r2a = reshape(R(sys.posv{3}),sys.size{3});
    r2b = reshape(R(sys.posv{4}),sys.size{4});
    r2c = reshape(R(sys.posv{5}),sys.size{5});
    
    for a = 1:size(r1a,1)
        for i = 1:size(r1a,2)
            if abs(r1a(a,i)) > thresh
                fprintf('%dA -> %dA      %4.8f\n',i,a+Nocc_a,r1a(a,i))
            end
        end
    end
    
   for a = 1:size(r1b,1)
        for i = 1:size(r1b,2)
            if abs(r1b(a,i)) > thresh
                fprintf('%dB -> %dB      %4.8f\n',i,a+Nocc_b,r1b(a,i))
            end
        end
   end
    
   for a = 1:size(r2a,1)
        for b = 1:size(r2a,2)
            for i = 1:size(r2a,3)
                for j = 1:size(r2a,4)
                    if abs(r2a(a,b,i,j)) > thresh
                        fprintf('%dA  %dA -> %dA  %dA    %4.8f\n',i,j,a+Nocc_a,b+Nocc_a,r2a(a,b,i,j))
                    end
                end
            end
        end
   end
    
   for a = 1:size(r2b,1)
        for b = 1:size(r2b,2)
            for i = 1:size(r2b,3)
                for j = 1:size(r2b,4)
                    if abs(r2b(a,b,i,j)) > thresh
                        fprintf('%dA  %dB -> %dA  %dB    %4.8f\n',i,j,a+Nocc_a,b+Nocc_b,r2b(a,b,i,j))
                    end
                end
            end
        end
   end
    
   for a = 1:size(r2c,1)
        for b = 1:size(r2c,2)
            for i = 1:size(r2c,3)
                for j = 1:size(r2c,4)
                    if abs(r2c(a,b,i,j)) > thresh
                        fprintf('%dA  %dB -> %dA  %dB    %4.8f\n',i,j,a+Nocc_b,b+Nocc_b,r2c(a,b,i,j))
                    end
                end
            end
        end
    end

end

