function [xout] = asym(xin,spincase)


    rank = ndims(xin)/2;
    
    switch rank
        
        case 1
            
            xout = xin;
        
        case 2
            
            switch spincase 
                
                case 'a'
                    
                    t2a = zeros(size(xin));
                    
                    for a = 1:size(xin,1)
                        for b = a+1:size(xin,2)
                            for i = 1:size(xin,3)
                                for j = i+1:size(xin,4)
                                    t2a(a,b,i,j) = xin(a,b,i,j);
                                    t2a(b,a,i,j) = -t2a(a,b,i,j);
                                    t2a(a,b,j,i) = -t2a(a,b,i,j);
                                    t2a(b,a,j,i) = t2a(a,b,i,j);
                                end
                            end
                        end
                    end
                    
                    xout = t2a;
                    
                case 'b'
                    
                    xout = xin;
                    
                case 'c'
                    
                    t2c = zeros(size(xin));
                    
                    for a = 1:size(xin,1)
                        for b = a+1:size(xin,2)
                            for i = 1:size(xin,3)
                                for j = i+1:size(xin,4)
                                    t2c(a,b,i,j) = xin(a,b,i,j);
                                    t2c(b,a,i,j) = -t2c(a,b,i,j);
                                    t2c(a,b,j,i) = -t2c(a,b,i,j);
                                    t2c(b,a,j,i) = t2c(a,b,i,j);
                                end
                            end
                        end
                    end
                    
                    xout = t2c;
                    
            end
            
        case 3
            
            switch spincase
                
                case 'a'
                    
                        t3a = zeros(size(xin));

                        for a = 1:size(xin,1)
                            for b = a+1:size(xin,2)
                                for c = b+1:size(xin,3)
                                    for i = 1:size(xin,4)
                                        for j = i+1:size(xin,5)
                                            for k = j+1:size(xin,6)

                                                % (1)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]  
                                                t3a(a,b,c,i,j,k) = xin(a,b,c,i,j,k);
                                                t3a(a,b,c,k,i,j) = t3a(a,b,c,i,j,k);
                                                t3a(a,b,c,j,k,i) = t3a(a,b,c,i,j,k);
                                                t3a(a,b,c,i,k,j) = -t3a(a,b,c,i,j,k);
                                                t3a(a,b,c,j,i,k) = -t3a(a,b,c,i,j,k);
                                                t3a(a,b,c,k,j,i) = -t3a(a,b,c,i,j,k);

                                                % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3a(b,a,c,i,j,k) = -t3a(a,b,c,i,j,k);
                                                t3a(b,a,c,k,i,j) = -t3a(a,b,c,i,j,k);
                                                t3a(b,a,c,j,k,i) = -t3a(a,b,c,i,j,k);
                                                t3a(b,a,c,i,k,j) = t3a(a,b,c,i,j,k);
                                                t3a(b,a,c,j,i,k) = t3a(a,b,c,i,j,k);
                                                t3a(b,a,c,k,j,i) = t3a(a,b,c,i,j,k);

                                                % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3a(a,c,b,i,j,k) = -t3a(a,b,c,i,j,k);
                                                t3a(a,c,b,k,i,j) = -t3a(a,b,c,i,j,k);
                                                t3a(a,c,b,j,k,i) = -t3a(a,b,c,i,j,k);
                                                t3a(a,c,b,i,k,j) = t3a(a,b,c,i,j,k);
                                                t3a(a,c,b,j,i,k) = t3a(a,b,c,i,j,k);
                                                t3a(a,c,b,k,j,i) = t3a(a,b,c,i,j,k);

                                                % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3a(c,b,a,i,j,k) = -t3a(a,b,c,i,j,k);
                                                t3a(c,b,a,k,i,j) = -t3a(a,b,c,i,j,k);
                                                t3a(c,b,a,j,k,i) = -t3a(a,b,c,i,j,k);
                                                t3a(c,b,a,i,k,j) = t3a(a,b,c,i,j,k);
                                                t3a(c,b,a,j,i,k) = t3a(a,b,c,i,j,k);
                                                t3a(c,b,a,k,j,i) = t3a(a,b,c,i,j,k);

                                                % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3a(b,c,a,i,j,k) = t3a(a,b,c,i,j,k);
                                                t3a(b,c,a,k,i,j) = t3a(a,b,c,i,j,k);
                                                t3a(b,c,a,j,k,i) = t3a(a,b,c,i,j,k);
                                                t3a(b,c,a,i,k,j) = -t3a(a,b,c,i,j,k);
                                                t3a(b,c,a,j,i,k) = -t3a(a,b,c,i,j,k);
                                                t3a(b,c,a,k,j,i) = -t3a(a,b,c,i,j,k);

                                                % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3a(c,a,b,i,j,k) = t3a(a,b,c,i,j,k);
                                                t3a(c,a,b,k,i,j) = t3a(a,b,c,i,j,k);
                                                t3a(c,a,b,j,k,i) = t3a(a,b,c,i,j,k);
                                                t3a(c,a,b,i,k,j) = -t3a(a,b,c,i,j,k);
                                                t3a(c,a,b,j,i,k) = -t3a(a,b,c,i,j,k);
                                                t3a(c,a,b,k,j,i) = -t3a(a,b,c,i,j,k);
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        xout = t3a;
                    
                case 'b'
                    
                        t3b = zeros(size(xin));
                        
                        for a = 1:size(xin,1)
                           for b = a+1:size(xin,2)
                               for c = 1:size(xin,3)
                                   for i = 1:size(xin,4)
                                       for j = i+1:size(xin,5)
                                           for k = 1:size(xin,6)
                                                t3b(a,b,c,i,j,k) = xin(a,b,c,i,j,k);
                                                t3b(b,a,c,i,j,k) = -t3b(a,b,c,i,j,k);
                                                t3b(a,b,c,j,i,k) = -t3b(a,b,c,i,j,k);
                                                t3b(b,a,c,j,i,k) = t3b(a,b,c,i,j,k);
                                           end
                                       end
                                   end
                               end
                            end
                        end
                       
                        xout = t3b;
            
                case 'c'
                    
                         t3c = zeros(size(xin));                        
                         for a = 1:size(xin,1)
                            for b = 1:size(xin,2)
                                for c = b+1:size(xin,3)
                                    for i = 1:size(xin,4)
                                        for j = 1:size(xin,5)
                                            for k = j+1:size(xin,6)
                                                t3c(a,b,c,i,j,k) = xin(a,b,c,i,j,k);
                                                t3c(a,c,b,i,j,k) = -t3c(a,b,c,i,j,k);
                                                t3c(a,b,c,i,k,j) = -t3c(a,b,c,i,j,k);
                                                t3c(a,c,b,i,k,j) = t3c(a,b,c,i,j,k);
                                            end
                                        end
                                    end
                                end
                            end
                         end
                        
                         xout = t3c;
                         
                case 'd'
                    
                         t3d = zeros(size(xin));                        
                         for a = 1:size(xin,1)
                            for b = a+1:size(xin,2)
                                for c = b+1:size(xin,3)
                                    for i = 1:size(xin,4)
                                        for j = i+1:size(xin,5)
                                            for k = j+1:size(xin,6)

                                                % (1)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]  
                                                t3d(a,b,c,i,j,k) = xin(a,b,c,i,j,k);
                                                t3d(a,b,c,k,i,j) = t3d(a,b,c,i,j,k);
                                                t3d(a,b,c,j,k,i) = t3d(a,b,c,i,j,k);
                                                t3d(a,b,c,i,k,j) = -t3d(a,b,c,i,j,k);
                                                t3d(a,b,c,j,i,k) = -t3d(a,b,c,i,j,k);
                                                t3d(a,b,c,k,j,i) = -t3d(a,b,c,i,j,k);

                                                % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3d(b,a,c,i,j,k) = -t3d(a,b,c,i,j,k);
                                                t3d(b,a,c,k,i,j) = -t3d(a,b,c,i,j,k);
                                                t3d(b,a,c,j,k,i) = -t3d(a,b,c,i,j,k);
                                                t3d(b,a,c,i,k,j) = t3d(a,b,c,i,j,k);
                                                t3d(b,a,c,j,i,k) = t3d(a,b,c,i,j,k);
                                                t3d(b,a,c,k,j,i) = t3d(a,b,c,i,j,k);

                                                % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3d(a,c,b,i,j,k) = -t3d(a,b,c,i,j,k);
                                                t3d(a,c,b,k,i,j) = -t3d(a,b,c,i,j,k);
                                                t3d(a,c,b,j,k,i) = -t3d(a,b,c,i,j,k);
                                                t3d(a,c,b,i,k,j) = t3d(a,b,c,i,j,k);
                                                t3d(a,c,b,j,i,k) = t3d(a,b,c,i,j,k);
                                                t3d(a,c,b,k,j,i) = t3d(a,b,c,i,j,k);

                                                % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3d(c,b,a,i,j,k) = -t3d(a,b,c,i,j,k);
                                                t3d(c,b,a,k,i,j) = -t3d(a,b,c,i,j,k);
                                                t3d(c,b,a,j,k,i) = -t3d(a,b,c,i,j,k);
                                                t3d(c,b,a,i,k,j) = t3d(a,b,c,i,j,k);
                                                t3d(c,b,a,j,i,k) = t3d(a,b,c,i,j,k);
                                                t3d(c,b,a,k,j,i) = t3d(a,b,c,i,j,k);

                                                % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3d(b,c,a,i,j,k) = t3d(a,b,c,i,j,k);
                                                t3d(b,c,a,k,i,j) = t3d(a,b,c,i,j,k);
                                                t3d(b,c,a,j,k,i) = t3d(a,b,c,i,j,k);
                                                t3d(b,c,a,i,k,j) = -t3d(a,b,c,i,j,k);
                                                t3d(b,c,a,j,i,k) = -t3d(a,b,c,i,j,k);
                                                t3d(b,c,a,k,j,i) = -t3d(a,b,c,i,j,k);

                                                % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                                t3d(c,a,b,i,j,k) = t3d(a,b,c,i,j,k);
                                                t3d(c,a,b,k,i,j) = t3d(a,b,c,i,j,k);
                                                t3d(c,a,b,j,k,i) = t3d(a,b,c,i,j,k);
                                                t3d(c,a,b,i,k,j) = -t3d(a,b,c,i,j,k);
                                                t3d(c,a,b,j,i,k) = -t3d(a,b,c,i,j,k);
                                                t3d(c,a,b,k,j,i) = -t3d(a,b,c,i,j,k);
                                            end
                                        end
                                    end
                                end
                            end
                         end
                        
                         xout = t3d;
            end
            
        otherwise
            fprintf('Tensor of rank %d not supported\n',rank)
            
    end
                    


    

  