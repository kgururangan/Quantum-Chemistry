function [xout] = asym_act_t3(T3,spincase)

    switch spincase
        
        case 'a'
            
            % 1 - PPPHHH
            xout.PPPHHH = asym(T3.PPPHHH,spincase);
            
            % 2 - PPpHHH
            x0 = T3.PPpHHH;
            x1 = zeros(size(x0));
            for a = 1:size(x0,1)
                for b = a+1:size(x0,2)
                    for c = 1:size(x0,3)
                        for i = 1:size(x0,4)
                            for j = i+1:size(x0,5)
                                for k = j+1:size(x0,6)
                                    x1(a,b,c,i,j,k) = x0(a,b,c,i,j,k);
                                    x1(a,b,c,j,i,k) = -x0(a,b,c,i,j,k);
                                    x1(a,b,c,i,k,j) = -x0(a,b,c,i,j,k);
                                    x1(a,b,c,k,j,i) = x0(a,b,c,i,j,k);
                                    x1(a,b,c,k,i,j) = -x0(a,b,c,i,j,k);
                                    x1(a,b,c,j,k,i) = x0(a,b,c,i,j,k);
                                    
                                    x1(b,a,c,i,j,k) = -x0(a,b,c,i,j,k);
                                    x1(b,a,c,j,i,k) = x0(a,b,c,i,j,k);
                                    x1(b,a,c,i,k,j) = x0(a,b,c,i,j,k);
                                    x1(b,a,c,k,j,i) = -x0(a,b,c,i,j,k);
                                    x1(b,a,c,k,i,j) = x0(a,b,c,i,j,k);
                                    x1(b,a,c,j,k,i) = -x0(a,b,c,i,j,k);
                                end
                            end
                        end
                    end
                end
            end
            xout.PPpHHH = x1;
            
            % 3 - PPPhHH
            x0 = T3.PPPhHH;
            x1 = zeros(size(x0));
            for a = 1:size(x0,1)
                for b = a+1:size(x0,2)
                    for c = b+1:size(x0,3)
                        for i = 1:size(x0,4)
                            for j = 1:size(x0,5)
                                for k = j+1:size(x0,6)
                                    x1(a,b,c,i,j,k) = x0(a,b,c,i,j,k);
                                    x1(b,a,c,i,j,k) = -x0(a,b,c,i,j,k);
                                    x1(a,c,b,i,j,k) = -x0(a,b,c,i,j,k);
                                    x1(c,b,a,i,j,k) = x0(a,b,c,i,j,k);
                                    x1(c,a,b,i,j,k) = -x0(a,b,c,i,j,k);
                                    x1(c,b,a,i,j,k) = x0(a,b,c,i,j,k);
                                    
                                    x1(a,b,c,i,k,j) = -x0(a,b,c,i,j,k);
                                    x1(b,a,c,i,k,j) = x0(a,b,c,i,j,k);
                                    x1(a,c,b,i,k,j) = x0(a,b,c,i,j,k);
                                    x1(c,b,a,i,k,j) = -x0(a,b,c,i,j,k);
                                    x1(c,a,b,i,k,j) = x0(a,b,c,i,j,k);
                                    x1(c,b,a,i,k,j) = -x0(a,b,c,i,j,k);               
                                end
                            end
                        end
                    end
                end
            end
            xout.PPPhHH = x1;
            
            % 4 - PPphHH
            x0 = T3.PPphHH;
            x1 = zeros(size(x0));
            for a = 1:size(x0,1)
                for b = a+1:size(x0,2)
                    for c = 1:size(x0,3)
                        for i = 1:size(x0,4)
                            for j = 1:size(x0,5)
                                for k = j+1:size(x0,6)
                                    x1(a,b,c,i,j,k) = x0(a,b,c,i,j,k);
                                    x1(b,a,c,i,j,k) = -x0(a,b,c,i,j,k);
                                    x1(a,b,c,i,k,j) = -x0(a,b,c,i,j,k);
                                    x1(b,a,c,i,k,j) = x0(a,b,c,i,j,k);
                                end
                            end
                        end
                    end
                end
            end
            xout.PPphHH = x1;
            
            
        case 'b'
            
            % 1 - PPPHHH
            xout.PPPHHH = asym(T3.PPPHHH,spincase);
            
            % 2 - PPpHHH
            x0 = T3.PPpHHH;
            x1 = zeros(size(x0));
            for a = 1:size(x0,1)
                for b = a+1:size(x0,2)
                    for c = 1:size(x0,3)
                        for i = 1:size(x0,4)
                            for j = i+1:size(x0,5)
                                for k = 1:size(x0,6)
                                    x1(a,b,c,i,j,k) = x0(a,b,c,i,j,k);
                                    x1(b,a,c,i,j,k) = -x0(a,b,c,i,j,k);
                                    x1(a,b,c,j,i,k) = -x0(a,b,c,i,j,k);
                                    x1(b,a,c,j,i,k) = x0(a,b,c,i,j,k);
                                end
                            end
                        end
                    end
                end
            end
            xout.PPpHHH = x1;
            
            % 3 - PpPHHH
            x0 = T3.PpPHHH;
            x1 = zeros(size(x0));
            for a = 1:size(x0,1)
                for b = 1:size(x0,2)
                    for c = 1:size(x0,3)
                        for i = 1:size(x0,4)
                            for j = i+1:size(x0,5)
                                for k = j+1:size(x0,6)
                                    x1(a,b,c,i,j,k) = x0(a,b,c,i,j,k);
                                    x1(a,b,c,j,i,k) = -x0(a,b,c,i,j,k);
                                end
                            end
                        end
                    end
                end
            end
            xout.PpPHHH = x1;
            
            % 4 - PPPHHh
            x0 = T3.PPPHHh;
            x1 = zeros(size(x0));
            for a = 1:size(x0,1)
                for b = a+1:size(x0,2)
                    for c = 1:size(x0,3)
                        for i = 1:size(x0,4)
                            for j = i+1:size(x0,5)
                                for k = 1:size(x0,6)
                                    x1(a,b,c,i,j,k) = x0(a,b,c,i,j,k);
                                    x1(b,a,c,i,j,k) = -x0(a,b,c,i,j,k);
                                    x1(a,b,c,j,i,k) = -x0(a,b,c,i,j,k);
                                    x1(b,a,c,j,i,k) = x0(a,b,c,i,j,k);
                                end
                            end
                        end
                    end
                end
            end
            xout.PPPHHh = x1;
            
            % 5 - PPPhHH
            
            % 6 - PPphHH
            
            % 7 - PPpHHh
            
            % 8 - PpPhHH
            
            % 9 - PpPHHh
            
            
        case 'c'
            
            % 1 - PPPHHH
            
            % 2 - PPpHHH
            
            % 3 - pPPHHH
            
            % 4 - PPPHHh
            
            % 5 - PPPhHH
            
            % 6 - PPphHH
            
            % 7 - PPpHHh
            
            % 8 - PpPhHH
            
            % 9 - PpPHHh
            
            
        case 'd'
            
            % 1 - PPPHHH
            
            % 2 - PPpHHH
            
            % 3 - PPPhHH
            
            % 4 - PPphHH
            
            
        otherwise
            
            disp('Enter valid spin case!')
            
    end