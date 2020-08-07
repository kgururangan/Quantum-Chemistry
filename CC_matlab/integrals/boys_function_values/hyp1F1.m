function [hg] = hyp1F1(a,b,x)

    % Confluent hypergoemetric function 1F1(a,b,x) numerical implementation
    % copied from FORTRAN77 code found somewhere... I forget exactly where
    
    

	
       a0 = a;
       a1 = a;
       x0 = x;
       hg = 0.0;

       if (b == 0.0 || b == -abs(b)) 
               hg = 1.0e+300;
       elseif (a==0.0 || x==0.0) 
               hg = 1.0;
       elseif (a==-1.0) 
               hg = 1.0 - x/b;
       elseif (a==b) 
               hg = exp(x);
       elseif (a-b==1.0) 
               hg = (1.0 + x/b)*exp(x);
       elseif (a==1.0 && b==2.0) 
               hg = (exp(x)-1.0)/x;
       elseif (a==int32(a) && a<0.0) 
               m = int32(-a);
               r = 1.0;
               hg = 1.0;
               for k = 1:m
                  r = r*(a+k-1.0)/k/(b+k-1.0)*x;
                  hg = hg + r;
               end
       end
       
%        if any(isinf([a0, a1, x0, hg, b, x]))
%            fprintf('NaN detected!')
%        end

       if (hg ~= 0.0) 
              return
       end

       if (x<0.0) 
                a = b - a;
                a0 = a;
                x = abs(x);
       end

       if (a < 2.0) 
              n1 = 0;
       end

       if (2.0 <= a) 
              n1 = 1;
              la = int32(a);
              a = a - la - 1.0;
       end
       
       

       for n = 0:n1
           %x
           
           if (2.0 <= a0) 
                  a = a + 1.0;
           end

          if (x <= 30.0 + abs(b) || a < 0.0 ) 
                  hg = 1.0;
                  rg = 1.0;
                  for j = 1:500
                        rg = rg*(a+j-1.0)...
                             /(j*(b+j-1.0))*x;
                        hg = hg + rg;
                        if (abs(rg/hg) < eps) 
                                break
                        end
                  end
          else
              
                   
                  
                  ta = gamma(a);
                  tb = gamma(b);
                  xg = b-a;
                  tba = gamma(xg);
                  sum1 = 1.0;
                  sum2 = 1.0;
                  r1 = 1.0;
                  r2 = 1.0;
                  for i = 1:8
                        r1 =-r1*(a+i-1.0)*(a-b+i)/(x*i);
                        r2 =-r2*(b-a+i-1.0)*(a-i)/(x*i);
                        sum1 = sum1+r1;
                        sum2 = sum2+r2;
                  end
                  hg1 = tb/tba*power(x,-a)*cos(pi*a)*sum1;
                  
%                   EXPX = exp(x);
%                   if isinf(EXPX)
%                       EXPX = realmax('double') * 1e-1;
%                   end
                  
                  hg2 = tb/ta*exp(x)*power(x,(a-b))*sum2;
                  
%                   if isinf(hg2)
% %                       hg2 = realmax('double');
%                       hg2 = 0.0; 
% %                       hg = 0.0;
%                       return
% %                       x
% %                       loghg2 = log(tb)-log(ta)+x+(a-b)*log(x)+log(sum2);
% %                       hg2 = exp(loghg2)
%                   end
                  
                  hg=hg1+hg2;
                  
                  %%%% THERE CAN BE A NAN PROBLEM HERE BECAUSE hg2 TAKES ON
                  %%%% THE VALUE OF INF. THIS IS BECAUSE x CAN BE A LARGER
                  %%%% NUMBER (E.G. 7.842068532021935e+02) AND exp(x)
                  %%%% RETURNS INF AS A RESULT
                  
%                   if isinf(hg2)    
%                       disp('NaN detected!')
%                   end

          end
       end

         if (2.0 <= a0) 
                 for i=1:la-1
                        hg=((2.0*a-b+x)*y1+(b-a)*y0)/a;
                        y0=y1;
                        y1=hg;
                        a=a+1.0;
                 end
                 
         end
         
%          if any(isnan([x0,hg]))
%              disp('NaN detected!')
%          end
             
%         hg

         if (x0<0.0) 
                 hg=hg*exp(x0);
         end

         a=a1;
         x=x0;

         return

end 
