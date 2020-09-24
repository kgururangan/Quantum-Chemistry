function [xmin, fval, flag_conv] = golden_section(cost_fcn, xrange, epsilon, iter) 

    if nargin  == 3
        iter = 50;
    end
    
    if nargin == 2
        epsilon = 1e-6;
        iter = 50;
    end
    
    a = xrange(1); b = xrange(2);   % search interval range        
    tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
    k=0;                            % number of iterations

    x1=a+(1-tau)*(b-a);             % computing x values
    x2=a+tau*(b-a);

    f_x1=feval(cost_fcn,x1);        % computing values in x points
    f_x2=feval(cost_fcn,x2);
    
    while ((abs(b-a)>epsilon) && (k<iter))
        k=k+1;
        if(f_x1<f_x2)
            b=x2;
            x2=x1;
            x1=a+(1-tau)*(b-a);

            f_x1=feval(cost_fcn,x1);
            f_x2=feval(cost_fcn,x2);

        else
            a=x1;
            x1=x2;
            x2=a+tau*(b-a);

            f_x1=feval(cost_fcn,x1);
            f_x2=feval(cost_fcn,x2);

        end

        k=k+1;
    end

    % choose minimum point
    if(f_x1<f_x2)
        xmin = x1;
        fval = f_x1;   
    else
        xmin = x2;
        fval = f_x2;      
    end
    
    if (abs(b-a)>epsilon)
        flag_conv = true;
    end

end