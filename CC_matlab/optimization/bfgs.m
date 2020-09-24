function [xsol] = bfgs(f,grad_f,x0,H,opts)

    fprintf('\n=====================++Entering BFGS Routine++=====================\n')

    if nargin == 3
        H = eye(size(x0,1));
        tol = 1e-6;
        maxit = 100;
        alpha_range = [0, 2];
    end
    
    if nargin == 4
        tol = 1e-6;
        maxit = 100;
        alpha_range = [0, 2];
    end
    
    if nargin == 5
        tol = opts.tol;
        maxit = opts.maxit;
        alpha_range = opts.alpha;
    end


    c0 = feval(grad_f,x0);
    L = chol(H,'lower');

    it = 0;
    while norm(c0) > tol && it <= maxit

        %d = -H\c0;
        d = -L'\(L\c0);

        fun = @(alpha) feval(f,x0 + alpha*d);
        
        [alpha,~,flag_conv] = golden_section(fun,alpha_range);
        
        while ~flag_conv
            alpha_range = [alpha_range(1)*0.75, alpha_range(2)*1.25];
            printf('Golden section did not converge! Restarting with bounds [%4.4f, %4.4f]...\n',alpha_range(1),alpha_range(2))
            [alpha,~,flag_conv] = golden_section(fun,alpha_range);
        end

        x = x0 + alpha*d;

        c = feval(grad_f,x);

        dc = c - c0;
        dx = x - x0;

        D = (dc*dc')/(dc'*dx);
        E = (c0*c0')/(c0'*d);
        H = H + D + E;
        L = chol(H,'lower');

        c0 = c;
        x0 = x;

        fprintf('Iter - %d:      FVAL = %4.8f      OPTIM = %4.8f\n',it, f(x), norm(c))

        it = it + 1;

    end
    
    xsol = x;
        
    if norm(c) <= tol
        fprintf('BFGS successfuly converged in %d iterations\n',it-1)
    else
        fprintf('BFGS not converged\n')
    end
    
    
end

