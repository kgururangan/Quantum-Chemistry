function [T] = AiryTMM(Potential,Range,Energy,Mass,V_Left,V_Right,NumLines)

%     svd_thresh = 1e-10;
%     cond_thresh = 1e4;
    
    [x_interp, V_interp, Vdata] = PiecewiseLinearInterpolate(Potential,Range,NumLines);

    x_interior = x_interp(2:end-1);

    M = cell(1,NumLines);
    Minv = cell(1,NumLines);
    for i = 1:NumLines 
        M{i} = @(x,e) airyMatrix(Vdata{i}.slope,x,Vdata{i}.ystart,Vdata{i}.xstart,e,Mass);
        Minv{i} = @(x,e) airyInverseMatrix(Vdata{i}.slope,x,Vdata{i}.ystart,Vdata{i}.xstart,e,Mass);
    end

    T = zeros(1,length(Energy));
    for i = 1:length(Energy)
        
        %fprintf('E = %4.4f\n',Energy(i))

        e = Energy(i);

        kL = sqrt(2*Mass*(e-V_Left));
        kR = sqrt(2*Mass*(e-V_Right));


        % left end x = x_interp(1)
        ML = [  exp(1i*kL*x_interp(1)),       exp(-1i*kL*x_interp(1));
              1i*kL*exp(1i*kL*x_interp(1)), -1i*kL*exp(-1i*kL*x_interp(1))];

        % right end at x = x_interp(end)
        MR = [exp(1i*kR*x_interp(end)-x_interp(1)),     0;
             1i*kR*exp(1i*kR*x_interp(end)-x_interp(1)), 0];

        Mtr = eye(2);
        for j = 1:NumLines-1
            Mtr = Mtr * (Minv{j}(x_interior(j)-x_interp(1),e)*M{j+1}(x_interior(j)-x_interp(1),e));
        end
        Mtot = (ML\M{1}(x_interp(1),e)) * Mtr * (M{end}(x_interp(end)-x_interp(1),e)\MR);

        T(i) = 1/Mtot(1,1) * kR / kL;
        T(i) = conj(T(i))*T(i);

    end

    subplot(211)
    semilogy(Energy,T,'Linewidth',2)
    xlabel('E / Hartree')
    ylabel('Transmission |T|^2')
    axis([-Inf,Inf,0,5])
    set(gca,'FontSize',14,'Linewidth',2,'Box','off')
    grid on
       
    subplot(212)
    domain_length = abs(Range(2) - Range(1));
    hold on
    x0 = linspace(x_interp(1)-0.15*domain_length, x_interp(1), 20);
    plot(x0,V_Left*ones(1,length(x0)),'k-','Linewidth',3)
    for i = 1:length(Vdata)
        y = linspace(x_interp(i),x_interp(i+1),20);
        v1 = Vdata{i}.slope*y + Vdata{i}.intercept;
        v2 = Potential(y);
        plot(y,v1,'r-','Linewidth',2)
        plot(y,v2,'k-','Linewidth',1)
    end
    xN = linspace(x_interp(end),x_interp(end)+0.5*domain_length,20);
    plot(xN,V_Right*ones(1,length(xN)),'k-','Linewidth',3)
    plot(x_interp,V_interp,'ro')
    hold off
    xlabel('x / a.u.')
    ylabel('Potential / Hartree')
    set(gca,'FontSize',14,'Linewidth',2,'Box','off')
    Vmax = max([V_Right,V_Left,V_interp]);
    Vmin = min([V_Right,V_Left,V_interp]);
    potential_length = abs(Vmax-Vmin);
    axis([-Inf,Inf,Vmin-0.10*potential_length,Vmax+0.10*potential_length])
    grid on

end

% Functions
function [M] = airyMatrix(alpha,x,y_start,x_start,e,m)

    if abs(alpha) > eps

        l0 = airylengthscale(alpha,m);

        z = @(x) 1/l0*(x - x_start - (e - y_start)/alpha);

        M = [airy(0,z(x)), airy(2,z(x));
             1/l0*airy(1,z(x)), 1/l0*airy(3,z(x))];
         
    else
        
        k = sqrt(2*m*(e-y_start));
        
        M = [exp(1i*k*x), exp(-1i*k*x);
            1i*k*exp(1i*k*x), -1i*k*exp(-1i*k*x)];
        
    end
     
end

function [Minv] = airyInverseMatrix(alpha,x,y_start,x_start,e,m)

    if abs(alpha) > eps
        
        W = 1/pi;
        
        l0 = airylengthscale(alpha,m);

        z = @(x) 1/l0*(x - x_start - (e - y_start)/alpha);

        Minv = l0/(W) * [1/l0*airy(3,z(x)), -airy(2,z(x));
                       -1/l0*airy(1,z(x)),  airy(0,z(x))];
         
    else
        
        k = sqrt(2*m*(e-y_start));
        
        Minv = 1/(2i*k)*[-1i*k*exp(-1i*k*x), -exp(-1i*k*x);
                      -1i*k*exp(1i*k*x), exp(1i*k*x)];
        
    end
     
end

function [l0] = airylengthscale(alpha,m)
    
    l0 = zeros(1,length(alpha));
    for i = 1:length(alpha)
        l0(i) = (2*m*abs(alpha(i))).^(-1/3);
        if alpha(i) < 0
            l0(i) = -l0(i);
        end
    end
    
end

function [Ainv] = svd_inverse(A,thresh)

    [U,s,V] = svd(A,'econ'); s = diag(s);
    idx = find(s > thresh);
    U = U(:,idx); V = V(:,idx); s = s(idx);
    Ainv = V*diag(s.^-1)*U';
        
end

