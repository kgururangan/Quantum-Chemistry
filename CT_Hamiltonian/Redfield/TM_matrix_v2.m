function [TM,Ev] = TM_matrix_v2(Kt,tau_3,trans,type)

% Kt is an input relaxation supermatrix. Kt = Ne^2 x Ne^2

if trans == 1
    Kt = Kt';
end

N = length(Kt); % N = Ne^2

TM = zeros(N^2,length(tau_3));

if ~strcmp(type,'RK4')
    
Ev = zeros(N,N,length(tau_3));

for i = 1:length(tau_3)
    %if isempty(V_S)
    A = expm(Kt*tau_3(i));
    Ev(:,:,i) = A;
    %else
    %    VSS = kron(V_S,V_S);
    %    VSST = kron(V_S',V_S');
    %    Ev = expm(VSS*Kt*VSST*tau_3(i) );
    %end
    TM(:,i) = A(:);
end

else
    
    % if RK4, Kt should be N x N x N3 matrix
    
    
    %R = reshape(Red_exc,N^2,tau_3);
    
    L = @(ti) interp1(tau_3,Kt,ti);
    
    Nteval = 2*N3 + 1;
    t0 = tau_3(1);
    tf = tau_3(end);
    dt_eval = 1/Nteval*(tf-t0);
    teval = 0:dt_eval:(Nteval-1)*dt_eval;
    x0 = [0,0,0,1]';
    xeval = zeros(length(x0),Nteval);
    
    for i = 1:Nteval
        
        % Forward Euler implicit
        %x0 = x0 + dt_ps*L(t(i))*x0;
        %rho(:,i) = x0;

        % Fourth-order Runge-Kutta (RK4)
        k1 = dt*L(teval(i))*x0;
        k2 = dt*L(teval(i)+0.5*dt)*(x0 + 0.5*k1);
        k3 = dt*L(teval(i)+0.5*dt)*(x0+0.5*k2);
        k4 = dt*L(teval(i)+dt)*(x0+k3);
        
        x0 = x0 + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        
        xeval(i,:) = x0;

    end

end

