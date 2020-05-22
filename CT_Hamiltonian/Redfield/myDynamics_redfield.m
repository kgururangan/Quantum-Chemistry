function TM = myDynamics_redfield(par,omega_T,pos,Redfield_pop)

%GammaPop = par.GammaPop;
Gammadab = par.Gammadab;
tau_3 = par.tau_3;
omega_g = par.omega_g;
bw_Elec = par.bw_Elec;
n_quanta = par.n_quanta;
Ne = par.Ne;

%omega_T = myBuiltOmegaT(omega_g,n_quanta,bw_Elec);

GammaM = zeros(Ne^2,Ne^2);

for a = 1:Ne
    for b = 1:Ne
        
        J = ind2D(Ne,a,b);
        
        if a ~= b
            wab = omega_T(a) - omega_T(b);
            GammaM(J,J) = -1i*wab - Gammadab;
        else
            for k = 1:Ne
                K = ind2D(Ne,k,k);
                GammaM(J,K) = Redfield_pop(a,k);
                %GammaM(J,K) = GammaPop(k,a);
                %GammaM(J,K) = Redfield_pop(J,K);
            end
        end
        
    end
end

N3 = length(tau_3);
xout = zeros(N3,size(GammaM,1),size(GammaM,2));
TM = zeros(Ne^4,length(tau_3));

for i = 1:N3
    Ev = expm(GammaM*tau_3(i));
    xout(i,:,:) = Ev;
    TM(:,i) = Ev(:);
end

if ~isempty(pos)
    TM = TM(pos,:);
end

end