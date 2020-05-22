function [Smat, Tmat, Vmat, VVmat] = spatial_integrals_v2(orbs,Rat0,Z)

Norb = length(orbs);
Smat = zeros(Norb); Vmat = zeros(Norb); Tmat = zeros(Norb); VVmat = zeros(Norb,Norb,Norb,Norb);

for alpha = 1:Norb
    for beta = 1:Norb
        
        gA = orbs{alpha}; gB = orbs{beta};
        
        for i = 1:length(gA.coeff)
            for j = 1:length(gB.coeff)
                cA = gA.coeff(i); cB = gB.coeff(j);
                xiA = gA.exps(i); xiB = gB.exps(j);
                rA = gA.origin; rB = gB.origin;
                
                Smat(alpha,beta) = Smat(alpha,beta) + cA*cB*sAB(xiA,rA,xiB,rB);
                Tmat(alpha,beta) = Tmat(alpha,beta) + cA*cB*tAB(xiA,rA,xiB,rB);
                
                for k = 1:length(Z)
                    Vmat(alpha,beta) = Vmat(alpha,beta) + cA*cB*vAB(xiA,rA,xiB,rB,Z(k),Rat0(k,:));
                end
                
                for gamma = 1:Norb
                    for delta = 1:Norb
                        
                        gC = orbs{gamma}; gD = orbs{delta};
                        
                        for l = 1:length(gC.coeff)
                            for m = 1:length(gD.coeff)
                                cC = gC.coeff(l); cD = gD.coeff(m);
                                xiC = gC.exps(l); xiD = gD.exps(m);
                                rC = gC.origin; rD = gD.origin;
                                
                                VVmat(alpha,beta,gamma,delta) = VVmat(alpha,beta,gamma,delta) +...
                                                                cA*cB*cC*cD*vABCD(xiA,rA,xiB,rB,xiC,rC,xiD,rD);
                            end
                        end
                        
                    end
                end
                
            end
        end
        
    end
    
end

VVmat = permute(VVmat,[1,3,2,4]);


end