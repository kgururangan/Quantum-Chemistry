function [energy] = lcc_energy(L,LH,t1,t2,sys)

% LCC energy is given by
% E_lcc = <0|LH|0> = <0|L_cl H_cl|0> + <0|L_op H_cl|0> + ...
%                    <0|L_cl H_op|0> + <0|L_op H_op|0>
%                  = <0|L_cl H_cl|0> + <0|L_op H_op|0>
%                  = 1*<0|H_cl|0> + <0|L_op H_op|0>

% <0|H_cl|0> = E_cc -> the ground-state CC correlation energy
% <0|L_op H_op = LH vector we are computing in LCC
% At convergence, LH = omega*L and and since |L| != 1 usually, 
% the measure |LH|/|L| + E_corr = omega + E_corr for any state, thus lcc energy recovers
% CC energy upon convergence.

    energy = sqrt(sum(LH.^2))./sqrt(sum(L.^2)) + cc_energy(t1,t2,sys);

end

