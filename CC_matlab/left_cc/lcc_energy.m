function [energy] = lcc_energy(LH,t1,t2,VM,FM,occ,unocc)

% LCC energy is given by
% E_lcc = <0|LH|0> = <0|L_cl H_cl|0> + <0|L_op H_cl|0> + ...
%                    <0|L_cl H_op|0> + <0|L_op H_op|0>
%                  = <0|L_cl H_cl|0> + <0|L_op H_op|0>
%                  = 1*<0|H_cl|0> + <0|L_op H_op|0>

% <0|H_cl|0> = E_cc -> the ground-state CC correlation energy
% <0|L_op H_op = LH vector we are computing in LCC
% At convergence, LH = omega*L and |L| = 1 (NOT TRUE!!!), therefore the measure
% |LH| + E_corr = E_corr + omega for any state, thus lcc energy recovers
% CC energy upon convergence.

    energy = sqrt(sum(LH.^2)) + cc_energy(t1,t2,VM,FM,occ,unocc);

end

