function [energy] = lucc_energy(L,LH,cc_t,sys)
    energy = sqrt(sum(LH.^2))./sqrt(sum(L.^2)) + ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);
end

