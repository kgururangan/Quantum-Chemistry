function [Lvec,resid,cc_t] = lefteomuccsd(omega,Rvec,HBar_t,cc_t,sys,opts)

    if opts.solver == 1
        [Lvec,resid,cc_t] = leftcc_solver1(omega,Rvec,HBar_t,cc_t,sys,opts);
    else
        [Lvec,resid,cc_t] = leftcc_solver2(omega,Rvec,HBar_t,cc_t,sys,opts);
    end

end

