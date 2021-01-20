function [Lvec,resid,cc_t,omega_left] = lefteomuccsd(omega,Rvec,HBar_t,cc_t,sys,opts)

    if opts.solver == 1
        [Lvec,resid,cc_t,omega_left] = leftcc_solver1(omega,Rvec,HBar_t,cc_t,sys,opts);
    else
        [Lvec,resid,cc_t,omega_left] = leftcc_solver2(omega,Rvec,HBar_t,cc_t,sys,opts);
    end

end

