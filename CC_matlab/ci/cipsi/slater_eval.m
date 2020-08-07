function [onebody, twobody] = slater_eval(z, v, det1p, det2p, sign, excit_rank)


%     [det1p, det2p, sign, excit_rank] = reorder_max_coincidence(det1, det2);
    
    onebody = 0.0; twobody = 0.0;

    if excit_rank == 0
        for i = 1:length(det1p)
            onebody = onebody + sign*z(det1p(i), det1p(i));
            for j = 1:length(det1p)
                twobody = twobody + 0.5*sign*( v(det1p(i),det1p(j),det1p(i),det1p(j)) -...
                                          v(det1p(i),det1p(j),det1p(j),det1p(j)) );
            end

        end
    elseif excit_rank == 1
        onebody = sign*z(det1p(1),det2p(1));
        for i = 1:length(det1p)
            twobody = twobody + sign*v(det1p(1),det1p(i),det2p(1),det2p(i)) -...
                                sign*v(det1p(1),det1p(i),det2p(i),det2p(1));
        end
    elseif excit_rank == 2
        twobody = sign*v(det1p(1),det1p(2),det2p(1),det2p(2)) - ...
                  sign*v(det1p(1),det1p(2),det2p(2),det1p(1));
    else
        twobody = 0.0; onebody = 0.0;
    end
    
        
        

end
      