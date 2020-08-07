function [val] = R(t,u,v,n,p,PCx,PCy,PCz,RPC)

    T = p*RPC*RPC;
    val = 0.0;
    if t == u && u == v && t == 0
        %boysval = boys_debug(n,T);   
        boysval = boys(n,T);
        val = val + power(-2*p,n)*boysval;
    elseif t == u && t == 0
        if v > 1
            val = val + (v-1)*R(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC);
        end
        val = val + PCz*R(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC);
    elseif t == 0
        if u > 1
            val = val + (u-1)*R(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC);
        end
        val = val + PCy*R(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC);
    else
        if t > 1
            val = val + (t-1)*R(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC);
        end
        val = val + PCx*R(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC);
    end
    
end

