function [HC_0] = build_HC_0(c1a,c1b,sys)
    
    HC_0 = einsum_kg(sys.fa_ov,c1a,'ia,ai->') + einsum_kg(sys.fb_ov,c1b,'ia,ai->');

end

