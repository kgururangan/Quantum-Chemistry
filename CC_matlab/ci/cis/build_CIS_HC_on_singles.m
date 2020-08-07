function [sigma] = build_CIS_HC_on_singles(c1,sys)

    HC1 = einsum_kg(sys.e1int(

    sigma = HC1;
    
end

