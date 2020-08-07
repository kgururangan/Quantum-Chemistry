function [sigma] = build_CISD_HC_on_singles(c1,c2,sys)

    HC1 = einsum_kg(sys.fvv,c1,'ae,ei->ai')...
          -einsum_kg(sys.foo,c1,'mi,am->ai')...
          +einsum_kg(sys.Vvoov,c1,'amie,em->ai');

    HC2 = 0.5*einsum_kg(sys.Vvovv,c2,'anef,efin->ai')...
          -0.5*einsum_kg(sys.Vooov,c2,'mnif,afmn->ai')...
          +einsum_kg(sys.fov,c2,'me,aeim->ai');

    sigma = HC1 + HC2;
    
end

