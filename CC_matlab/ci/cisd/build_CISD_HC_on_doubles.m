function [sigma] = build_CISD_HC_on_doubles(c1,c2,sys)

    HC1 = einsum_kg(sys.Vvvvo,c1,'abej,ei->abij')-einsum_kg(sys.Vvvvo,c1,'abei,ej->abij')...
          -einsum_kg(sys.Vovoo,c1,'mbij,am->abij')+einsum_kg(sys.Vovoo,c1,'maij,bm->abij');
    
    HC2 = einsum_kg(sys.fvv,c2,'ae,ebij->abij')-einsum_kg(sys.fvv,c2,'be,eaij->abij')...
          -einsum_kg(sys.foo,c2,'mi,abmj->abij')+einsum_kg(sys.foo,c2,'mj,abmi->abij')...
          +0.5*einsum_kg(sys.Voooo,c2,'mnij,abmn->abij')...
          +0.5*einsum_kg(sys.Vvvvv,c2,'abef,efij->abij')...
          +einsum_kg(sys.Vvoov,c2,'amie,ebmj->abij')...
          -einsum_kg(sys.Vvoov,c2,'bmie,eamj->abij')...
          -einsum_kg(sys.Vvoov,c2,'amje,ebmi->abij')...
          +einsum_kg(sys.Vvoov,c2,'bmje,eami->abij');

    sigma = HC1 + HC2;

end

