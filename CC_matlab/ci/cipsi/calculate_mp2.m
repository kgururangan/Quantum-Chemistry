function [eMP2] = calculate_mp2(Alpha,Wavefunction,sys)

    eMP2 = 0.0;
    
    n_internal = size(Wavefunction.Dets,1);
    
    for I = 1:n_internal
        for J = 1:n_internal
            
            [val1] = slater_eval(sys.ZM,sys.VM,Alpha,Wavefunction.Dets(I,:));
            [val2] = slater_eval(sys.ZM,sys.VM,Wavefunction.Dets(J,:),Alpha);
            [e_alpha] = slater_eval(sys.ZM,sys.VM,Alpha,Alpha);
            
            val = Wavefunction.Coef(I)*Wavefunction.Coef(J)*val1*val2;
            
            eMP2 = eMP2 + val/(Wavefunction.Energy - e_alpha);
        end
    end
    
end