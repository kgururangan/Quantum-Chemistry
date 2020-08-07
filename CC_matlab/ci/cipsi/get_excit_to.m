function [excit_unocc] = get_excit_to(Dref,D)
    
    [excit_unocc, ~] = setdiff(D,Dref,'stable');
    
end

