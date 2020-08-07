function [excit_occ] = get_excit_from(Dref,D)
    
    [excit_occ, ~] = setdiff(Dref,D,'stable');
    
end

