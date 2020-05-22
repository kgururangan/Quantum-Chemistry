function [states] = get_states(num_sites)

    % states are represented by occupation number vectors with length
    % equal to the number of sites in the model
    % i.e. this is a minimum basis set problem (spinless) and we are 
    % building the FCI space

    max_excit = 1; % |0> or |1>
    num_basis = (max_excit+1)^num_sites;
    states = zeros(num_sites, num_basis);
    
    excits = 0:max_excit;
    
    odo = zeros(1,num_sites);
    for i = 1:num_basis
        states(:,i) = excits(odo + 1);
        for j = num_sites:-1:1
            odo(j) = odo(j) + 1;
            if odo(j) == max_excit+1
                odo(j) = 0;
            else
                break;
            end
        end
    end

end

