function [e1int, e2int, Vnuc, Norb] = load_integrals(varargin)

    % varargin{1} = work_path where onebody.inp and twobody.inp are located
    % varargin{2} = Number of spatial orbitals

    if length(varargin) == 1
        %path = '/Users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/';
        onebody = [varargin{1},'onebody.inp'];
        twobody = [varargin{1},'twobody.inp'];

    %     if length(varargin) == 2
    %         [e1int, Norb] = read_onebody_ints(onebody,varargin{2});
    %     else
        [e1int,Norb] = read_onebody_ints(onebody);
    %     end
        [e2int, Vnuc] = read_twobody_ints(twobody, Norb);
    else
        onebody = varargin{1}; twobody = varargin{2};
        [e1int,Norb] = read_onebody_ints(onebody);
        [e2int, Vnuc] = read_twobody_ints(twobody, Norb);
    end
    
end



