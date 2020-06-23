clear all
 
Nunocc = 6;
Nocc = 3;

Gamma = zeros(Nunocc,Nunocc,Nunocc);

i = 1; j = 2; k = 3;

A = rand(Nocc,Nunocc);
B = rand(Nocc,Nocc,Nunocc,Nunocc);

% for a = 1:Nunocc
%     for b = 1:Nunocc
%         for c = 1:Nunocc
%             Gamma(a,b,c) = A(i,a)*B(j,k,b,c);
%         end
%     end
% end

Gamma = 0.0;
for i = 1:Nocc
    for a = 1:Nunocc
        Gamma = Gamma + A(i,a)*B(1,i,1,a);
    end
end

Gamma2 = einsum_kg(A,B(1,:,1,:),'ia,ia->');

Gamma-Gamma2

%%

s = 'a,a->';

%split input and output indices
s=split(s,'->');

%split input indices
in=s{1};
out=s{2};

in=split(in,',');

idxA=in{1};
idxB=in{2};

final_permutation = [];
iA_con = [];
iB_con = [];

for i=1:length(idxA)
    
    j=find(idxB==idxA(i));
    
    if isempty(j)   % i is an output index
        
        j=find(out==idxA(i));
        final_permutation(end+1)=j;
        
    else            % i is contracted
        
        iA_con(end+1)=i; 
        iB_con(end+1)=j; 
        
    end
end

for i=1:length(idxB)
    
    j=find(idxB(i)==out);
    
    if ~isempty(j)   % i is an output index
        final_permutation(end+1)=j;
    end
    
end


[~, final_permutation]=sort(final_permutation);