function [C] = einsum_kg(A, B, contract_scheme)

% we put this to let us perform contractions between sliced tensors where
% the slicing could result in singleton dimensions, i.e.
% einsum_kg(v(:,:,2,:)*t(:,:,:,6),'pqr,rqp->') as is used in CR-CC(2,3) codes.
% but it's dangerous when you e.g. have active orbitals that leave a singleton dimension. 
% we need to implement some fix that catches both cases e.g.
% squeeze only if # contraction indices = dim(A) - # singleton dimensions

% s=split(contract_scheme,'->');
% in=s{1}; out=s{2};
% 
% %split input indices
% in=split(in,',');
% idxA=in{1}; idxB=in{2};
% 
% if length(size(A)) ~= length(idxA)
%     A = squeeze(A);
% end
% if length(size(B)) ~= length(idxB)
%     B = squeeze(B);
% end

%A = squeeze(A); B = squeeze(B); 

szA=size(A);  szB=size(B);

[iA_con, iB_con, final_permutation] = parse_contraction(contract_scheme);

if size(iA_con)~=size(iB_con)
    error('number of dimensions to contract should be equal')
end

for i=1:length(iA_con)
    %fprintf('Contr %d: A(%d) with B(%d) - %d x %d\n', ...
    %         i,iA_con(i),iB_con(i),size(A,iA_con(i)),size(B,iB_con(i)))
    if size(A,iA_con(i))~=size(B,iB_con(i))
        error(['cannot contract dimension %d of 1st argument (length=%d)'...
            ' with dimension %d of 2nd argument (length=%d)'],...
            iA_con(i),size(A,iA_con(i)),iB_con(i),size(B,iB_con(i)))
    end
end
if length(iA_con)~=length(unique(iA_con)) || length(iB_con)~=length(unique(iB_con))
    error('each dimension should appear only once.')
end

dimA_unc=setdiff(1:length(szA),iA_con);
dimB_unc=setdiff(1:length(szB),iB_con);

A=permute(A, [dimA_unc iA_con]);
B=permute(B, [iB_con dimB_unc]);

A=reshape(A, [], prod(szA(iA_con)));
B=reshape(B, prod(szB(iB_con)), []);

C=A*B;

output_shape=[szA(dimA_unc),szB(dimB_unc)];

if length(output_shape)>1
    
    C = squeeze(reshape(C,[szA(dimA_unc),szB(dimB_unc)]));

    
    if final_permutation
        C=permute(C,final_permutation);
    end
    
end

end



function [iA_con, iB_con, final_permutation] = parse_contraction(s)

        %split input and output indices
        s=split(s,'->');
        in=s{1}; out=s{2};
        
        %split input indices        
        in=split(in,',');
        idxA=in{1}; idxB=in{2};

        final_permutation = [];
        iA_con = [];
        iB_con = [];

        for i=1:length(idxA)

            j=find(idxB==idxA(i));

            if isempty(j)   % i is an output index

                j=find(out==idxA(i));
                final_permutation(end+1)=j;%#ok<AGROW>

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

%     end
end

function C = OuterProduct(A, B)  % version 5
    C = reshape(A(:) * B(:).', [size(A), size(B)]);
end