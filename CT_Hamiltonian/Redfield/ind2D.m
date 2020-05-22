function J = ind2D(Ne,a,b)

% returns row/col index of superoperatr matrix Ne^2 x Ne^2

va = zeros(1,Ne);
va(a) = 1;
vb = zeros(1,Ne);
vb(b) = 1;

Vtemp = kron(va,vb); [~,J] = find(Vtemp == 1);

end