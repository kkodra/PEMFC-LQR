function [Tordered,U] = OrderedSchur(A)
%Calculate ordered Schur form of a real non-symmetric matrix.
%   Ref: L. K. Balyan, "An Algorithm of Ordered Schur
%   Factorization for Real Nonsymmetric Matrix," World Aca. of Sci., 2010.

[U,T] = schur(A);

location = 0;
diagonal = diag(T);
for i = 1:(length(diagonal) - 1)
    if diagonal(i)>diagonal(i+1)
        location = i;
    end
end 
    
T11 = T(location,location);
T12 = T(location,location+1);
T22 = T(location+1,location+1);

Tswap = [T11 T12;0 T22];
X = lyap(T11,-T22,-T12);

qr_Matrix = [-X;1];

[Q,~] = qr(qr_Matrix);

Tnew = Q'*Tswap*Q;

T(location:location+1,location:location+1) = Tnew;
Tordered = T;
end

