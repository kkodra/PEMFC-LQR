function [As,Af,Bs,Bf,Cs,Cf,Test1norm,Test2norm,L,H] = DecoupleSys(A,B,C,DimTS,epsilon)

% If permutation matrix used
% A1 = A(length(A)-DimTS +1:end,length(A)-DimTS +1:end);
% A2 = A(length(A)-DimTS +1:end,1:length(A)-DimTS);
% A3 = A(1:length(A)-DimTS,length(A)-DimTS +1:end);
% A4 = A(1:length(A)-DimTS,1:length(A)-DimTS);

% From bottom to top [1 4]
A1 = A(DimTS +1:end,DimTS +1:end);
A2 = A(DimTS +1:end,1:DimTS);
A3 = A(1:DimTS,DimTS +1:end);
A4 = A(1:DimTS,1:DimTS);

B1 = B(DimTS +1:end,:);
B2 = B(1:DimTS,:);
% % 
C1 = C(:,DimTS +1:end);
C2 = C(:,1:DimTS);

[L,H] = eval_L_H(A1,A2,A3,A4,epsilon,'recursive');

As = (A1 - A2*L);
Af = (A4 + epsilon*L*A2);

Bs = (B1 - H*B2 - epsilon*H*L*B1);
Bf = (B2 + epsilon*L*B1);
% % 
Cs = C1 - C2*L;
Cf = C2 - epsilon*C2*L*H + epsilon*C1*H;

Test1 = A4*L-epsilon*L*A1+epsilon*L*A2*L-A3;
Test2 = epsilon*(A1 - A2*L)*H-H*(A4 +epsilon*L*A2)+A2;

Test1norm = norm(Test1);
Test2norm = norm(Test2);

end

