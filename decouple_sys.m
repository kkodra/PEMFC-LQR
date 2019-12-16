function [slow_sys, fast_sys, LH_test, L, H] = decouple_sys(A,B,C,dim,epsilon)

% TODO: Create loop. Run len(dim) times and create subsystems
% NEEDS modification to agree with algorithm in paper

% From bottom to top [1 4]
A1 = A(1:dim,1:dim);
A2 = A(1:dim,dim+1:end);
A3 = A(dim+1:end,1:dim);
A4 = A(dim+1:end,dim+1:end);

B1 = B(1:dim,:);
B2 = B(dim+1:end,:);

C1 = C(:,1:dim);
C2 = C(:,dim+1:end);

[L,H] = eval_L_H(A1,A2,A3,A4,epsilon,'recursive');

As = (A1 - A2*L);
Af = (A4 + epsilon*L*A2);

Bs = (B1 - H*B2 - epsilon*H*L*B1);
Bf = (B2 + epsilon*L*B1);

Cs = C1 - C2*L;
Cf = C2 - epsilon*C2*L*H + epsilon*C1*H;

Test1 = A4*L-epsilon*L*A1+epsilon*L*A2*L-A3;
Test2 = epsilon*(A1 - A2*L)*H-H*(A4 +epsilon*L*A2)+A2;

%% Create structs
% Slow and fast subsystems
slow_sys.A = As;
slow_sys.B = Bs;
slow_sys.C = Cs;

fast_sys.A = Af;
fast_sys.B = Bf;
fast_sys.C = Cf;

% L and H test (see if they are correct; Test1 and TEst2 should evaluate to 0)
LH_test.Test1 = Test1;
LH_test.Test2 = Test2;

end

