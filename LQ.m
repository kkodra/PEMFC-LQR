% TO BE MODIFIED

'DefaultTextFontSize',18,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',18,...
'DefaultLineLineWidth',2,...
'DefaultLineMarkerSize',7.75)
A = A_schur_ord7;
B = T_schur_ord7*B_aug;
C = C_aug*inv(T_schur_ord7);

% Weights are predetermined. Add them here.

% Input due to blower
B_blow = B(:,1);
B_valve = B(:,2);

% Select B
t = 0:0.01:10;
u = 7*ones(1,length(t));

% Full system
Q = eye(size(A));
[P,Lcare,G] = care(A,B_blow,Q,R);
K_opt = inv(R)*B_blow'*P;
BK_opt = B_blow*K_opt;

A_FB = A - BK_opt;

sys = ss(A_FB,B_blow,C,0);

[Y,~,X] = lsim(sys,u,t); 

%% Slowest subsystem
A_1 = A(1:9,1:9);
B_1 = B(1:9,:);
C_1 = C(:,1:9);

B_1_blow = B_1(:,1);
B_1 = B_1_blow;

% Full system
Q = eye(size(A_1));
[P,Lcare,G] = care(A_1,B_1,Q,R);
K_opt = inv(R)*B_1'*P;
B_1_K_opt = B_1*K_opt;

A_1_FB = A_1 - B_1_K_opt;

sys = ss(A_1_FB,B_1,C_1,0);

[Y,~,X] = lsim(sys,u,t);
