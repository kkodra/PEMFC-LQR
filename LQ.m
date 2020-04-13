% LQ Energies
set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',18,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',18,...
'DefaultLineLineWidth',2,...
'DefaultLineMarkerSize',7.75)

% Get augmented model
[augSys,eig_Aug,epsilon,TS_size] = PEMFC_FPS_Model;

% Get Schur ordered system
[T_schur_ord,ordSys] = ordered_Schur(augSys);

A = ordSys.A;
B = ordSys.B;
C = ordSys.C;

% Input due to blower
B_blow = B(:,1);
B_valve = B(:,2);

% Select B
t = 0:0.01:15;
u = 1*ones(1,length(t));

%% Slowest subsystem

load('Q1.mat')
load('R1.mat')
A_1 = A(1:9,1:9);
B_1 = B(1:9,:);
C_1 = [0 0 0 1 0 0 0 0 0]; %C(4,1:9);
D_1 = zeros(1);

B_1_blow = B_1(:,1);
B_1 = B_1_blow;

sys_ol = ss(A_1,B_1,C_1,0);

[P,~,~] = care(A_1,B_1,Q1,R1);
K_opt = inv(R1)*B_1'*P;  

Nbar = rscale(sys_ol,K_opt);
B_1_K_opt = B_1*K_opt;

A_1_FB = A_1 - B_1_K_opt;

sys = ss(A_1_FB,B_1,C_1,D_1);

[~,~,X_cl] = lsim(sys,Nbar*u,t);

figure(1)
plot(t,X_cl(:,4)); hold on% mass of oxygen
xlabel('Time (s)');grid on
ylabel('Mass O_2 [kg]')


%% 2nd subsystem
% 
load('Q2.mat')
load('R2.mat')

t = 0:0.001:0.2;
u = 7*ones(1,length(t));

A_2 = A(10:15,10:15);
B_2 = B(10:15,:);
C_2 =[0 0 0 1 0 0]; % C(:,10:15);

B_2_blow = B_2(:,1);
B_2_valve = B_2(:,2);
B_2 = B_2_valve;

sys_ol = ss(A_2,B_2,C_2,0);

[P,~,~] = care(A_2,B_2,Q2,R2);
K_opt = inv(R2)*B_2'*P;
B_2_K_opt = B_2*K_opt;

Nbar = rscale(sys_ol,K_opt);

A_2_FB = A_2 - B_2_K_opt;

sys = ss(A_2_FB,B_2,C_2,0);

[~,~,X] = lsim(sys,Nbar*u,t);

figure(2)
plot(t,X(:,4)); hold on% pressure gas supply manifold
xlabel('Time (s)');grid on
ylabel('Gas Pressure in SM [kPa]')

%% Fastest subsystem
t = 0:0.001:0.05;
u = 0.5*ones(1,length(t));
A_3 = A(16:end,16:end);
B_3 = B(16:end,:);
C_3 = [0 1 0]; %C(:,16:end);
% 
B_3_blow = B_3(:,1);
B_3 = B_3_blow;

load('Q3.mat')
load('R3.mat')

sys_ol = ss(A_3,B_3,C_3,0);

[P,Lcare,G] = care(A_3,B_3,Q3,R3);
K_opt = inv(R3)*B_3'*P;
B_3_K_opt = B_3*K_opt;

Nbar = rscale(sys_ol,K_opt);
A_3_FB = A_3 - B_3_K_opt;

sys = ss(A_3_FB,B_3,C_3,0);

[~,~,X] = lsim(sys,Nbar*u,t);

figure(3)
plot(t,X(:,2)); hold on% mass of nitrogen
xlabel('Time (s)');grid on
ylabel('Mass N_2 [kg]')

