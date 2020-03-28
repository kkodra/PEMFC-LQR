% LQ Energies

set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',18,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',18,...
'DefaultLineLineWidth',2,...
'DefaultLineMarkerSize',7.75)
A = A_schur_ord7;
B = T_schur_ord7*B_aug;
C = C_aug*inv(T_schur_ord7);

R = 1;

% Input due to blower
B_blow = B(:,1);
B_valve = B(:,2);

% Select B
t = 0:0.01:10;
u = ones(1,length(t));

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

plot(t,X(:,3)); hold on  % speed of blower
plot(t,X(:,4)); % mass of oxygen
plot(t,X(:,8)) % pressure on anode
legend('\omega_{blower}','m_{O_2}','p_{an}')
xlabel('Time (s)');grid on
ylabel('Amplitude')


%% 2nd subsystem
A_2 = A(10:15,10:15);
B_2 = B(10:15,:);
C_2 = C(:,10:15);

B_2_blow = B_2(:,1);
B_2_valve = B_2(:,2);
B_2 = B_2_valve;

% Full system
Q = eye(size(A_2));
[P,Lcare,G] = care(A_2,B_2,Q,R);
K_opt = inv(R)*B_2'*P;
B_2_K_opt = B_2*K_opt;

A_2_FB = A_2 - B_2_K_opt;

sys = ss(A_2_FB,B_2,C_2,0);

[Y,~,X] = lsim(sys,u,t);

figure
% plot(t,X(:,2));   
plot(t,X(:,2));hold on %pressure return manifold
plot(t,X(:,3));  % mass of water anode
%plot(t,X(:,5)) % pressure hydrogen anode
%plot(t,X(:,6)) % compressor speed
legend('p_{rm}','m_{w,an}')
xlabel('Time (s)');grid on
ylabel('Amplitude')
xlim([0 1])
%% Fastest subsystem
A_3 = A(16:end,16:end);
B_3 = B(16:end,:);
C_3 = C(:,16:end);

B_3_blow = B_3(:,1);
B_3 = B_3_blow;

% Full system
Q = eye(size(A_3));
[P,Lcare,G] = care(A_3,B_3,Q,R);
K_opt = inv(R)*B_3'*P;
B_3_K_opt = B_3*K_opt;

A_3_FB = A_3 - B_3_K_opt;

sys = ss(A_3_FB,B_3,C_3,0);

[Y,~,X] = lsim(sys,u,t);

figure
%plot(t,X(:,1)); hold on  % Tcpox
plot(t,X(:,2)); % mass of nitrogen
% plot(t,X(:,3)); % mass air in supply manifold
ylim([0.01 0.017])
xlim([0 1])
legend('m_{N_2}')
xlabel('Time (s)');grid on
ylabel('Amplitude')
