% Main script to run algorithms

clear
clc

% Import model
[sys_aug, sys_eig] = PEMFC_FPS_Model;

% SP system matrices
A = sys_SP.A;
B = sys_SP.B;
C = sys_SP.C;

% Plot magnitude of eigenvalues of the system
figure
bar(sort(abs(sys_eig)))w grid on
xlabel('State','interpreter','latex');
ylabel('$|\lambda_i|$','interpreter','latex')


% Create a SP model using developed algorithm
% Add Schur decomposed model
[T_ordered, ordSys] = ordered_Schur(augSys)

A_schur_ord = ordSys.A;
B_schur_ord = ordSys.B;
C_schur_ord = ordSys.C;

% TODO: Maybe a loop to get N different time-scales that the system has
% See algorithm in paper
% COMPLETE THIS!
num_TS = 3;
for i = 1:num_TS-1
	[slow_sys, fast_sys, LH_test, L, H] = decouple_sys(A_schur_ord,B_schur_ord,C_schur_ord,dim,epsilon)
end

% Create a decoupled model using Chang
[slow_sys, fast_sys, LH_test, L, H] = decouple_sys(A_schur_ord,B_schur_ord,C_schur_ord,dim,epsilon)

% Ensure L and H have been solved correctly. Set a threshold
if norm(LH_test.Test1) > 10e-6 || norm(LH_test.Test2) > 10e-6
        error('L or H are not correct.')
end


%% Controller design
% Controllability matrices 
SP_cont = ctrb(A,B);
slow_cont = ctrb(slow_sys.A, slow_sys.B);
fast_cont = ctrb(fast_sys.A, fast_sys.B);

% LQR for slow and fast models
if rank(slow_cont) < rank (slow_sys.As)
	error('Controllability matrix of the slow subsystem is singular.');
elseif rank(fast_cont) < rank(fast_sys.Af)
	error('Controllability matrix of the fast subsystem is singular.');
end

% Define time and input functions
t = 0:0.01:10;
r = ones(size(t));

% LQR controller for overall augmented system
Q_aug = sys_aug.C'*sys_aug.C;
R_aug = 1;
[K_aug, sys_FB_aug, y_aug, x_aug] = LQR_control(Q_aug, R_aug, sys_aug, r, t)

% LQR controller for slow subsystem
Q_slow = slow_sys.C'*slow_sys.C;
R_slow = 1;

[K_slow, sys_FB_slow, y_slow, x_slow] = LQR_control(Q_slow, R_slow, sys_slow, r, t)

% LQR controller for fast subsystem
Q_fast = fast_sys.C'*fast_sys.C;
R_fast = 1;

[K_fast, sys_FB_fast, y_fast, x_fast] = LQR_control(Q_fast, R_fast, sys_fast, r, t)

% Plot results. TODO check & complete
% Slow subsystem state
figure
plot(t,x_slow(1,:));grid on
xlabel('Time (s)'); ylabel('States')

% Fast subsystem state
figure
plot(t,x_fast(1,:));grid on
xlabel('Time (s)'); ylabel('States')

%% Controller for Schur-decomposed system
% Obtain Schur decomposed system
[Tordered,U] = ordered_Schur(A)
A_schur = U;
B_schur = Tordered*B;
C_schur = C*Tordered;

% Get slow and fast subsystems
A_schur_slow = A_schur(1:dim,1:dim);
A_schur_fast = A_schur(dim+1:end,dim+1:end);

B_schur_slow = B_schur(1:dim,:);
B_schur_fast = B_schur(dim+1:end,:);

C_schur_slow = C_schur(:,1:dim);
C_schur_slow = C_schur(:,dim+1:end);
