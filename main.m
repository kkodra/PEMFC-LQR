% Main script to run algorithms
% (Work in progress)

clear
clc

% Import model
[sys_aug, sys_SP, sys_eig, epsilon, dim] = PEMFC_FPS_Model;

% SP system matrices
A = sys_SP.A_SP;
B = sys_SP.B_SP;
C = sys_SP.C_SP;

% Create a decoupled model using Chang
[slow_sys, fast_sys, LH_test, L, H] = decouple_sys(A,B,C,dim,epsilon)

% Ensure L and H have been solved correctly. Set a threshold
if norm(LH_test.Test1) > 10e-6 || norm(LH_test.Test2) > 10e-6 
	error('L or H are not correct.')
end

% Create a SP model using developed algorithm
% TODO

% Controller design
% TODO
