% Main script to run algorithms
% (Work in progress)

clear
clc

% Import model
[sys_aug, sys_SP, sys_eig, epsilon, dim] = PEMFC_FPS_Model;

% Create a decoupled model using Chang
[slow_sys, fast_sys, LH_test, L, H] = decouple_sys(A,B,C,dim,epsilon)

% Create a SP model using developed algorithm
% TODO

% Controller design
% TODO
