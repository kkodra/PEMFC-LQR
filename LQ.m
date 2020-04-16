function LQ()

set(0,'DefaultTextFontName','Times','DefaultTextFontSize',18,...
'DefaultAxesFontName','Times','DefaultAxesFontSize',18,...
'DefaultLineLineWidth',2,'DefaultLineMarkerSize',7.75)

% Get augmented model
[augSys,eig_Aug,epsilon,TS_size] = PEMFC_FPS_Model;

% Get Schur ordered system
[T_schur_ord,ordSys] = ordered_Schur(augSys);

A = ordSys.A;
B = ordSys.B;
C = ordSys.C;

% Input due to blower and valve
B_blow = B(:,1);
B_valve = B(:,2);

% Select input
blow = false;

% Weights (design)
weights = cell(numel(TS_size),2);
weights{1,1} = diag([3613,7417,7059,7009,62,3743,9015,3183,5971]);
weights{2,1} = diag([8223940,251505,4144289,7314075,7813740,3672859]);
weights{3,1} = diag([131830,123500,190903]);
weights{1,2} = 298;
weights{2,2} = 1;
weights{3,2} = 15;

% Run time and input
t = cell(1,numel(TS_size)); u = cell(1,numel(TS_size));
t{1} = 0:0.01:15;      u{1} = 1*ones(1,length(t{1}));
t{2} = 0:0.001:0.2;    u{2} = 7*ones(1,length(t{2}));
t{3} = 0:0.001:0.05;   u{3} = 0.5*ones(1,length(t{3}));

% State to be plotted
state_ind = [4 4 2];

% Y labels for the plots
y_lab = {'Mass O_2 [kg]','Gas Pressure in SM [kPa]','Mass N_2 [kg]'};

for i = 1:3
    if i == 1
        ind_start = 1;
    end
    
    ind_end = ind_start+TS_size(i)-1;
    A_1 = A(ind_start:ind_end,ind_start:ind_end);   % DO NOT START AT 1 EVERY ITERATION!
    B_1 = B(ind_start:ind_end,:);
    C_1 = zeros(1,TS_size(i)); C_1(state_ind(i)) = 1; 
    D_1 = 0;
    
    if blow    
        B_1 = B_1(:,1);
    else
        B_1 = B_1(:,2);
    end
    
    sys_ol = ss(A_1,B_1,C_1,D_1);

    [P,~,~] = care(A_1,B_1,weights{i,1},weights{i,2});
    K_opt = inv(weights{i,2})*B_1'*P;  

    Nbar = rscale(sys_ol,K_opt);
    B_1_K_opt = B_1*K_opt;

    A_1_FB = A_1 - B_1_K_opt;

    sys = ss(A_1_FB,B_1,C_1,D_1);
    % Add Nbar function
    [~,~,X_cl] = lsim(sys,Nbar*u{i},t{i});
    
    figure(i)
    plot(t{i},X_cl(:,state_ind(i))); hold on
    xlabel('Time (s)');grid on
    ylabel(y_lab{i})
    
    ind_start = sum(TS_size(1:i)) + 1;

end
