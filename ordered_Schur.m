function [T_ordered, A_schur, B_schur, C_schur] = ordered_Schur(augSys)
% TODO: Automate the script. Solution correct but not efficient.

[T_schur,A_schur] = schur(A_aug);
% Current ordering according to the E below. Change indexing accordingly
E = ordeig(A_schur);
E_aug = [abs(E) [1:18]'];
E_aug_sort = sortrows(E_aug,1);

% ============================
% Add while loop to automate...
% ============================

select = zeros(1,18);
select(13:14) = 1;

% Swap 1
[T_schur_ord,A_schur_ord] = ordschur(T_schur,A_schur,select);

E1 = ordeig(A_schur_ord);
E1_aug = [abs(E1) [1:18]'];
E1_aug_sort = sortrows(E1_aug,1);

select = zeros(1,18);
select(1:2) = 1;
select(10) = 1;
select(17) = 1;

% Swap 2
[T_schur_ord1,A_schur_ord1] = ordschur(T_schur_ord,A_schur_ord,select);

E2 = ordeig(A_schur_ord1);
E2_aug = [abs(E2) [1:18]'];
E2_aug_sort = sortrows(E2_aug,1);

select = zeros(1,18);
select(1:4) = 1;
select(11) = 1;
select(16:17) = 1;

% Swap 3
[T_schur_ord2,A_schur_ord2] = ordschur(T_schur_ord1,A_schur_ord1,select);
E3 = ordeig(A_schur_ord2);
E3_aug = [abs(E3) [1:18]'];
E3_aug_sort = sortrows(E3_aug,1);

select = zeros(1,18);
select(1:7) = 1;
select(13) = 1;
select(18) = 1;

% Swap 4
[T_schur_ord3,A_schur_ord3] = ordschur(T_schur_ord2,A_schur_ord2,select);
E4 = ordeig(A_schur_ord3);
E4_aug = [abs(E4) [1:18]'];
E4_aug_sort = sortrows(E4_aug,1);

select = zeros(1,18);
select(1:9) = 1;
select(18) = 1;

% Swap 5
[T_schur_ord4,A_schur_ord4] = ordschur(T_schur_ord3,A_schur_ord3,select);
E5 = ordeig(A_schur_ord4);
E5_aug = [abs(E5) [1:18]'];
E5_aug_sort = sortrows(E5_aug,1);

select = zeros(1,18);
select(1:10) = 1;
select(15) = 1;

% Swap 6
[T_schur_ord5,A_schur_ord5] = ordschur(T_schur_ord4,A_schur_ord4,select);
E6 = ordeig(A_schur_ord5);
E6_aug = [abs(E6) [1:18]'];
E6_aug_sort = sortrows(E6_aug,1);

select = zeros(1,18);
select(1:11) = 1;
select(15) = 1;
select(18) = 1;

% Swap 7
[T_schur_ord6,A_schur_ord6] = ordschur(T_schur_ord5,A_schur_ord5,select);
E7 = ordeig(A_schur_ord6);
E7_aug = [abs(E7) [1:18]'];
E7_aug_sort = sortrows(E7_aug,1);

select = zeros(1,18);
select(1:12) = 1;
select(16) = 1;

% Swap 8
[T_schur_ord7,A_schur_ord7] = ordschur(T_schur_ord6,A_schur_ord6,select);
E8 = ordeig(A_schur_ord7);
E8_aug = [abs(E8) [1:18]'];
E8_aug_sort = sortrows(E8_aug,1);

select = zeros(1,18);
select(1:14) = 1;
select(16) = 1;
select(18) = 1;

% Swap 9
[T_schur_ord8,A_schur_ord8] = ordschur(T_schur_ord7,A_schur_ord7,select);

E9 = ordeig(A_schur_ord8);
E9_aug = [abs(E9) [1:18]'];
E9_aug_sort = sortrows(E9_aug,1);

T_ordered = T_schur_ord8;

A_schur = A_schur_ord8;
B_schur = T_ordered*B_aug;
C_schur = C/T_ordered; 
