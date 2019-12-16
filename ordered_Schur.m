function [T_schur_ord,ordSys] = ordered_Schur(augSys)

sysSize = length(augSys.A);

[T_schur,A_schur_ord] = schur(augSys.A);

E = ordeig(A_schur_ord);
E_aug = [abs(E) [1:sysSize]'];
E_aug_sort = sortrows(E_aug,1);
