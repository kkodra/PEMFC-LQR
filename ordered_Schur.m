function [T_schur_ord,ordSys] = ordered_Schur(augSys)
% Function evaluates the ordered Schur realization of the input state space model

sysSize = length(augSys.A);

[T_schur,A_schur_ord] = schur(augSys.A);

E = ordeig(A_schur_ord);
E_aug = [abs(E) [1:sysSize]'];
E_aug_sort = sortrows(E_aug,1);

while ~issorted(diag(A_schur_ord),'descend')
    sel = E_aug_sort(:,2);

    count = 1;
    sel_vec = [];
    while (sel(count+1) > sel(count))
        sel_vec(count) = sel(count);
        count = count+1;
    end
    sel_vec(length(sel_vec)+1) = sel(count);
    
    select = zeros(1,sysSize);
    select(sel_vec) = 1;


    [T_schur,A_schur_ord] = ordschur(T_schur,A_schur_ord,select);
    
    E = ordeig(A_schur_ord);
    E_aug = [abs(E) [1:sysSize]'];
    E_aug_sort = sortrows(E_aug,1);

end

T_schur_ord = T_schur;
B_schur_ord = T_schur_ord*augSys.B;
C_schur_ord = augSys.C/T_schur_ord;

ordSys.A = A_schur_ord;
ordSys.B = B_schur_ord;
ordSys.C = C_schur_ord;

