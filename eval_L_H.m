function [L,H] = eval_L_H(A1,A2,A3,A4,epsilon,method)
% Function evaluates matrices L and H of the Chang transformation using:
%           - Recursive method: 
%              Ref: Gajic and Grodt, "The recursive reduced-order numerical solution of the singularly perturbed matrix differential equation," IEEE TAC, 1988.
%           - Eigenvectors method:
%              Ref: Kecman, Bingulac, and Gajic, "Eigenvector approach for order-reduction of singularly perturbed LQ optimal control problem," Automatica, 1999.
% Syntax: [L,H] = eval_L_H(A1,A2,A3,A4,epsilon,method)
% User can select 'recursive' or 'eigenvector' for method.

switch method
case 'recursive'
    Li = (A4)\A3;
    for j = 1:1000
        D1 = A4 +epsilon*Li*A2;
        D2 = -epsilon*(A1 - A2*Li);
        Qi = A3 +epsilon*Li*A2*Li;
        Li_1 = lyap(D1,D2,-Qi);
        Hi_1 = lyap(D2,D1,-A2);
        Lnew = Li_1;
        Hnew = Hi_1;
    end
    L = Lnew;
    H = Hnew;
case 'eigenvector'
    R = [-epsilon*A1 epsilon*A2; A3 -A4];
    [V,~] = eig(R);
otherwise
    error('Method options are \'recursive\' and \'eigenvector\'.');
    M = V;
    i = 1;
    Mnew = zeros(length(M));
    while i  <= length(R)
        if isreal(M(:,i))== 0
            Mnew(:,i) = real(M(:,i));
            Mnew(:,i+1) = imag(M(:,i));
            if i == length(R)
                break
            end
            i = (i+1) + 1;
        else
            Mnew(:,i) = M(:,i);
            i = i + 1;
        end
    end

    M1 = Mnew(:,1:length(A1));
    M11 = M1(1:length(A1),:);
    M21 = M1(length(A1)+ 1:end,:);
    L = M21/M11;

    % Find L
    % Test_L = A4*L - epsilon*L*A1 - A3 + L*epsilon*A2*L   % Test is good

    Asyl = -epsilon*(A1-A2*L);
    Bsyl = A4 + epsilon*L*A2;
    Csyl = -A2;

    H = lyapkr(Asyl,Bsyl,Csyl);
end
end

