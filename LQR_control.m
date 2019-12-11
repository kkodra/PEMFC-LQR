function [K, sys_FB, y, x] = LQR_control(Q, R, sys, r, t)
% LQR controller design function

if nargin < 5
	error('Number of arguments is less than accepted.');
end

A = sys.A;
B = sys.B;
C = sys.C;

K = lqr(A, B, Q, R);

A_fb = A - B*K;
B_fb = B;
C_fb = C;
D_fb = 0;

sys_FB.A = A_fb;
sys_FB.B = B_fb;
sys_FB.C = C_fb;
sys_FB.D = D_fb;

sys_ss = ss(A_fb, B_fb, C_fb, D_fb);

[y,~,x] = lsim(sys_ss, r, t);



