function[Nbar]=rscale(model_ss,K)
         % Function rscale(sys,K) 
         % finds the scale factor N which will
         % eliminate the steady-state error to a step reference
         % for a continuous-time, single-input system
         % with full-state feedback.
         % Based on implementation by J. Yook and D. Tilbury (University of Michigan)
         [A,B,C,D] = ssdata(model_ss);
         
         % compute Nbar
         s = size(A,1);
         Z = [zeros([1,s]) 1];
         N = inv([A,B;C,D])*Z';
         Nx = N(1:s);
         Nu = N(1+s);
         Nbar=Nu + K*Nx;
