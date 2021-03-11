%function to calculate analytical travelling wave speed
%following suggested Gerlee's numerical method
function [c] = AnalyticalSpeed(alpha,v,mu,q_p,q_m)
    
    %initialize c
    c=0;
%Jacobian matrix
    J=[0 0 1 0; 0 0 0 1;2/alpha*(q_m+mu-alpha) -2*q_p/alpha -2*c/alpha 0;-2*q_m/v 2*(q_p+mu)/v 0 -2*c/v];
%slight increase for c
    increase=0.0001;
    
    while isreal(eig(J))==0%while there is at least one complex eigenvalue
        J=[0 0 1 0; 0 0 0 1;2/alpha*(q_m+mu-alpha) -2*q_p/alpha -2*c/alpha 0;-2*q_m/v 2*(q_p+mu)/v 0 -2*c/v];
        c=c+increase;%change c slightly
    end
end