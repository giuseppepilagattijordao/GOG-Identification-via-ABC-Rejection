function [U]=FisherKPP_AsympSol(x,t,c)
%Asumptotic solution for Fisher KPP equation
%for epsilon = 1/c^2 <= 0.25,
z=(x-c*t);
% z=x;
epsilon=1/c^2;

U=(1+exp(z/c)).^(-1)+((1/c^2)*exp(z/c).*(1+exp(z/c)).^(-2))*log((4*exp(z/c))/((1+exp(z/c)).^2));

end