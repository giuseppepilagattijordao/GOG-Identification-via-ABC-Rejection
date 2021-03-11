
%ABC rejection for comparing Continuum Approximations
clear all
clc
tic
winner=[];%register winner for all trials
LR=[];
AIC=[];
BIC=[];
%-----------------------------------------------------------
%GO OR GROW
%lists to store accepted values of each parameter
alpha_vector=[];
v_vector=[];
mu_vector=[];

v=1; 
alpha=0.01;
mu=0.001;
% beta=10;
q_p=1;
q_m=1;
%INITIAL CONDITIONS
dt=0.5;
dx=0.01;
x_inf=50;%25 for short term (50x50 lattice) / 50 for long-term (100x100 lattice)
x=0:dx:x_inf;
%SS
% a=ones(1,length(x));
% for i=1:length(a)
% %     if i<=find(x==x_inf/2)-10/dx || i>=find(x==x_inf/2)+10/dx
%       if i>=find(x==10)
%         a(i)=0;
%     end
% end
%Short-term
% a=exp(-1*(x));
%Long-term
a=exp(-0.5*(x));
%GOG
% m=0.5*ones(length(x),1);
% p=0.5*ones(length(x),1);
% for i=1:length(x)
% %     if i<=find(x==x_inf/2)-10/dx || i>=find(x==x_inf/2)+10/dx
%     if i>=find(x==10)
%         p(i)=0;
%         m(i)=0;
%     end
% end
%Short-term
% m=exp(-1*(x))*0.5;
% p=exp(-1*(x))*0.5;
%Long-term
m=exp(-0.5*(x))*0.5;
p=exp(-0.5*(x))*0.5;


figure
plot(x,a,'black','LineWidth',2)
ylim([0 1])
xlabel('x')
ylabel('a(x, t = 0)')
title('Single-species')
figure
plot(x,m,'blue','LineWidth',2)
ylim([0 1])
xlabel('x')
ylabel('m(x, t = 0)')
title('GOG: Motile cells')
figure
plot(x,p,'red','LineWidth',2)
ylim([0 1])
xlabel('x')
ylabel('p(x, t = 0)')
title('GOG: Proliferative cells')
% 
% toc