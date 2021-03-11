%INITIAL CONDITIONS
dt=0.5;
dx=0.01;
x_inf=25;%25 for short term (50x50 lattice) / 50 for long-term (100x100 lattice)
x=0:dx:x_inf;
%SS
%Short-term
a=exp(-1*(x));
%Long-term
% a=exp(-0.5*(x));
%GOG
%Short-term
m=exp(-1*(x))*0.5;
p=exp(-1*(x))*0.5;
%Long-term
% m=exp(-0.5*(x))*0.5;
% p=exp(-0.5*(x))*0.5;


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
