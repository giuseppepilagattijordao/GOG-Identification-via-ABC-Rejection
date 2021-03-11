%Phase portrait analysis to study point stability for Single-Species model
clear all 
clc
%Fisher KPP (example from Murray book
alpha=1;
nu=5;
D_alpha=alpha/2;
D_v=nu/2;
mu=0.001;
% c=sqrt(2*(alpha-mu)*(alpha+nu)); %for (0,0)
c=1
% c=sqrt(2*(mu+alpha)*(mu+nu));%for (1-mu/alpha,0)
[x, y, z] = meshgrid(-2:0.1:2);
u=y;
v=(-c*y-alpha*x.*(1-x)+mu*x)./(D_alpha*(1-x)+D_v);
w = z;      % z-component of vector field
h=streamslice(x, y, z, u, v, w, 2, 2, -1.5);
set( h, 'Color', [0 0 0] )
box on;
axis([-0.5 2 -1 2 -1.5 2]);
% Plot X and Y axis along y=0 and x=0.
line(xlim(), [0,0], 'Linewidth',1,'Color', 'k');
hold on
plot([0 0], ylim,'black','Linewidth',1)
hold on
%plot equilibrium points
scatter([0,1-mu/alpha],[0,0],200,'filled')
% grid on;
% view(-30, 60);
xlabel('u_1')
ylabel('u_2')
zlabel('z')
pause(0.1)



