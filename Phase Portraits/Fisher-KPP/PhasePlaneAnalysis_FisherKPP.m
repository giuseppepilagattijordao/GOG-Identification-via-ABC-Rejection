%Phase portrait analysis to study point stability for FisherKPP model
clear all 
clc
%Fisher KPP (example from Murray book)
c=1;
[x, y, z] = meshgrid(-2:0.1:2);
%system derived from Fisher-KPP equation:
u=y;
v=-c*y-x.*(1-x); 
w = z;      % z-component of vector field
h=streamslice(x, y, z, u, v, w, 2, 2, -1.5);%create streamlines
set( h, 'Color', [0 0 0] )
box on;
axis([-0.5 2 -1 2 -1.5 2]);
% Plot X and Y axis along y=0 and x=0.
line(xlim(), [0,0], 'Linewidth',1,'Color', 'k');
hold on
plot([0 0], ylim,'black','Linewidth',1)
hold on
%plot equilibrium points
scatter([0,1],[0,0],200,'filled')
% grid on;
% view(-30, 60);
xlabel('u_1')
ylabel('u_2')
zlabel('z')
pause(0.1)



