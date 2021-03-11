%GERLEE CONTINUUM MODEL (EXPLICIT)

clear all,close all,clc

%NON DIMNESIONALIZE (TIME AND SPACE)
C=[]; %for all wave speeds
C_analytical=[];
Xf=[];
Xi=[];
%PARAMETERS
v=5;
alpha=1;
mu=1e-3;
beta=10;
q_p=10;
q_m=0;
%CALCULATE ANALYTICAL TRAVELLING WAVE SPEEDS
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];

%NUMERICAL PARAMETERS
h=1;
D_alpha=(h^2)*alpha/2;
D_v=(h^2)*v/2;
dx=0.1;
dt=(dx^2/2)/100;
% dt=5e-5;
r=dt/dx^2;
% dt=D_alpha^(-1)*1e-1;
x_inf=100;
x=0:dx:x_inf;
CellCycles=50;
epsilon=1e-8;
%initial conditions
p=exp(-beta*x);
m=zeros(1,length(x));
% m(:,1)=zeros(length(x),1);
p_next=p;
m_next=m;
    for j=1:CellCycles/dt
        [j r CellCycles/dt]
        for i=2:length(x)-1
            %Crank-Nicolson (Average between implicit and explicit
            d2pdx2=(1/2)*((p_next(i-1)+p_next(i+1)-2*p_next(i))/dx^2+(p(i-1)+p(i+1)-2*p(i))/dx^2);
            d2mdx2=(1/2)*((m_next(i-1)+m_next(i+1)-2*m_next(i))/dx^2+(m(i-1)+m(i+1)-2*m(i))/dx^2);
%             d2pdx2=(p(i-1)+p(i+1)-2*p(i))/dx^2;
%             d2mdx2=(m(i-1)+m(i+1)-2*m(i))/dx^2;
            %explicit algebraic equation
            p_next(i)=p(i)+dt*((D_alpha)*(1-p(i)-m(i))*(d2pdx2)+alpha*p(i)*(1-p(i)-m(i))-(q_m+mu)*p(i)+q_p*m(i));
            m_next(i)=m(i)+dt*((D_v)*((1-p(i))*(d2mdx2)+m(i)*(d2pdx2))-(q_p+mu)*m(i)+q_m*p(i));
        end
        %Bondary conditions (no flux)
        p_next(1)=p_next(2);
        m_next(1)=m_next(2);
        p_next(length(x))=p_next(length(x)-1);
        m_next(length(x))=m_next(length(x)-1);
        %animation
%         plot(x,p,x,m,'--','LineWidth',2)
%         ylim([0 1])
%         xlim([0 100])
%         pause(0.0000001)
        %update
        p=p_next;
        m=m_next;
end
%     
plot(x,p,'black','LineWidth',2)
hold on
plot(x,m,'--black','LineWidth',2)
legend('Proliferative','motile')
        ylim([0 1])
        xlim([0 100])
        xlabel('x')
ylabel('Probability')