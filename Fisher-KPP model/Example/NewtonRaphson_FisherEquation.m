clear all,close all,clc

%Newton-Raphson Method
C=[]; %to study the speed of the wave
for D=1:1
    
D=1;
k=1;
dx=1e-2;
dt=1e-2;
x_inf=50;
x=0:dx:x_inf;
epsilon=1e-8;
%initial condition
% u=zeros(length(x),1);
% 
% for i=1:length(x)
%     if x(i)<1
%         u(i)=0.5;
%     end
% end
%initial condition
u=0.5*heaviside(-x+1);
% u=1./(2.*sqrt(pi.*D.*(1/32))).*exp(-x.^2./(4.*D.*(1/32)));
u=u';
% u2=u;
plot(x,u,'Linewidth',2)
xlabel('x')
ylabel('u(x,0)')
% hold on
up=u;
% up2=u2;
T=14;
F=zeros(1,length(x));
% F2=zeros(1,length(x));
d1=zeros(1,length(x)-1);
d2=zeros(1,length(x));
d3=zeros(1,length(x)-1);
tic
for j=1:T/dt%time iterations
    t=dt*j
    

    m_u_up=0;
    du=0.0001*ones(length(x),1);
%     du2=0.0001*ones(length(x),1);
    while norm(du,inf)>epsilon %|| norm(du2,inf)>epsilon
                %---------------------------
       %Start Boundary conditions for i=1
        F(1)=u(2)-u(1);
%         F2(1)=u2(2)-u2(1);
        d1(1)=0;%dF/du_iminus1=0
        d2(1)=-1;%dF/du_i=-1
        d3(1)=1;%dF/du_iplus1=1
 
        g=0;
        for i=2:length(x)-1
           % function for each row of system F
            F(i)=-(u(i)-up(i))/dt+D*((u(i-1)-2*u(i)+u(i+1))/(dx^2))+k*u(i)*(1-u(i));

            %fill main diagonal of Jacobian (tridiagonal matrix)            
            d1(i)=D/(dx^2);
            d2(i)=-1/dt-2*(D/(dx^2))+k-2*k*up(i);
            d3(i)=D/(dx^2);
        end
        % End Boundary Conditions
        d1(end)=0;
        d2(end)=-1;
        d3(end)=1;
        F(end)=u(length(x))-u(length(x)-1);
%         F2(length(x))=u2(length(x))-u2(length(x)-1);
        %-------------------------------------

        %-----------------------------------------
        %construct tridiagonal matrix
        d1=[0 d1];
        d3=[d3 0];
        du=tridiag(d1',d2',d3',-F);
        %-------------------------------------------
        u=u+du;

    end   
    
     %keep record of previous u
     up=u;
%      up2=u;
%%%---------------------------------------------
%     plot(x,u); hold on
%     xlim([-10 10])
%     ylim([-0.5 2])
%     xlabel('x')
%     ylabel('u(x,t)')
    
        %ANIMATION -----------------------
%            plot(x,u)
% %            hold on
% %            plot(x,u2,'Linewidth',2)
% %            hold off
%            ylim([0 1])
% %            legend('du=-J\F','du=tridiag(d1,d2,d3,-F)')
%            pause(0.00001)
%         c=2*sqrt(D*k);
%         u_asymp=FisherKPP_AsympSol(x,dt*j,c);
%         plot(x,u,x,u_asymp)
%         legend('Numerical solution','Asymptotic solution')
%         xlabel('x')
%         ylabel('u(x,t)')
%         title('Fisher-KPP equation')
%         xlim([-50 50]) 
%         ylim([-0.05 2])
%         pause(0.1)
        

        %---------------------------------
%         %Plots for different times
        if t==0 || t==2||t==4||t==6||t==8||t==10||t==12||t==14
%             plot(x,u_asymp)
%             legend(strcat('t=',num2str(1*dt)),strcat('t=',num2str(50*dt)),strcat('t=',num2str(100*dt)))            
%             hold on
%             legend('Asymptotic solution','Numerical solution')
            plot(x,u,'Linewidth',2)
            hold on
            xlim([0 50])
            ylim([0 1])
            xlabel('x')
            ylabel('u(x,t)')
            
        end
%         
        
end
toc
end
