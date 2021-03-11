clear all,close all,clc

%Newton-Raphson Method
C=[]; %to study the speed of the wave
% for D=1:1
%     D
D=6;
k=1;
dx=1e-1;
dt=1e-1;
x_inf=50;
x=-50:dx:x_inf;
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
% u=0.5*heaviside(-x+1);
u=1./(2.*sqrt(pi.*D.*(1/32))).*exp(-(x).^2./(4.*D.*(1/32)));
u=u';
% u2=u;
% plot(x,u,'Linewidth',2)
% xlabel('x')
% ylabel('u(x,0)')
% hold on
up=u;
% up2=u2;
T=10;
F=zeros(1,length(x));
% F2=zeros(1,length(x));
d1=zeros(1,length(x)-1);
d2=zeros(1,length(x));
d3=zeros(1,length(x)-1);
tic
for j=1:T/dt%time iterations
    t=dt*j
  
%     [t D]
%     if t==4 %evaluate wave's speed using two snap shots
%         ui=u;
%         ui_point=max(find(abs(ui-0.5)<=0.01));%using right-most wave
%     elseif t==7
%         uf=u;
%         uf_point=max(find(abs(uf-0.5)<=0.01));%using right-most wave
%         c=dx*abs(uf_point-ui_point)/(7-4);%multiply by dx to take into account discretization of rod
%         C=[C c];
%     end
    m_u_up=0;
    du=0.0001*ones(length(x),1);
%     du2=0.0001*ones(length(x),1);
    while norm(du,inf)>epsilon %|| norm(du2,inf)>epsilon
%         norm(du,Inf)
%         norm(du2,inf)
%         m_u_up=m_u_up+norm(abs(u-up));
                %---------------------------
       %Start Boundary conditions for i=1
        F(1)=u(2)-u(1);
%         F2(1)=u2(2)-u2(1);
        d1(1)=0;%dF/du_iminus1=0
        d2(1)=-1;%dF/du_i=-1
        d3(1)=1;%dF/du_iplus1=1
 
%         if size(d1,2)>1
%             d1=d1';
%             d2=d2';
%             d3=d3';
%         end
        g=0;
        for i=2:length(x)-1
            %define algebraic nonlinear equation to solve
            %f=@(u1,u2,u3)(-u2+up(i))/dt+(D/(dx^2))*(u1-2*u2+u3)+k*u2*(1-u2);
           % function for each row of system F
            %F(i)=f(u(i-1),u(i),u(i+1));%system F is a vector of each function
            F(i)=-(u(i)-up(i))/dt+D*((u(i-1)-2*u(i)+u(i+1))/(dx^2))+k*u(i)*(1-u(i));
%             F(i)=-(up(i)-u(i))/dt+D*((up(i-1)-2*up(i)+up(i+1))/(dx^2))+k*up(i)*(1-up(i));
%             F2(i)=-(u2(i)-up2(i))/dt+D*((u2(i-1)-2*u2(i)+u2(i+1))/(dx^2))+k*u2(i)*(1-u2(i));

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
%         J=full(gallery('tridiag',length(F2),D/(dx^2),-1/dt-2*(D/(dx^2))+k-2*k*u2(i),D/(dx^2)));
%         J(1,1)=-1;
%         J(2,1)=0;
%         J(1,2)=1;
%         J(end,end)=-1;
%         J(end,end-1)=0;
%         J(end-1,end)=1;
%         du2=-J\F2';%solve linear system
        d1=[0 d1];
        d3=[d3 0];
        du=tridiag(d1',d2',d3',-F);
        %-------------------------------------------
        u=u+du;
%         u2=u2+du2;
       
%         plot(x,du,x,u)
%         xlim([-1.5 1.5])
%         ylim([-0.05 0.7])
%         legend('du','u')
%         pause(0.01)
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
%            hold on
%            plot(x,u2,'Linewidth',2)
%            hold off
%            ylim([0 1])
%            legend('du=-J\F','du=tridiag(d1,d2,d3,-F)')
%            pause(0.00001)
        c=2*sqrt(D*k);
        time_lag=0; %time delay of numerical solution to reach travelling wave shape from initial condition
        space_lag=11; %spacial delay of numerical travelling wave due to time delay
        u_asymp=FisherKPP_AsympSol(x+space_lag,t-time_lag,c);
%         plot(x,u,x,u_asymp)
%         legend('Numerical solution','Asymptotic solution')
%         xlabel('x')
%         ylabel('u(x,t)')
%         title('Fisher-KPP equation')
%         xlim([-50 50]) 
%         ylim([-0.05 2])
%         pause(0.001)
        

        %---------------------------------
%         %Plots for different times
            if t==5 || t==7.5 || t==10
                plot(x,u,'Linewidth',2)
                hold on
                plot(x,u_asymp)
                xlim([-50 50])
                ylim([0 1.5])
                xlabel('x')
                ylabel('u(x,t)')
                legend('t=5','t=7.5','t=10')
            end

%         if t==0 || t==2||t==4||t==6||t==8||t==10||t==12||t==14
% %             plot(x,u_asymp)
% %             legend(strcat('t=',num2str(1*dt)),strcat('t=',num2str(50*dt)),strcat('t=',num2str(100*dt)))            
% %             hold on
% %             legend('Asymptotic solution','Numerical solution')
%             plot(x,u,'Linewidth',2)
%             hold on
%             xlim([0 50])
%             ylim([0 1])
%             xlabel('x')
%             ylabel('u(x,t)')
%             
%         end
%         
        
        %Travelling wave ansazt
%         c=2;
%         [t,y]=ode45(@TravWaveAnsazt,[0 dt],[u(1);F(1)]);
%         if mod(length(y(:,1)),2)==0 %if even
%             xansazt=-50:(100)/length(y(:,1)):50;
%         else
%             xansazt=-50:(100)/length(y(:,1)):50-1;
%         end
%         plot(xansazt',y(:,1));hold on %get only solution (y1=solution and y2=derivative)
%         pause(0.01)
%         xlabel('x')
%         ylabel('u(x,t)')
        
        %FOURIER TRANSFORM
%         A=FKPP_Fourier(u,x_inf,1e-8,D,x);
%         norm(abs(A-u))
%         plot(x,A)
%         legend('Newton','Fourier')
%         pause(0.1)
%         hold off
%         hold off
%         
end
toc
% end
%Wave Speed plot
% plot(1:20,C)


% figure
%  plot(x,u)
%     xlim([-50 50])
%     ylim([-0.5 2])
%     xlabel('x')
%     ylabel('u(x,t)')
%     legend(strcat('t=',num2str(j)))
    
