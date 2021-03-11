%GERLEE CONTINUUM MODEL

clear all%close all,clc

%NON DIMNESIONALIZE (TIME AND SPACE)
C=[]; %for all wave speeds
C_analytical=[];
Xf=[];
Xi=[];
DX=[0.1 0.01 0.001];
Dt=[0.1];
figure
% for dx=DX
% for dx=DX
for q_m=0:30

%PARAMETERS
v=5;
alpha=1;
mu=1e-3;
beta=10;
q_p=22;
% q_m=10;

%CALCULATE ANALYTICAL TRAVELLING WAVE SPEEDS
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];

%NUMERICAL PARAMETERS
h=1;
D_alpha=(h^2)*alpha/2;
D_v=(h^2)*v/2;
dx=0.01;
dt=0.01;
% dt=0.001;%dt=0.01 means that for T=50 we need 50/0.01=5000 iterations

% dt=D_alpha^(-1)*1e-1;
x_inf=100;
x=0:dx:x_inf;

epsilon=1e-8;
%initial condition
% m=0.5*ones(length(x),1);
% p=0.5*ones(length(x),1);
% for i=1:length(x)
%     if i<=find(x==200)-20/dx || i>=find(x==200)+20/dx
%         p(i)=0;
%         m(i)=0;
%     end
% end
p=exp(-beta*x);
p=p';
m=zeros(length(x),1);
% plot(x,p);hold on
% plot(x,m)
% xlabel('x')
% ylabel('p(x,0)')
p_previous=p;
m_previous=m;
%     SS1=[];
%     SS2=[];
%     SS3=[];
%     SS4=[];
%     SS5=[];
%     SS6=[];
P=zeros(1,length(x));
M=zeros(1,length(x));
d1_p=zeros(1,length(x)-1);
d2_p=zeros(1,length(x));
d3_p=zeros(1,length(x)-1);
d1_m=zeros(1,length(x)-1);
d2_m=zeros(1,length(x));
d3_m=zeros(1,length(x)-1);
CellCycles=56;
j50_check=1;%for the snapshots to be at the right time point
j56_check=1;
for j=0:dt:CellCycles %time iterations
    t=dt*j;
    j
%     clear d1_p d2_p d3_p d1_m d2_m d3_m P M %clear diagonal values at each cycle
    if j==50 %evaluate wave's speed using two snap shots
        j50_check=abs(j-50);
        p_i=p;
        m_i=m;
        pm_i=p_i+m_i; %sum of m and p
    elseif j==56
        j56_check=j;
        p_f=p;
        m_f=m;
        pm_f=p_f+m_f; %sum of m and p
        %get point close to 1/2
        x_i=min(find(abs(pm_i-max(pm_i)/2)<=(0.005)));
        x_f=min(find(abs(pm_f-max(pm_f)/2)<=(0.005)));
        c=(dx*x_f-dx*x_i)/(56-50);%wave speed
        Xf=[Xf x_f];
        Xi=[Xi x_i];
        C=[C c];
    end

%     j
    m_u_up=0;
    dp=0.0001*ones(length(x),1);
    dm=0.0001*ones(length(x),1);
    %vectors to check derivatives
%     S1=[];
%     S2=[];
%     S3=[];
%     S4=[];
%     S5=[];
%     S6=[];

    while norm(dp,inf)>epsilon || norm(dm,inf)>epsilon
%         [norm(dp,inf),norm(dm,inf)]
       %---------------------------
       %Start Boundary conditions for i=1
        P(1)=p(2)-p(1);
        M(1)=m(2)-m(1);
        d1_p(1)=0;%dF/du_iminus1=0
        d2_p(1)=-1;%dF/du_i=-1
        d3_p(1)=1;%dF/du_iplus1=1
        d1_m(1)=0;%dF/du_iminus1=0
        d2_m(1)=-1;%dF/du_i=-1
        d3_m(1)=1;%dF/du_iplus1=1

        for i=2:length(x)-1
            
            %define algebraic nonlinear equation to solve
            %f=@(u1,u2,u3)(-u2+up(i))/dt+(D/(dx^2))*(u1-2*u2+u3)+k*u2*(1-u2);
           % function for each row of system F
            %F(i)=f(u(i-1),u(i),u(i+1));%system F is a vector of each function
%             P(i)=-(p(i)-p_previous(i))/dt+(D_alpha)*(1-p(i)-m(i))*(p(i-1)-2*p(i)+p(i+1))/(dx^2)+alpha*p(i)*(1-p(i)-m(i))-(q_m+mu)*p(i)+q_p*m(i);
%             M(i)=-(m(i)-m_previous(i))/dt+(D_v)*((1-p(i))*(m(i-1)-2*m(i)+m(i+1))/(dx^2)+m(i)*(p(i-1)-2*p(i)+p(i+1))/(dx^2))-(q_p+mu)*m(i)+q_m*p(i);
            
            %Crank-Nicolson (Average between implicit and explicit
%             d2pdx2=(1/2)*((p_previous(i-1)+p_previous(i+1)-2*p_previous(i))/dx^2+(p(i-1)+p(i+1)-2*p(i))/dx^2);
%             d2mdx2=(1/2)*((m_previous(i-1)+m_previous(i+1)-2*m_previous(i))/dx^2+(m(i-1)+m(i+1)-2*m(i))/dx^2);
%             IMPLICIT
            d2pdx2=(p(i-1)+p(i+1)-2*p(i))/dx^2;
            d2mdx2=(m(i-1)+m(i+1)-2*m(i))/dx^2;
            %Algebraic equation
            P(i)=-(p(i)-p_previous(i))/dt+(D_alpha)*(1-p(i)-m(i))*(d2pdx2)+alpha*p(i)*(1-p(i)-m(i))-(q_m+mu)*p(i)+q_p*m(i);
            M(i)=-(m(i)-m_previous(i))/dt+(D_v)*((1-p(i))*(d2mdx2)+m(i)*(d2pdx2))-(q_p+mu)*m(i)+q_m*p(i);
            if isnan(P(i))||isnan(M(i))
                display('stop')
            end
            %numerical derivatives
%             diff_piminus1=@(p_i,m_i)1/2-p_i/2-m_i/2;
%             diff_pi=@(p_iplus1,p_iminus1)-p_iplus1/2-p_iminus1/2-30001/1000;
%             diff_piplus1=@(p_i,m_i)1/2-p_i/2-m_i/2;
%             diff_miminus1=@(p_i)5/2-(5*p_i)/2;
%             diff_mi=@(p_iplus1,p_iminus1)(5*p_iplus1)/2+(5*p_iminus1)/2-25001/1000;
%             diff_miplus1=@(p_i)5/2-(5*p_i)/2;
            %fill main diagonal of Jacobian (tridiagonal matrix)            
%             d1_p(i)=(D_alpha/dx^2)*(1-p(i)-m(i));
%             d2_p(i)=-1/dt+(D_alpha/dx^2)*(-p(i-1)+4*p(i)-p(i+1)-2*m(i)-2)+alpha*(1-2*p(i)-m(i))-(q_m+mu);
%             d3_p(i)=(D_alpha/dx^2)*(1-p(i)-m(i));
%             d1_m(i)=(D_v/dx^2)*(1-p(i));
%             d2_m(i)=-1/dt+(D_v/dx^2)*(-2+p(i-1)+p(i+1))-(q_p+mu);
%             d3_m(i)=(D_v/dx^2)*(1-p(i));
            
            d1_p(i)=(D_alpha/(dx^2))*(1-p(i)-m(i));
            d2_p(i)=-1/dt+(D_alpha/(dx^2))*(-p(i-1)+4*p(i)-p(i+1)+2*m(i)-2)+alpha*(1-2*p(i)-m(i))-(q_m+mu);
            d3_p(i)=(D_alpha/(dx^2))*(1-p(i)-m(i));
            d1_m(i)=(D_v/(dx^2))*(1-p(i));
            d2_m(i)=-1/dt+(D_v/(dx^2))*(-2+p(i-1)+p(i+1))-(q_p+mu);
            d3_m(i)=(D_v/(dx^2))*(1-p(i));
            
% %             
            %TEST------------
%             d1_p(i)=diff_piminus1(p(i),m(i));
%             d2_p(i)=diff_pi(p(i+1),p(i-1));
%             d3_p(i)=diff_piplus1(p(i),m(i));
%             d1_m(i)=diff_miminus1(p(i));
%             d2_m(i)=diff_mi(p(i+1),p(i-1));
%             d3_m(i)=diff_miplus1(p(i));
            %END OF TEST-----
            %comparing analytical derivatives and numerical ones
%             piminus1_comparison=abs(diff_piminus1(p(i),m(i))-d1_p(i));
%             pi_comparison=abs(diff_pi(p(i+1),p(i-1))-d2_p(i));
%             piplus1_comparison=abs(diff_piplus1(p(i),m(i))-d3_p(i));
%             miminus1_comparison=abs(diff_miminus1(p(i))-d1_m(i));
%             mi_comparison=abs(diff_mi(p(i+1),p(i-1))-d2_m(i));
%             miplus1_comparison=abs(diff_miplus1(p(i))-d3_m(i));
%             S1=[S1 piminus1_comparison];
%             S2=[S2 pi_comparison];
%             S3=[S3 piplus1_comparison];
%             S4=[S4 miminus1_comparison];
%             S5=[S5 mi_comparison];
%             S6=[S6 miplus1_comparison];
           
        end
        % End Boundary Conditions

        d1_p(end)=0;
        d2_p(end)=-1;
        d3_p(end)=1;
        d1_m(end)=0;
        d2_m(end)=-1;
        d3_m(end)=1;
        P(end)=p(end)-p(end-1);
        M(end)=m(end)-m(end-1);
%         length(x)
%         length(d1_p)
%         length(P)
        %-------------------------------------
%           solve Jdu=-F
%         Jp=zeros(length(x),length(x));
%         Jm=zeros(length(x),length(x));
%         Jp(1,1)=-1;
%         Jp(2,1)=0;
%         Jp(1,2)=1;
%         Jm(1,1)=-1;
%         Jm(2,1)=0;
%         Jm(1,2)=1;
%         for o=2:length(x)-1
%           
%           Jp(o+1,o)=(D_alpha/(dx^2))*(1-p(i)-m(i));
%           Jp(o,o)=-1/dt+(D_alpha/(dx^2))*(-p(i-1)+4*p(i)-p(i+1)+2*m(i)-2)+alpha*(1-2*p(i)-m(i))-(q_m+mu);
%           Jp(o,o+1)=(D_alpha/(dx^2))*(1-p(i)-m(i));
%           Jm(o+1,o)=(D_v/(dx^2))*(1-p(i));
%           Jm(o,o)=-1/dt+(D_v/(dx^2))*(-2+p(i-1)+p(i+1))-(q_p+mu);
%           Jm(o,o+1)=(D_v/(dx^2))*(1-p(i));
%           
%         end
%           Jp(length(x),length(x))=-1;
%           Jp(length(x)-1,length(x))=1;
%           Jp(length(x),length(x)-1)=0;
%           Jm(length(x),length(x))=-1;
%           Jm(length(x)-1,length(x))=1;
%           Jm(length(x),length(x)-1)=0;
%         
%          dp=-Jp\P';%solve linear system
%          dm=-Jm\M';
        %-----------------------------------------
        d1_p=[0 d1_p];
        d3_p=[d3_p 0];
        d1_m=[0 d1_m];
        d3_m=[d3_m 0];
        dp=tridiag(d1_p',d2_p',d3_p',-P);
        dm=tridiag(d1_m',d2_m',d3_m',-M);
%         dp=Thomas(d1_p',d2_p',d3_p',-P);
%         dm=Thomas(d1_m',d2_m',d3_m',-M);
        if isnan(dp(1))||isnan(dm(1))
            display('stop')
        end
        %-------------------------------------------
        p=p+dp;
        m=m+dm;
        
        
    end   
      %keep record of previous u
     p_previous=p;
     m_previous=m;
%%%---------------------------------------------
%     
%         plot(x,p,x,m,'--','LineWidth',2)
%         ylim([0 1])
%         pause(0.0000001)
%         xlabel('x')
%         ylabel('probability')
%         legend('proliferating','moving')
% %         hold off
%         xlim([0 100])
%         ylim([-0.01 1])
%         pause(0.01)
% %     if j==4000
% %         figure()
%         plt=plot(x,p,x,m,'--','LineWidth',2)
%         plt(1).LineWidth=1.5;
%         plt(1).Color='black';
%         plt(2).Color='black';
%         hold on
%         xlabel('x')
%         ylabel('probability')
%         legend(plt,'proliferating','moving')
% %         hold off
%         xlim([0 100])
%         ylim([-0.01 1])
% %         pause(0.001)
%     elseif j==5000
%         plt=plot(x,p,x,m,'--','LineWidth',2)
%         plt(1).LineWidth=1.5;
%         plt(1).Color='red';
%         plt(2).Color='red';
% %         hold on
%         xlabel('x')
%         ylabel('probability')
%         legend(plt,'proliferating','moving')
% %         hold off
%         xlim([0 100])
%         ylim([-0.01 1])
%     end
%     hold off
%     hold off
%         SS1=[SS1 mean(S1)];
%         SS2=[SS2 mean(S2)];
%         SS3=[SS3 mean(S3)];
%         SS4=[SS4 mean(S4)];
%         SS5=[SS5 mean(S5)];
%         SS6=[SS6 mean(S6)];
end

% if dt==DeltaT(1)
%     plot(x,p_i,'b')
%     hold on
% elseif dt==DeltaT(2)
%     plot(x,p_i,'r')
%     hold on
% elseif dt==DeltaT(3)
%     plot(x,p_i,'g')
%     hold on
% elseif dt==DeltaT(4)
%     plot(x,p_i,'black')
%     hold 
% plot(x,p,'blue','Linewidth',2)
% hold on
% plot(x,m,'blue--','LineWidth',2)
% ylim([0 1])
% plot(x,p_i,'green','Linewidth',2)
% xlabel('x')
% ylabel('Probability')
% hold on
% plot(x,m_i,'green--','LineWidth',2)
% legend(num2str(beta))
end
%Plot
plot(C)
xlabel('q_m (cell cycle^{-1})')
ylabel('wave speed c')

% legend(['dx=',num2str(dx),';dt=',num2str(DeltaT(1))],['dx=',num2str(dx),';dt=',num2str(DeltaT(2))],['dx=',num2str(dx),';dt=',num2str(DeltaT(3))],['dx=',num2str(dx),';dt=',num2str(DeltaT(4))],['dx=',num2str(dx),';dt=',num2str(DeltaT(5))])
% 

 
%  plot(x,m_i,'Linewidth',2);hold on

%EVALUATE DIFFERENT DT
% if dt==Dt(1)
%     p_f_dt1=p_f;
%     m_f_dt1=m_f;
% elseif dt==Dt(2)
%     p_f_dt2=p_f;
%     m_f_dt2=m_f;
% elseif dt==Dt(3)
%     p_f_dt3=p_f;
%     m_f_dt3=m_f;
% elseif dt==Dt(4)
%     p_f_dt4=p_f;
%     m_f_dt4=m_f;
% elseif dt==Dt(5)
%     p_f_dt5=p_f;
%     m_f_dt5=m_f;
% end

%EVALUATE DIFFERENT DX
% if dx==DX(1)
%     p_f_dx1=p_f;
%     m_f_dx1=m_f;
% %   plot(x,p_f,'Linewidth',2);hold on
% elseif dx==DX(2)
%     p_f_dx2=p_f;    
%     m_f_dx2=m_f;
% elseif dx==DX(3)
%     p_f_dx3=p_f;    
%     m_f_dx3=m_f;
% end
%     clear error_independence
%     error_independence=zeros(1,length(p_f_dx1)); %length of the shortest one
%  plot(x,p_f,'Linewidth',2)
%     xlim([0 100])
%     ylim([0 1])
%     xlabel('x')
%     ylabel('u(x,t)')
%     hold off
%     if dt==DeltaT(1)
%         for i=1:length(p_f_dx1)%get vector of common points of different grids for comparison
%             error_independence(i)=abs(p_f_dx1(i)-p_f_dx2(2*(i-1)+1));
%         end
%     elseif dt==DeltaT(2)
%         for i=1:length(p_f_dx1)%get vector of common points of different grids for comparison
%             error_independence(i)=abs(p_f_dx1(i)-p_f_dx2(2*(i-1)+1));
%         end
%     end
%     plot(error_independence)
%     hold on
%     legend(['dt=',num2str(DeltaT(1))],['dt=',num2str(DeltaT(2))])
%     title('absolute error in common points')
% elseif dx==DX(3)
%     p_f_dx3=p_f;
%     m_f_dx3=m_f;
% end

%     legend(strcat('t=',num2str(j)))

%c=1.48, 1.88 and 1.63
%(q_p,q_m)=(10,0),(20,10) and (10,20)
%WAVE SPEED
% c=(cp_f-cp_i)/(500-400);
% C=[C c]; %all wave speeds for parameter analysis
% end
% plot(0:100/(length(p_f_dt1)-1):100,p_f_dt1)
% hold on
% plot(0:100/(length(p_f_dt2)-1):100,p_f_dt2)
% hold on
% plot(0:100/(length(p_f_dt3)-1):100,p_f_dt3)
% hold on
% plot(0:100/(length(m_f_dt1)-1):100,m_f_dt1,'--')
% hold on
% plot(0:100/(length(m_f_dt2)-1):100,m_f_dt2,'--')
% hold on
% plot(0:100/(length(m_f_dt3)-1):100,m_f_dt3,'--')
% ylim([0 1])
% legend('dt=0.005','dt=0.001','dt=0.0005')
% plot(0:100/(length(p_f_dx1)-1):100,p_f_dx1)
% hold on
% plot(0:100/(length(p_f_dx2)-1):100,p_f_dx2)
% hold on
% plot(0:100/(length(p_f_dx3)-1):100,p_f_dx3)
% hold on
% plot(0:100/(length(m_f_dx1)-1):100,m_f_dx1,'--')
% hold on
% plot(0:100/(length(m_f_dx2)-1):100,m_f_dx2,'--')
% hold on
% plot(0:100/(length(m_f_dx3)-1):100,m_f_dx3,'--')

%PLOT SUM OF P AND M
% plot(x,pm_f,'Linewidth',2)
% xlabel('x')
% ylabel('Probability')
%---------------PLOT P OR M---------------------------
% plot(x,p+m,'Linewidth',2)
% xlabel('x')
% ylabel('Probability')
% plot(x,p,'blue','Linewidth',2)
% hold on
% plot(x,m,'--blue','LineWidth',2)
% ylim([0 1])
% plot(x,p_i,'green','Linewidth',2)
% xlabel('x')
% ylabel('Probability')
% hold on
% plot(x,m_i,'green--','LineWidth',2)
% ylim([0 1])
% figure()
% plot(0:30,C_analytical)
% hold on
% plot(0:30,C,'--')
% legend('Analytical','Numerical')
% xlabel('q_p (cell cycle^{-1})')
% ylabel('Wave speed c')
% %-----------------------------------
% %SENSITIVITY ANALYSIS
% %clear all,close all,clc
% 
% v=5;
% alpha=1;
% mu=1e-3;
% beta=10;
% q_p=22;
% 
% figure
% plot(0:0.1:10,C)%numerical wave's speed vs q_m
% hold on
% plot(0:0.1:10,Ca)%numerical wave's speed vs q_m
% xlabel('q_m')
% ylabel('wave speed c')
% %q_m=6.3;
% % for q_p=0:0.1:10
% %     i=round(q_p*10+1,0);
% %     C(i)=GerleeContinuumModelFun(q_m,q_p,alpha,v,mu);
% % end
% % figure
% % plot(0:0.1:10,C)%wave's speed vs q_m
% % xlabel('q_m')
% % ylabel('wave speed c')
% 
