function [c,C_analytical,p,m] = GerleeContinuumFunction(x,dx,dt,p,m,CellCycles,alpha,v,mu,q_p,q_m)
%GERLEE CONTINUUM MODEL

%NON DIMNESIONALIZE (TIME AND SPACE)
C=[]; %for all wave speeds
c=0;
C_analytical=[];
Xf=[];
Xi=[];
beta=10;

if CellCycles==0
    p=p';
    m=m';
    return
end

%CALCULATE ANALYTICAL TRAVELLING WAVE SPEEDS
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];
%initial conditions
%    Chain=zeros(N1,N2);
%     Chain(round(N1/2),round(N2/2))=1; %single P-cell

%NUMERICAL PARAMETERS
h=1;
D_alpha=(h^2)*alpha/2;
D_v=(h^2)*v/2;

%intial conditions2
% m=0.5*ones(1,length(x));
% p=0.5*ones(1,length(x));
% for i=1:length(x)
%     if i<=find(x==200)-20/dx || i>=find(x==200)+20/dx
%         p(i)=0;
%         m(i)=0;
%     end
% end
p0=p;
m0=m;
epsilon=1e-8;
% p=exp(-beta*x);
% p=0.5*heaviside(-x+1);
if size(p,2)>1
    p=p';
    m=m';
end
% m=zeros(length(x),1);
p_previous=p;
m_previous=m;
j40_check=1;%for the snapshots to be at the right time point
j50_check=1;
for j=0:dt:CellCycles %time iterations
%     t=dt*j;
    j
%     clear d1_p d2_p d3_p d1_m d2_m d3_m P M %clear diagonal values at each cycle
%     if abs(j-40)<j40_check %evaluate wave's speed using two snap shots
%         j40_check=abs(j-40);
%         p_i=p;
%         m_i=m;
%         pm_i=p_i+m_i; %sum of m and p
%     elseif abs(j-50)<abs(j-50)
%         j50_check=j;
%         p_f=p;
%         m_f=m;
%         pm_f=p_f+m_f; %sum of m and p
%         %get point close to 1/2
%         x_i=min(find(abs(pm_i-max(pm_i)/2)<=(0.005)));
%         x_f=min(find(abs(pm_f-max(pm_f)/2)<=(0.005)));
%         c=(dx*x_f-dx*x_i)/(50-40);%wave speed
%         Xf=[Xf x_f];
%         Xi=[Xi x_i];
%         C=[C c];
%     end

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
                %return
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

        d1_p(length(x))=0;
        d2_p(length(x))=-1;
        d3_p(length(x))=1;
        d1_m(length(x))=0;
        d2_m(length(x))=-1;
        d3_m(length(x))=1;
        P(length(x))=p(length(x))-p(length(x)-1);
        M(length(x))=m(length(x))-m(length(x)-1);
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

%WAVE SPEED---------------------------------
p_i=p0;
m_i=m0;
pm_i=p_i+m_i; %sum of m and p

p_f=p;
m_f=m;
pm_f=p_f+m_f; %sum of m and p
%get point close to 1/2
x_i=1800;
x_f=min(find(abs(pm_f-max(pm_f)/2)<=(0.01)));
c=abs((dx*x_f-dx*x_i))/CellCycles;%wave speed
if isempty(c)==1
   c=0; 
end
%Wave Slope---------------------------------
%using the final snapshot to get the slope

%gradient()


end
