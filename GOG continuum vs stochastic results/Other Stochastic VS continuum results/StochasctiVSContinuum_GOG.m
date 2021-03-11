%GERLEE CONTINUUM MODEL-----------------------------------------------
%---------------------------------------------------------------------
clear all
clc
%PARAMETERS
v=1;
alpha=1;
mu=0.001;
beta=10;
q_p=10;
q_m=0;
N1=200;
N2=200;
%initial conditions
%    Chain=zeros(N1,N2);
%     Chain(round(N1/2),round(N2/2))=1; %single P-cell

%NUMERICAL PARAMETERS
h=1;
D_alpha=(h^2)*alpha/2;
D_v=(h^2)*v/2;
dx=0.1;
dt=0.1;
x_inf=100;
x=0:dx:x_inf;
%intial conditions2
% m=0.5*ones(1,length(x));
% p=0.5*ones(1,length(x));
% for i=1:length(x)
%     i
%     if i<=find(x==200)-20/dx || i>=find(x==200)+20/dx
%         p(i)=0;
%         m(i)=0;
%     end
% end
% p0=p;
% m0=m;
epsilon=1e-6;
p=exp(-beta*x);
% p=0.5*heaviside(-x+1);
m=zeros(1,length(x));
p=p';
m=m';
p_previous=p;
m_previous=m;
P=zeros(1,length(x));
M=zeros(1,length(x));
d1_p=zeros(1,length(x)-1);
d2_p=zeros(1,length(x));
d3_p=zeros(1,length(x)-1);
d1_m=zeros(1,length(x)-1);
d2_m=zeros(1,length(x));
d3_m=zeros(1,length(x)-1);
CellCycles=20;
for j=0:dt:CellCycles %time iterations
    j
    if j-dt<=10 && j+dt>=10
        p50=p;
        m50=m;
    end
    if j-dt<=50 && j+dt>=50
        p100=p;
        m100=m;
    end
%     clear d1_p d2_p d3_p d1_m d2_m d3_m P M %clear diagonal values at each cycle
    m_u_up=0;
    dp=0.0001*ones(length(x),1);
    dm=0.0001*ones(length(x),1);
    
%     if length(P)<length(x)
%         somethingtostop
%     end
    
    while norm(dp,Inf)>epsilon || norm(dm,Inf)>epsilon
%         norm(dp)
%         norm(dm)
       %---------------------------
       %Start Boundary conditions for i=1
        P(1)=(p(2)-p(1));
        M(1)=(m(2)-m(1));
        d1_p(1)=0;%dF/du_iminus1=0
        d2_p(1)=-1;%dF/du_i=-1
        d3_p(1)=1;%dF/du_iplus1=1
        d1_m(1)=0;%dF/du_iminus1=0
        d2_m(1)=-1;%dF/du_i=-1
        d3_m(1)=1;%dF/du_iplus1=1

        for i=2:length(x)-1
            
            %define algebraic nonlinear equation to solve
% %             IMPLICIT
            d2pdx2=(p(i-1)+p(i+1)-2*p(i))/dx^2;
            d2mdx2=(m(i-1)+m(i+1)-2*m(i))/dx^2;
%             %Algebraic equation
            P(i)=-(p(i)-p_previous(i))/dt+(D_alpha)*(1-p(i)-m(i))*(d2pdx2)+alpha*p(i)*(1-p(i)-m(i))-(q_m+mu)*p(i)+q_p*m(i);
            M(i)=-(m(i)-m_previous(i))/dt+(D_v)*((1-p(i))*(d2mdx2)+m(i)*(d2pdx2))-(q_p+mu)*m(i)+q_m*p(i);
            if isnan(P(i))||isnan(M(i))
                display('stop')
            end
            
            d1_p(i)=(D_alpha/(dx^2))*(1-p(i)-m(i));
            d2_p(i)=-1/dt+(D_alpha/(dx^2))*(-p(i-1)+4*p(i)-p(i+1)+2*m(i)-2)+alpha*(1-2*p(i)-m(i))-(q_m+mu);
            d3_p(i)=(D_alpha/(dx^2))*(1-p(i)-m(i));
            d1_m(i)=(D_v/(dx^2))*(1-p(i)-m(i));
            d2_m(i)=-1/dt+(D_v/(dx^2))*(-2+p(i-1)+p(i+1))-(q_p+mu);
            d3_m(i)=(D_v/(dx^2))*(1-p(i)-m(i));
            
%             P(i)=-(p(i)-p_previous(i))/dt+alpha*p(i)*(1-p(i)-m(i))-(q_m+mu)*p(i)+q_p*m(i);
%             M(i)=-(m(i)-m_previous(i))/dt+(D_v/dx^2)*((1-p(i))*(m(i-1)+m(i+1)-2*m(i))+m(i)*(p(i-1)+p(i+1)-2*p(i)))-(q_p+mu)*m(i)+q_m*p(i);
%             if isnan(p_previous(i))||isnan(m_previous(i))
%                 display('stop')
%                 break
%             end            
%             d1_p(i)=0;
%             d2_p(i)=-1/dt+alpha*(1-2*p(i)-m(i))-(q_m+mu);
%             d3_p(i)=0;
%             d1_m(i)=(D_v/(dx^2))*(1-p(i));
%             d2_m(i)=-1/dt+(D_v/(dx^2))*(-2+p(i-1)+p(i+1))-(q_p+mu);
%             d3_m(i)=(D_v/(dx^2))*(1-p(i));
        end
        % End Boundary Conditions

        d1_p(end)=0;
        d2_p(end)=-1;
        d3_p(end)=1;
        d1_m(end)=0;
        d2_m(end)=-1;
        d3_m(end)=1;
        P(end)=(p(end)-p(end-1));
        M(end)=(m(end)-m(end-1));
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
            disp('stop')
           break
        end
        %-------------------------------------------
        p=p+dp;
        m=m+dm;
        
        
    end   
      %keep record of previous u
     p_previous=p;
     m_previous=m;
     
%      plot(x,p+m,'blue')
%      xlabel('x')
%      ylabel('p + m')
%      hold on
%      ylim([0 1])
%      xlim([min(x) max(x)])
%      hold off
%      pause(0.00001)
     
%%%---------------------------------------------
%GERLEE STOCHASTIC MODEL-----------------------------------------------
%---------------------------------------------------------------------
Pcells_FinalTime=[];%make list of final population
Mcells_FinalTime=[];
N1=200;
N2=200;
%for p and m
Chain_avg=zeros(N1,N2);%Chain to get the overall average of simulations
Chain_avg50=zeros(N1,N2);%Chain to get the overall average of simulations
Chain_avg100=zeros(N1,N2);%Chain to get the overall average of simulations
%for p
Chain_avg_p=zeros(N1,N2);%Chain to get the overall average of simulations
Chain_avg_p0=zeros(N1,N2);%Chain to get the overall average of simulations
Chain_avg_p50=zeros(N1,N2);%Chain to get the overall average of simulations
Chain_avg_p100=zeros(N1,N2);%Chain to get the overall average of simulations
%for m
Chain_avg_m0=zeros(N1,N2);%Chain to get the overall average of simulations
Chain_avg_m=zeros(N1,N2);%Chain to get the overall average of simulations
Chain_avg_m50=zeros(N1,N2);%Chain to get the overall average of simulations
Chain_avg_m100=zeros(N1,N2);%Chain to get the overall average of simulations

Total_cell_avg=zeros(20,20);%matrix to study q rates with average number of cells
Total_cells=zeros(20,20);%matrix to study q_p and q_m parameters with total cells
WaveSpeed_vector=[];
    C_analytical=[];
    I=100; %number of sims
    WaveSpeed_avg=0;
    %     [j*dt i]
 density_cellcycles_avg=zeros(1,round(j*dt)+1);
for i=1:I
% i
t=0;
% dim=size(Chain);
M_vector=[];
P_vector=[];
random_indices_P=[];
random_indices_M=[];
%initial conditions
    Chain=zeros(N1,N2);
    Chain(N1/2,N2/2)=1;
%     for o=1:size(Chain,1)
%         for j=180:220
%             r=rand();
%             if r<0.5
%                 Chain(o,j)=1;
%             else
%                 Chain(o,j)=2;
%             end
%         end
%     end
[M2,P2,Chain_out,Chain50,Chain100,Chain200,Q_ij,WaveSpeed,density_cellcycles,WaveSpeed_vect]=Gerlee_Main_loop(Chain,alpha,v,q_p,q_m,mu,j,0);
% imagesc(Chain)
% title(['Sim=' num2str(i)])
% pause(0.1)
% density_cellcycles_avg=density_cellcycles_avg+density_cellcycles;
Mcells_FinalTime=[Mcells_FinalTime M2];
Pcells_FinalTime=[Pcells_FinalTime P2];

%prepare matrices for the marginal results of m and p
Chain_p0=zeros(N1,N2);
Chain_m0=zeros(N1,N2);
Chain_p50=zeros(N1,N2);
Chain_m50=zeros(N1,N2);
Chain_p100=zeros(N1,N2);
Chain_m100=zeros(N1,N2);
Chain_p=zeros(N1,N2);%final result
Chain_m=zeros(N1,N2);%final result
for l=1:numel(Chain_out)%Construt marginal population plots for densities
    %for initial time
%     if Chain(l)==1%proliferative
%         Chain_p0(l)=1;
%         Chain_m0(l)=0;
%     elseif Chain(l)==2%motile
%         Chain_p0(l)=0;
%         Chain_m0(l)=1;
%     elseif Chain(l)==0
%         Chain_p0(l)=0;
%         Chain_m0(l)=0;
%     end
%     %for final time
%     if Chain_out(l)==1%proliferative
%         Chain_p(l)=1;
%         Chain_m(l)=0;
%     elseif Chain_out(l)==2%motile
%         Chain_p(l)=0;
%         Chain_m(l)=1;
%     elseif Chain_out(l)==0
%         Chain_p(l)=0;
%         Chain_m(l)=0;
%     end
%     %for t=50
%     if Chain50(l)==1%proliferative
%         Chain_p50(l)=1;
%         Chain_m50(l)=0;
%     elseif Chain50(l)==2%motile
%         Chain_p50(l)=0;
%         Chain_m50(l)=1;
%     elseif Chain50(l)==0
%         Chain_p50(l)=0;
%         Chain_m50(l)=0;
%     end
%     %for t=100
%     if Chain100(l)==1%proliferative
%         Chain_p100(l)=1;
%         Chain_100(l)=0;
%     elseif Chain100(l)==2%motile
%         Chain_p100(l)=0;
%         Chain_m100(l)=1;
%     elseif Chain100(l)==0
%         Chain_p100(l)=0;
%         Chain_m100(l)=0;
%     end
end
%for p
Chain_avg_p0=Chain_avg_p0+Chain_p0;
Chain_avg_p=Chain_avg_p+Chain_p;
Chain_avg_p50=Chain_avg_p50+Chain_p50;
Chain_avg_p100=Chain_avg_p100+Chain_p100;
%for m
Chain_avg_m0=Chain_avg_m0+Chain_m0;
Chain_avg_m=Chain_avg_m+Chain_m;
Chain_avg_m50=Chain_avg_m50+Chain_m50;
Chain_avg_m100=Chain_avg_m100+Chain_m100;
%for both
Chain_avg=Chain_avg+Matrix2Binary(Chain_out);%Matrix2Binary is required since there are 1's and 2's on the matrix
Chain_avg50=Chain_avg50+Matrix2Binary(Chain50);
Chain_avg100=Chain_avg100+Matrix2Binary(Chain100);
% Total_cell_avg(find(Qp==q_p),find(Qm==q_m))=Total_cell_avg(find(Qp==q_p),find(Qm==q_m))+M2+P2;
WaveSpeed_avg=WaveSpeed_avg+WaveSpeed;
end
density_cellcycles_avg=density_cellcycles_avg/I;
% Chain=Chain_out;
Chain_avg=Chain_avg/I;
Chain_avg50=Chain_avg50/I;
Chain_avg100=Chain_avg100/I;
%for p
Chain_avg_p0=Chain_avg_p0/I;
Chain_avg_p=Chain_avg_p/I;
Chain_avg_p50=Chain_avg_p50/I;
Chain_avg_p100=Chain_avg_p100/I;
%for m
Chain_avg_m0=Chain_avg_m0/I;
Chain_avg_m=Chain_avg_m/I;
Chain_avg_m50=Chain_avg_m50/I;
Chain_avg_m100=Chain_avg_m100/I;
% imagesc(Chain_avg)
% title(['Sim=' num2str(i) ' ;T=' num2str(j*dt)])
% colorbar()
% pause(0.1)
% Chain=Chain_avg;
% Total_cell_avg(find(Qp==q_p),find(Qm==q_m))=Total_cell_avg(find(Qp==q_p),find(Qm==q_m))/I;
%DENSITY PLOTS
g=Chain_avg(100,:);
plot(0:100,g(100:200),'black','Linewidth',2)
hold on
plot(x,p+m,'--','LineWidth',2)
hold off
ylim([0 1])
legend('Stochastic','Continuum')
pause(0.00001)
end
%DENSITY PLOTS
% plot(density_cellcycles_avg)
% g=Chain_avg(:,100);
g=mean(Chain_avg);
g0=mean(Chain);
g50=mean(Chain_avg50);
g100=mean(Chain_avg100);
%for p
g_p=mean(Chain_avg_p);
g0_p=mean(Chain_avg_p0);
g50_p=mean(Chain_avg_p50);
g100_p=mean(Chain_avg_p100);
%for m
g_m=mean(Chain_avg_m);
g0_m=mean(Chain_avg_m0);
g50_m=mean(Chain_avg_m50);
g100_m=mean(Chain_avg_m100);
%Plots
% figure()
% plot(x,p0+m0,'blue','LineWidth',2)
% hold on
% plot(x,p50+m50,'blue','LineWidth',2)
% hold on
% plot(x,p100+m100,'blue','LineWidth',2)
% hold on
% plot(x,p+m,'blue','LineWidth',2)
% hold on
% plot(0:401/(length(g0)):400,g0,'red','Linewidth',1)
% hold on
% plot(0:401/(length(g50)):400,g50,'red','Linewidth',1)
% hold on
% plot(0:401/(length(g100)):400,g100,'red','Linewidth',1)
% hold on
% plot(0:401/(length(g)):400,g,'red','Linewidth',1)
% ylim([0 1])
% xlabel('x')
% ylabel('Probability of occupancy')
% title('Total population')
% hold off
% hold off
% hold off
% hold off
% hold off
% hold off
% hold off
%----------------------------------------------------------------
%plots for p and m separate
% %P cells 
% figure()
% L(1)=plot(x,p0,'blue','Linewidth',2);
% hold on
% L(2)=plot(x,p50,'blue','Linewidth',2);
% hold on
% L(3)=plot(x,p100,'blue','Linewidth',2);
% hold on
% L(4)=plot(x,p,'blue','Linewidth',2);
% hold on
% plot(0:401/(length(g0_p)):400,g0_p,'red','Linewidth',1)
% hold on
% plot(0:401/(length(g50_p)):400,g50_p,'red','Linewidth',1)
% hold on
% plot(0:401/(length(g100_p)):400,g100_p,'red','Linewidth',1)
% hold on
% plot(0:401/(length(g_p)):400,g_p,'red','Linewidth',1)
% xlabel('x')
% ylabel('Probability of occupancy')
% ylim([0 1])
% title('Proliferative cells')
% %For m cells
% figure()
% L(1)=plot(x,m0,'blue','LineWidth',2);
% hold on
% L(2)=plot(x,m50,'blue','LineWidth',2)
% hold on
% L(3)=plot(x,m100,'blue','LineWidth',2)
% hold on
% L(4)=plot(x,m,'blue','LineWidth',2)
% hold on
% plot(0:401/(length(g0_m)):400,g0_m,'red','Linewidth',1)
% hold on
% plot(0:401/(length(g50_m)):400,g50_m,'red','Linewidth',1)
% hold on
% plot(0:401/(length(g100_m)):400,g100_m,'red','Linewidth',1)
% hold on
% plot(0:401/(length(g_m)):400,g_m,'red','Linewidth',1)
% % legend(L,{'t=0','t=10','t=50','t=100'})
% xlabel('x')
% ylabel('Probability of occupancy')
% ylim([0 1])
% title('Motile cells')
