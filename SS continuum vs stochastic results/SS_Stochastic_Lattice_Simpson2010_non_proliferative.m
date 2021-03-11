%STOCHASTIC VS CONTINUUM - SINGLE-SPECIES

%CONTINUUM MODEL

clear all,close all,clc
%NON DIMNESIONALIZE (TIME AND SPACE)
% C=[]; %for all wave speeds
% C_analytical=[];
% C_analytical2=[];
% Xf=[];
% Xi=[];
% DT=[0.01, 0.005, 0.001];
% DX=[0.01,0.005,0.001];
% %     q_p
% % for v=0:8
% %PARAMETERS
% v=1;
% mu=0;
% alpha=0;
% %Initial condiitons for stochastic
% N1=20;
% N2=200;
% Chain=zeros(N1,N2);
% % Chain(round(N1/2),round(N2/2))=1;
% % C_0=0.6;
% % 
% % InitialField=Chain;
% % InitialField=Chain(:,181:220);
% % while nnz(InitialField)/numel(InitialField)<C_0
% %     Site=randi(numel(InitialField));
% %     InitialField(Site)=1;%get position in the area of initial cells
% %     start=sub2ind(size(Chain),1,180); %row should be 1
% %     Chain(Site+start)=1;%populate original matrix properly
% % end
% %CALCULATE ANALYTICAL TRAVELLING WAVE SPEEDS
% c_analytical=sqrt(2*(alpha-mu)*(alpha+v));%from point (0,0)
% C_analytical=[C_analytical c_analytical];
% % c_analytical2=sqrt(2*(mu+v)*(mu+alpha));%from point (1-mu/alpha,0)
% % C_analytical2=[C_analytical2 c_analytical2];
% %NUMERICAL PARAMETERS
% h=1;
% D_alpha=(h^2)*alpha/4;
% D_v=(h^2)*v/4;
% dx=0.1;
% dt=0.1;
% x_inf=200;
% x=0:dx:x_inf;
% epsilon=1e-6;
% %Continuum initial conditions
% % a=exp(-beta*x);
% D=1/4;
% x0=40;
% a=ones(1,length(x));
% for i=1:length(a)
%     if i<=find(x==100)-20/dx || i>=find(x==100)+20/dx
%         a(i)=0;
%     end
% end
% % a=heaviside(x+10)-heaviside(x-8);
% % a=(C_0/2)*((erf(-(x-x0/2)/sqrt(4*D*dt)))+(erf((x+x0/2)/sqrt(4*D*dt))));
% a0=a;
% a_previous=a;
% A=zeros(1,length(x));
% d1_a=zeros(1,length(x)-1);
% d2_a=zeros(1,length(x));
% d3_a=zeros(1,length(x)-1);
% CellCycles=1000;
% for j=0:dt:CellCycles %time iterations
%     [mu dt dx j]
% %     clear d1_p d2_p d3_p d1_m d2_m d3_m P M %clear diagonal values at each cycle
% 
%     if j-dt<500 && j+dt>500 && abs(floor(j-dt)-floor(j))==1
%         a1000=a;
%     end
%     if j-dt<1000 && j+dt>1000 && abs(floor(j-dt)-floor(j))==1
%         a2000=a;
%     end
% 
% 
%     if j==15 %evaluate wave's speed using two snap shots
%         a_i=a;
%     elseif j==CellCycles
%         a_f=a;
%         %get point close to 1/2
%         x_i=min(find(abs(a_i-max(a_i)/2)<=(0.01)));
%         x_f=min(find(abs(a_f-max(a_f)/2)<=(0.01)));
%         c=dx*(x_f-x_i)/(j-15);%wave speed
%         Xf=[Xf x_f];
%         Xi=[Xi x_i];
%         C=[C c];
%     end
% 
% %     j
%     da=0.0001*ones(length(x),1);
% 
%     while norm(da,Inf)>epsilon
% %            norm(da)
%         
%        %---------------------------
%        %Start Boundary conditions for i=1
%         A(1)=(a(2)-a(1));
%         d1_a(1)=0;%dF/da_iminus1=0
%         d2_a(1)=-1;%dF/da_i=-1
%         d3_a(1)=1;%dF/da_iplus1=1
%         
%         for i=2:length(x)-1
%             
%             %define algebraic nonlinear equation to solve
%             A(i)=-(a(i)-a_previous(i))/dt+(D_alpha/dx^2*(1-a(i))+D_v/dx^2)*(a(i-1)-2*a(i)+a(i+1))+alpha*a(i)*(1-a(i))-mu*a(i);
%             if isnan(A(1))
%                 display('stop')
%             end
%            
% %             d1_a(i)=(D_v/(dx^2));
%             d1_a(i)=(D_alpha/(dx^2))*(1-a(i))+(D_v/(dx^2));
%             d2_a(i)=-1/dt-2*((D_v/dx^2))+alpha*(1-a(i))-alpha*a(i)-mu;
%             d3_a(i)=(D_alpha/(dx^2))*(1-a(i))+(D_v/(dx^2));
% %             d3_a(i)=(D_v/(dx^2));
%             
%         end
%         
%         % End Boundary Conditions
%         d1_a(end)=0;
%         d2_a(end)=-1;
%         d3_a(end)=1;
%         A(end)=(a(end)-a(end-1));
%         %-----------------------------------------
%         d1_a=[0 d1_a];
%         d3_a=[d3_a 0];
%         da=tridiag(d1_a',d2_a',d3_a',-A);
% %        dm=Thomas(d1_a',d2_a',d3_a',-A);
%         if isnan(da(1))
%             display('stop')
%         end
%         %-------------------------------------------
%         a=a+da';       
%     end   
%       %keep record of previous u
%      a_previous=a;
% %      plot(x,a)
% %      ylim([0 1])
% %      xlim([min(x) max(x)])
% %      pause(0.00001)
% end


% C=[];
% C_analytical=[];
%%%--------------------------------------------
%SINGLE-SPECIES STOCHASTIC MODEL
%Single cell model
%numerical parameters
N1=20;
N2=200;
CellCycles=1000;
%model parameters
alpha=0;
v=1;
mu=0;
TumorSize_singleCell=[];

%multiple simulations
% T=25;
% % x=-25:25;
tau=1;
% Chain_avg=zeros(N1,N2);
% Chain_avg0=zeros(N1,N2);
% Chain_avg500=zeros(N1,N2);
% Chain_avg1000=zeros(N1,N2);
% Chain_avg2000=zeros(N1,N2);
% 
Chain_avg2=zeros(N1,N2);
Chain_avg02=zeros(N1,N2);
Chain_avg5002=zeros(N1,N2);
Chain_avg10002=zeros(N1,N2);
Chain_avg20002=zeros(N1,N2);
% 
% Chain_avg2=zeros(N1,N2);
% WaveSpeed_avg=0;
%     %Calculate Analytical Wave Speed
%     c_analytical=sqrt(2*(alpha-mu)*(alpha+v));%from point (0,0)
%     C_analytical=[C_analytical c_analytical];
%     Density_history_avg=zeros(1,T+1);
%     %initial conditions
%     Chain=zeros(N1,N2);
%     for i=1:size(Chain,1)
%         for j=81:120
%             Chain(i,j)=1;
%         end
%     end
    I=1;
    for i=1:I
        i
% % %set cells
% % Chain=zeros(N1,N2);
% % Chain(round(N1/2),round(N2/2))=1;
Chain=zeros(N1,N2);
% % InitialField=Chain;
C_0=0.6;
% Chain(round(N1/2),round(N2/2))=1;
% InitialField=Chain;
InitialField=Chain(:,81:120);
while nnz(InitialField)/numel(InitialField)<C_0
    Site=randi(numel(InitialField));
    InitialField(Site)=1;%get position in the area of initial cells
    start=sub2ind(size(Chain),1,80); %row should be 1
    Chain(Site+start)=1;%populate original matrix properly
end
% 
% 
% %models
% % [Chain_out,Chain02,Chain5002,Chain10002,Chain20002,Cells,TumorSize_cellcycles_singlecell,WaveSpeed,Q_ij,Cells_index_vect,Density_history,c1,C1]=SingleCellSim(T,alpha,v,mu,Chain);
[Chain_out2,Chain02,Chain5002,Chain10002,Chain20002,Cells_index_vect,WaveSpeed,Density_history,c2,C2,WaveSpeed_vect] = SingleCellSim2(CellCycles,alpha,v,Chain,tau);
% % TumorSize_singleCell=[TumorSize_singleCell nnz(Chain)];
% % Chain_avg=Chain_avg+Chain_out;
% % Chain_avg0=Chain_avg0+Chain;
% % Chain_avg500=Chain_avg500+Chain500;
% % Chain_avg1000=Chain_avg1000+Chain1000;
% % Chain_avg2000=Chain_avg2000+Chain2000;
% 
Chain_avg2=Chain_avg2+Chain_out2;
Chain_avg02=Chain_avg02+Chain;
Chain_avg5002=Chain_avg5002+Chain5002;
Chain_avg10002=Chain_avg10002+Chain10002;
Chain_avg20002=Chain_avg20002+Chain20002;
    end
% % Chain_avg=Chain_avg/I;
% % Chain_avg0=Chain_avg0/I;
% % Chain_avg500=Chain_avg500/I;
% % Chain_avg1000=Chain_avg1000/I;
% % Chain_avg2000=Chain_avg2000/I;
% 
Chain_avg2=Chain_avg2/I;
Chain_avg02=Chain_avg02/I;
Chain_avg5002=Chain_avg5002/I;
Chain_avg10002=Chain_avg10002/I;
Chain_avg20002=Chain_avg20002/I;
% % g0=mean(Chain_avg0);
% % g500=mean(Chain_avg500);
% % g1000=mean(Chain_avg1000);
% % g2000=mean(Chain_avg2000);
% % gfinal=mean(Chain_avg);
% 
g02=mean(Chain);
g5002=mean(Chain_avg5002);
g10002=mean(Chain_avg10002);
g20002=mean(Chain_avg20002);
gfinal2=mean(Chain_avg2);
% % g0=g0(100:300);
% % g500=g500(100:300);
% % g1000=g1000(100:300);
% % Chain_avg2=Chain_avg2/I;
% %Continuum for SS2
% tau=0.6;
% h=1;
% D=v*h^2/(4*tau);
% % 
% % D=D_v;
% C_0=0.01;
% x0=40;
% % g=Chain_avg(10,:);
% T=3000;
% x=-200:200;
% Continuum_approx_nonProlif=(C_0/2)*((erf(-(x-x0/2)/sqrt(4*D*T)))+(erf((x+x0/2)/sqrt(4*D*T))));
% 
% Continuum_approx0=(C_0/2)*((erf(-(x-x0/2)/sqrt(4*D*1000)))+(erf((x+x0/2)/sqrt(4*D*1000))));
% Continuum_approx500=(C_0/2)*((erf(-(x-x0/2)/sqrt(4*D*2000)))+(erf((x+x0/2)/sqrt(4*D*2000))));
% Continuum_approx1000=(C_0/2)*((erf(-(x-x0/2)/sqrt(4*D*3000)))+(erf((x+x0/2)/sqrt(4*D*3000))));

% %DENSITY PLOTS
% g2=mean(Chain_avg2);
% % g=g(100:200);
% g2=g2(100:200);

% plot(0:101/(length(g)):100,g,'Linewidth',2)
% hold on
% plot(0:101/(length(g2)):100,g2,'Linewidth',2)
% hold on
% plot(0:400/(length(a0)-1):400,a0,'blue','LineWidth',2)
% hold on
% plot(0:400/(length(a1000)-1):400,a1000,'blue','LineWidth',2)
% hold on
% plot(0:400/(length(a2000)-1):400,a2000,'blue','LineWidth',2)
% hold on
% plot(0:400/(length(a)-1):400,a,'blue','LineWidth',2)
% hold on
% plot(0:401/(length(g0)):400,g0,'red','LineWidth',1)
% hold on
% plot(0:401/(length(g1000)):400,g1000,'red','LineWidth',1)
% hold on
% plot(0:401/(length(g2000)):400,g2000,'red','LineWidth',1)
% hold on
% plot(0:401/(length(gfinal)):400,gfinal,'red','LineWidth',1)
% hold on

%Lattice plots
subplot(3,1,1)
imagesc(Chain)
set(gca,'XTick',0:20:200)
set(gca,'XTickLabel',-100:20:100)
colormap([1 1 1;0 0 0;1 0 0])
ylabel('y')
xlabel('x')
subplot(3,1,2)
imagesc(Chain5002)
set(gca,'XTick',0:20:200)
set(gca,'XTickLabel',-100:20:100)
colormap([1 1 1;0 0 0;1 0 0])
ylabel('y')
xlabel('x')
subplot(3,1,3)
imagesc(Chain10002)
set(gca,'XTick',0:20:200)
set(gca,'XTickLabel',-100:20:100)
colormap([1 1 1;0 0 0;1 0 0])
ylabel('y')
xlabel('x')

%Continuum
% plot(x,a0,'LineWidth',2)
% hold on
% plot(x,a1000,'LineWidth',2)
% hold on
% plot(x,a2000,'LineWidth',2)
% hold on
% plot(x,a,'LineWidth',2)
% hold on
% %Stochastic
% plot(0:201/(length(g02)):200,g02,'blue','LineWidth',1)
% hold on
% plot(0:201/(length(g10002)):200,g10002,'blue','LineWidth',1)
% hold on
% plot(0:201/(length(g20002)):200,g20002,'blue','LineWidth',1)
% hold on
% plot(0:201/(length(gfinal2)):200,gfinal2,'blue','LineWidth',1)
% xlabel('x')
% ylabel('Cell density')



% legend(L,{'Gillespie','Approximation'})
% [row col]=ind2sub([N1 N2],Cells_index_vect);
% scatter(col, row,'red','filled')
% plot(0:T,Density_history_avg,'Linewidth',2)
% % hold on
% plot(x,Continuum_approx0,'Linewidth',2)
% hold on
% plot(x,Continuum_approx500,'Linewidth',2)
% hold on
% plot(x,Continuum_approx1000,'Linewidth',2)
% hold on
% plot(x,g0)
% hold on
% plot(x,g500)
% hold on
% plot(x,g1000)
% xlabel('x')
% ylabel('Probability')
% legend('T=0: Continuum','T=500: Continuum','T=1000: Continuum','T=0: Stochastic','T=500: Stochastic','T=1000: Stochastic')
% plot(x,g,'Linewidth',2)

% title(['T=' num2str(j*dt)])
% ylabel('Probability')
% xlabel('t')
% legend('Stochastic','Continuum')
% ylim([0 1])
% hold off
% hold off
% g=Chain_avg(100,:);
% plot(0:100,g(100:200),'--','Linewidth',2)
% hold on
% plot(x,a,'LineWidth',2)
% legend('Stochastic','Continuum')
% hold off
% ylim([0 1])
% pause(0.00001)
% end
% g=Chain_avg(:,100);
% plot(0:100,g(100:200),'black','Linewidth',2)
% hold on
% plot(x,a,'--','LineWidth',2)
% hold off
% ylim([0 1])
