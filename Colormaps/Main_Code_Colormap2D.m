% profile on
%Gillespie algorithm for Gerlee's IB model
clear all
clc
%Initialize number of molecules in system, reaction constants and random
%number generators
Pcells_FinalTime=[];%make list of final population
Mcells_FinalTime=[];
N1=200;
N2=200;
Chain_avg=zeros(N1,N2);
Chain_avg2=zeros(N1,N2);%Chain to get the overall average of simulations
Chain=zeros(N1,N2);
% InitialField=Chain;
% C_0=0.01;
% while nnz(InitialField)/numel(InitialField)<C_0
%     
%     Site=randi(numel(InitialField));
%     InitialField(Site)=1;%get position in the area of initial cells
%     Chain(Site)=1;%populate original matrix properly
% end
Total_cell_avg=zeros(20,20);%matrix to study q rates with average number of cells
Total_cells=zeros(20,20);%matrix to study q_p and q_m parameters with total cells
WaveSpeed_vector=[];
figure
for k=1
    C=[];
    c=0;
    C_analytical=[];
    I=10; %number of simulations
    T=50;
    Qp=0:1.55:30;
    Qm=0:1.55:30;
%     Qp=20;
%     Qm=10;
for  q_m=Qm
for q_p=Qp
    WaveSpeed_avg=0;
    WaveSpeed=0;
    Density_history_avg=zeros(1,T+1);
for i=1:I
%     i
    [i q_p q_m WaveSpeed_avg WaveSpeed]
% figure()
% for q_m=1:30
Chain=zeros(N1,N2);

%circular initial conditions (Full circle)---------------------------------
% Chain=Full_circle(Chain,N1);
% Chain=2*Chain;%all m cells
%circular initial conditions (holed circle)-----------------------------
% Chain=Holed_circle(Chain,N);
Chain(round(N1/2),round(N2/2))=1; %single P-cell
% C=[10 10;110 30; 50 90]
% C=[round(N1/2),round(N1/2)];
% R=ones(size(C,1),1)*5;%radius of circles
% pre_Chain=MakeCircles(Chain,N1,N2,C,R);
% for i=1:numel(pre_Chain)
%     if pre_Chain(i)~=0
% %           Chain(i)=2;
%         Chain(i)=1+round(rand()); %randomly assign motile/proliferative to each position of initial condition
%     end
% end


%      figure(1)
% imagesc(Chain)
% colorbar('Ticks',[0,1,2],...
%          'TickLabels',{'Empty Space','Proliferating','Moving'})
%Right Strip
% Chain(:,round(N/2)+33:end)=1;
%Middle Strip
% Chain(:,round(N/2-8.5):round(N/2+8.5))=1; %middle strips
%Left-Right Strips (Full Scratch Assay)
% Chain(:,1:round(N/2)-41)=1; 
% Chain(:,round(N/2)+42:end)=1; 
%Left-Right Strips (Random Scratch Assay)
% Chain=HoledScratch(Chain,N1,N2); 

M=0;%initial number of M cells (random walk)
% P=0;
P=1;%initial number of P cells (stationary and proliferating)
alpha=1;
v=5;
% q_p=15;%rate from m to p
% q_m=15;%rate from p to m
mu=1e-3;%death rate
%Monte Carlo: Generate random numbers to determine the next reaction to occur as well as
%the time interval. The probability of a given reaction to be chosen is
%proportional to the number of substrate molecules, the time interval is
%exponentially distributed with mean q_p or q_m (depending on the species, both <24)
t=0;
% T=20;
dim=size(Chain);
M_vector=[];
P_vector=[];
random_indices_P=[];
random_indices_M=[];
[M,P,Chain_out,Chain50,Chain100,Chain200,Q_ij,WaveSpeed,TumorSize_cellcycles,WaveSpeed_vect]=Gerlee_Main_loop(Chain,alpha,v,q_p,q_m,mu,T,t);
%in case a cell cycle was missed for the density history, repeat the
%simulation
% while length(Density_history_avg)~=length(Density_history)
%     [M,P,Chain_out,Q_ij,WaveSpeed,Density_history]=Gerlee_Main_loop(Chain,alpha,v,q_p,q_m,mu,T,t);
% end
%     Density_history_avg=Density_history_avg+Density_history;
Mcells_FinalTime=[Mcells_FinalTime M];
Pcells_FinalTime=[Pcells_FinalTime P];
Chain_avg=Chain_avg+Matrix2Binary(Chain_out);
% Chain_avg=Chain_avg/max(max(Chain_avg));%normalize matrix
% Total_cells(find(Qp==q_p),find(Qm==q_m))=M+P;
    
Total_cell_avg(find(Qp==q_p),find(Qp==q_m))=Total_cell_avg(find(Qp==q_p),find(Qp==q_m))+M+P;

% [i,q_m,q_p,Total_cell_avg(find(Qp==q_p),find(Qm==q_m))]
   
% colormap jet
% imagesc(Total_cells);hold on
% xlabel('q_m')
% ylabel('q_p')
% colorbar()


%TRAVELLING WAVE SPEED
%to calculate the Average Margin velocity
% WaveSpeed2=AveragePropagationSpeed(Chain_avg,T)
% c=c+WaveSpeed2;
WaveSpeed_avg=WaveSpeed_avg+WaveSpeed;
end
% Total_cell_avg(find(Qp==q_p),find(Qm==q_m))=Total_cell_avg(find(Qp==q_p),find(Qm==q_m))/I;
%TRAVELLING WAVE SPEED
Chain_avg=Chain_avg/I;
Density_history_avg=Density_history_avg/I;
WaveSpeed_avg=WaveSpeed_avg/I;
WaveSpeed_vector=[WaveSpeed_vector WaveSpeed_avg];
% [i q_m q_p M+P Total_cell_avg(find(Qp==q_p),find(Qm==q_m))]
% c=c/I;%Average Margin Velocity
% C=[C c]
%Analytical velocity
% c_analytical=AnalyticalSpeed(alpha,v/2,mu,q_p,q_m);
% C_analytical=[C_analytical c_analytical];
Total_cell_avg(find(Qp==q_p),find(Qp==q_m))=Total_cell_avg(find(Qp==q_p),find(Qp==q_m))/I;

% imagesc(Chain_avg)
colormap jet
imagesc(Total_cell_avg)
colorbar()
pause(0.00001)
% plot(Mcells_FinalTime);hold on
% title(strcat('k=',num2str(k)))
% plot(Pcells_FinalTime)
% legend('M cells','P cells')
% hold off
% pause(0.01)
% colormap jet
% imagesc(Total_cells)
% xlabel('q_m')
% ylabel('q_p')
% colorbar()

end
end

% subplot(2,1,1)
% imagesc(Chain_avg)
% subplot(2,1,2)
%DENSITY PLOTS (FOR COMPARISON WITH CONTINUUM)
% for i=1:200
% Compare with Logistic curve
% x0=C_0;
% [t,Continuum_approx_Prolif]=ode45('logistic',[0 T],x0);
% plot(0:T,Density_history_avg,'Linewidth',2)
% hold on
% plot(t,Continuum_approx_Prolif,'Linewidth',2)
% ylabel('Probability')
% xlabel('t')
% legend('Stochastic','Continuum')
% ylim([0 1])
% g=Chain_avg(:,100);
% plot(0:100,g(100:200),'black','Linewidth',2)
% ylim([0 1])
% pause(0.1)
% end
end

% figure
% imagesc(Chain_avg)
% ylabel('i')
% xlabel('j')
% title('Probability of occupancy Q(i, j)')
% colorbar
% colormap(flipud(gray))
% figure
% plot(Chain_avg(100,:),'black','LineWidth',2)
% ylabel('Probability Q(i, j = 100)')
% xlabel('i (position)')

%Density plots
% Chain_avg=Chain_avg(100,:);
% Chain_avg=Chain_avg(100:200);
% plot(0:100,Chain_avg,'Linewidth',2)

% EXPERIMENTAL ANALYSIS ------------
% SINGLE-SPECIES:
% [Chain_out2,Chain0,Chain500,Chain1000,Chain2000,Cells_index_vect,WaveSpeed,Density_history,c2,c,WaveSpeed_vect] = SingleCellSim2(T,alpha,v,pre_Chain,1/500)
% 
% figure
% subplot(2,4,1)
% imagesc(Chain)
% colormap(gca,[1,1,1;0,0,1;1,0,0])
% ylabel('GOG')
% title('0 h')
% subplot(2,4,2)
% imagesc(Chain50)
% colormap(gca,[1,1,1;0,0,1;1,0,0])
% title('4 h')
% subplot(2,4,3)
% imagesc(Chain100)
% colormap(gca,[1,1,1;0,0,1;1,0,0])
% title('24 h')
% subplot(2,4,4)
% imagesc(Chain_out)
% colormap(gca,[1,1,1;0,0,1;1,0,0])
% title('36 h')
% colormap(gca,[1,1,1;0,0,1;1,0,0])
% subplot(2,4,5)
% imagesc(pre_Chain)
% colormap(gca,[1,1,1;0,0,0;0,0,0])
% ylabel('SS')
% subplot(2,4,6)
% imagesc(Chain500)
% colormap(gca,[1,1,1;0,0,0;0,0,0])
% subplot(2,4,7)
% imagesc(Chain1000)
% colormap(gca,[1,1,1;0,0,0;0,0,0])
% subplot(2,4,8)
% imagesc(Chain_out2)
% colormap(gca,[1,1,1;0,0,0;0,0,0])
% Chain_avg=Chain_avg/I;
% %for multiple simulations
% figure
% plot(1:T,Mcells_FinalTime);hold on
% plot(1:T,Pcells_FinalTime)
% hold off
% 
% %Study q_m q_p parameters
figure()
colormap jet
imagesc(Total_cell_avg)
xlabel('q_m')
ylabel('q_p')
colorbar()
figure()
colormap jet
imagesc(linspace(Qm(1),Qm(end),20),linspace(Qp(1),Qp(end),20),Total_cell_avg)
caxis([0 12500])
xlabel('q_m')
ylabel('q_p')
colorbar()



% % %TRAVELLING WAVE SPEED PLOTS
% figure()
% plot(0:30,C_analytical)
% hold on
% plot(0:30,C,'--')
% legend('Analytical','Numerical')
% xlabel('q_p (cell cycle^{-1})')
% ylabel('Wave speed c')

% profile viewer

%      figure()
% imagesc(Q_ij)
% colorbar('Ticks',[0,1,2],...
%          'TickLabels',{'Empty Space','Proliferating','Moving'})
%      
% plot(Chain_avg(50,1:end))

%surf(Chain_avg);hold on
% figure(2)
% plot(P_vector);hold on
% plot(M_vector)
% xlabel('time')
% ylabel('cell number')
% legend('Proliferating','Moving')
% figure(3)
% plot(random_indices_M);hold on
% plot(random_indices_P)
% xlabel('time')
% ylabel('indices when moving/diving')
% legend('Proliferating','Moving')
% figure(4)
% hist(random_indices_P,unique(random_indices_P))
% legend('Proliferating')
% figure(5)
% hist(random_indices_M,unique(random_indices_M))
% legend('moving')
% figure(6)
% hist([random_indices_M random_indices_P],unique([random_indices_M random_indices_P]))
% legend('all displacement')
%Iterate: Go back to Monte Carlo unless the number of reactants is zeroor
%the simulation time has been exceeded.