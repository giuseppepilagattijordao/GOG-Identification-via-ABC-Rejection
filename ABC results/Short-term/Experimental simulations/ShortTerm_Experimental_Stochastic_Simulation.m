% profile on
%Gillespie algorithm for Gerlee's IB model
clear all
clc
%Initialize number of molecules in system, reaction constants and random
%number generators
Pcells_FinalTime=[];%make list of final population
Mcells_FinalTime=[];
N1=50;
N2=50;
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
    I=1; %number of simulations
    T=2;
%     Qp=0:1.55:30;
%     Qm=0:1.55:30;
    Qp=1;
    Qm=1;
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
% Chain(round(N1/2),round(N2/2))=1; %single P-cell
% C=[10 10;110 30; 50 90]
C=[round(N1/2),round(N1/2)];
R=ones(size(C,1),1)*5;%radius of circles
pre_Chain=MakeCircles(Chain,N1,N2,C,R);
for i=1:numel(pre_Chain)
    if pre_Chain(i)~=0
%           Chain(i)=2;
        Chain(i)=1+round(rand()); %randomly assign motile/proliferative to each position of initial condition
    end
end

M=0;%initial number of M cells (random walk)
% P=0;
P=1;%initial number of P cells (stationary and proliferating)
alpha=0.001;
v=100;

mu=1e-3;%death rate
t=0;
% T=20;
dim=size(Chain);
M_vector=[];
P_vector=[];
random_indices_P=[];
random_indices_M=[];
[M,P,Chain_out,Chain50,Chain100,Q_ij,WaveSpeed,TumorSize_cellcycles]=Gerlee_Main_loop(Chain,alpha,v,q_p,q_m,mu,T,t);
Mcells_FinalTime=[Mcells_FinalTime M];
Pcells_FinalTime=[Pcells_FinalTime P];
Chain_avg=Chain_avg+Matrix2Binary(Chain_out);
    
Total_cell_avg(find(Qp==q_p),find(Qp==q_m))=Total_cell_avg(find(Qp==q_p),find(Qp==q_m))+M+P;

% [i,q_m,q_p,Total_cell_avg(find(Qp==q_p),find(Qm==q_m))]
   
WaveSpeed_avg=WaveSpeed_avg+WaveSpeed;
end

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


end
end
end


% EXPERIMENTAL ANALYSIS ------------
% SINGLE-SPECIES:
[Chain_out2,Chain0,Chain500,Chain1000,Chain2000,Cells_index_vect,WaveSpeed,Density_history,c2,c,WaveSpeed_vect] = SingleCellSim2(T,alpha,v,mu,pre_Chain,1/500)

figure
subplot(2,4,1)
imagesc(Chain)
colormap(gca,[1,1,1;0,0,1;1,0,0])
ylabel('GOG')
title('0 h')
subplot(2,4,2)
imagesc(Chain50)
colormap(gca,[1,1,1;0,0,1;1,0,0])
title('4 h')
subplot(2,4,3)
imagesc(Chain100)
colormap(gca,[1,1,1;0,0,1;1,0,0])
title('24 h')
subplot(2,4,4)
imagesc(Chain_out)
colormap(gca,[1,1,1;0,0,1;1,0,0])
title('36 h')
colormap(gca,[1,1,1;0,0,1;1,0,0])
subplot(2,4,5)
imagesc(pre_Chain)
colormap(gca,[1,1,1;0,0,0;0,0,0])
ylabel('SS')
subplot(2,4,6)
imagesc(Chain500)
colormap(gca,[1,1,1;0,0,0;0,0,0])
subplot(2,4,7)
imagesc(Chain1000)
colormap(gca,[1,1,1;0,0,0;0,0,0])
subplot(2,4,8)
imagesc(Chain_out2)
colormap(gca,[1,1,1;0,0,0;0,0,0])
