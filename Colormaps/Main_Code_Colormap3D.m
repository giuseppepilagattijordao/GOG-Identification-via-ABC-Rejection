% profile on
%Gillespie algorithm for Gerlee's IB model
clear all
%Initialize number of molecules in system, reaction constants and random
%number generators
Pcells_FinalTime=[];%make list of final population
Mcells_FinalTime=[];
N1=200;
N2=200;
N3=200;
Chain_avg=zeros(N1,N2);%Chain to get the overall average of simulations
Total_cell_avg=zeros(20,20);%matrix to study q rates with average number of cells
Total_cells=zeros(20,20);%matrix to study q_p and q_m parameters with total cells
WaveSpeed_vector=[];
for k=1
    C=[];
    c=0;
    C_analytical=[];
    I=10; %number of simulations
    Qp=0:1.55:30;
    Qm=0:1.55:30;
%     Qp=20;
%     Qm=10;
for  q_m=Qm
for q_p=Qp
    WaveSpeed_avg=0;
    
for i=1:I
%     [i q_p q_m]
% figure()
% for q_m=1:30
for z=1:N3
    Chain(:,:,z)=zeros(N1,N2);
end

Chain(round(N1/2),round(N2/2),round(N3/2))=1; %single P-cell

M=0;%initial number of M cells (random walk)
% P=0;
P=nnz(Chain);%initial number of P cells (stationary and proliferating)
alpha=1;
v=5;
% q_p=15;%rate from m to p
% q_m=15;%rate from p to m
mu=1e-3;%death rate
T=25;

t=0;
dim=size(Chain);
M_vector=[];
P_vector=[];
random_indices_P=[];
random_indices_M=[];
[M,P,Chain_out]=Gerlee_Main_loop3D(Chain,alpha,v,q_p,q_m,mu,T,t);
Mcells_FinalTime=[Mcells_FinalTime M];
Pcells_FinalTime=[Pcells_FinalTime P];

Chain_avg=Chain_avg+Matrix2Binary(Chain_out);

Total_cell_avg(find(Qp==q_p),find(Qm==q_m))=Total_cell_avg(find(Qp==q_p),find(Qm==q_m))+M+P;
   [i q_p q_m M+P]

end
Total_cell_avg(find(Qp==q_p),find(Qm==q_m))=Total_cell_avg(find(Qp==q_p),find(Qm==q_m))/I;
colormap jet
imagesc(Total_cell_avg);hold on
xlabel('q_m')
ylabel('q_p')
colorbar()
pause(0.0000001)
end
end

end
figure()
colormap jet
imagesc(linspace(Qm(1),Qm(end),20),linspace(Qp(1),Qp(end),20),Total_cell_avg)
caxis([0 70000])
xlabel('q_m')
ylabel('q_p')
colorbar()
