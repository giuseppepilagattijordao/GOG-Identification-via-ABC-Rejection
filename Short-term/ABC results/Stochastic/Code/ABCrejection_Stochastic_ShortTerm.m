%ABC rejection for Single Cell model
clear all
clc
Nsim=100; %number of trials for comparison 

winner=[];%register winner for all trials
PR_vect=[];
PR_analysis=[];
LR_analysis=[];
LR=[];
AIC=[];
BIC=[];
Count_rejections=[];
%for the calculation of averages
AIC_GOG_avg=0;
AIC_SS_avg=0;
BIC_GOG_avg=0;
BIC_SS_avg=0;
Rel_likelihood_avg=0;
LR_avg=0;
PR_avg=0;
%-----------------------------------------------------------
%GO OR GROW
%lists to store accepted values of each parameter
alpha_vector=[];
v_vector=[];
mu_vector=[];
%create y data
%PARAMETER VALUES FOR GOOD IDENTIFICATION according to number of rejections:
alpha=0.001;
v=100;
mu=0.001;
q_p=1;
q_m=1;

%-------------------------
%Initial condition filled circle
N1=50;
N2=50;
Chain=zeros(N1,N2);
Chain2=Chain;
C=[round(N1/2),round(N1/2)];
R=ones(size(C,1),1)*5;%radius of circles
Chain1=MakeCircles(Chain,N1,N2,C,R);
for i=1:numel(Chain1)
    if Chain1(i)==1
        Chain2(i)=1+round(rand()); %randomly assign motile/proliferative to each position of initial condition
    end
end


t=0;
dim=size(Chain1);
M_vector=[];
P_vector=[];
random_indices_P=[];
random_indices_M=[];
%GENERATE DATA FROM MODELS
Y_obs=[];
tau=1/500;
sim=10;
%to collect all results of all cell cycles
All_SS_Results1=[];
All_SS_Results2=[];
All_SS_Results3=[];
All_GOG_Results1=[];
All_GOG_Results2=[];
All_GOG_Results3=[];

ss_D1_ss=[];%tumor size
ss_D1_gog=[];
ss_D2_ss=[];%Wave speed
ss_D2_gog=[];
ss_D3_ss=[];%Wave slope
ss_D3_gog=[];
gog_D1_ss=[];%tumor size
gog_D1_gog=[];
gog_D2_ss=[];%Wave speed
gog_D2_gog=[];
gog_D3_ss=[];%Wave slope
gog_D3_gog=[];

Chain_avg1=zeros(N1,N2);
Chain_avg2=zeros(N1,N2);
WaveSpeed_avg1=0;
WaveSpeed_avg2=0;
tic
for model=1:2
    
for T=0:0.1:2
Y_mean=0;
display('Generating data, please wait...')
for j=1:sim
%     j
if model==1
    [Chain_out1,Chain0,Chain500,Chain1000,Chain2000,Cells_index_vect,WaveSpeed,Density_history,c2,c,WaveSpeed_vect] = SingleCellSim2(T,alpha,v,mu,Chain1,tau);
    Y_obs1=nnz(Chain_out1);%get tumor size from density, and take average accross all cell cycles

    Chain_avg1=Chain_avg1+Chain_out1;

elseif model==2
    [M,P,Chain_out2,Chain50,Chain100,chain200,Q_ij,WaveSpeed,TumorSize_cellcycles,WaveSpeed_vect]=Gerlee_Main_loop(Chain2,alpha,v,q_p,q_m,mu,T,t);
    Y_obs1=nnz(Chain_out2);%TumorSize_cellcycles is a density obtain by diving by the number of sites in the lattice 

    Chain_avg2=Chain_avg2+Chain_out2;

end
end
Chain_avg1=Chain_avg1/sim;
Chain_avg2=Chain_avg2/sim;
if model==1
    Y_obs2=AveragePropagationSpeed(Chain_avg1,T);
    Y_obs3=max(abs(gradient(Chain_avg1(:,N1/2))));
elseif model==2
    Y_obs2=AveragePropagationSpeed(Chain_avg2,T);
    Y_obs3=max(abs(gradient(Chain_avg2(:,N2/2))));
end

display('Data generated.')

Y_sim=[];
Y_sim2=[];
y_sim=0;
Count_avg=0;
ALPHA_SS=[];
NU_SS=[];
ALPHA_GOG=[];
NU_GOG=[];
%alternative way of comparison: counting number of rejections
count_SS_rejections=0;
count_GOG_rejections=0;
SS_win1=0;
SS_win2=0;
SS_win3=0;
GOG_win1=0;
GOG_win2=0;
GOG_win3=0;
draws1=0;
draws2=0;
draws3=0;
N=100;
D1_ss=[];
D2_ss=[];
D3_ss=[];
D1_gog=[];
D2_gog=[];
D3_gog=[];
for trial=1:N
    ['Time: ' num2str(T) ', ' num2str(trial)]
    ['Tumor Size: ' num2str(SS_win1) ', ' num2str(GOG_win1) ', ' num2str(draws1) newline 'Wave Speed: ' num2str(SS_win2) ', ' num2str(GOG_win2) ', ' num2str(draws2) newline 'Wave Slope: ' num2str(SS_win3) ', ' num2str(GOG_win3) ', ' num2str(draws3)]
        %Short-term
        alpha_hat1=0.002*rand(); %sample from prior for SS model
        v_hat1=100*rand()+50;
        alpha_hat2=0.002*rand(); %sample from prior for GoG model
        v_hat2=100*rand()+50;

        %parallel computations:
        y_sim1=cell(1,1); %sliced variable
        y_sim2=cell(1,1);
        y_sim3=cell(1,1);
        y_sim2_1=cell(1,1); %sliced variable
        y_sim2_2=cell(1,1);
        y_sim2_3=cell(1,1);
        parfor i = 1:2
            Chain_out_avg1=zeros(N1,N2);
            Chain_out_avg2=zeros(N1,N2);
            if i == 1
                for k=1:sim
                      [Chain_out1,Chain0,Chain500,Chain1000,Chain2000,Cells_index_vect,WaveSpeed,Density_history,c2,c,WaveSpeed_vect] = SingleCellSim2(T,alpha_hat1,v_hat1,mu,Chain1,tau)
                      y_sim1{i}=nnz(Chain_out1);%tumor size in number of cells
%                       y_sim2{i}=WaveSpeed;%wave speed
                        Chain_out_avg1=Chain_out_avg1+Chain_out1;
%                       y_sim3{i}=max(gradient(Chain_out1(:,25)));%max slope of wave
%                         y_sim3{i}=mean(mean(gradient(Matrix2Binary(Chain_out1))));
                end
                Chain_out_avg1=Chain_out_avg1/sim;
                y_sim2{i}=AveragePropagationSpeed(Chain_out_avg1,T);
                y_sim3{i}=max(abs(gradient(Chain_out_avg1(:,50))));
            else
                for k=1:sim
                      [M,P,Chain_out2,Chain50,Chain100,chain200,Q_ij,WaveSpeed,TumorSize_cellcycles,WaveSpeed_vect]=Gerlee_Main_loop(Chain2,alpha_hat2,v_hat2,q_p,q_m,mu,T,t);
                      y_sim2_1{i}=nnz(Chain_out2);%probability of occupancy
%                       y_sim2_2{i}=WaveSpeed;%wave speed
                        Chain_out_avg2=Chain_out_avg2+Matrix2Binary(Chain_out2);

                end
                Chain_out_avg2=Chain_out_avg2/sim;
                y_sim2_2{i}=AveragePropagationSpeed(Chain_out_avg2,T);
                y_sim2_3{i}=max(abs(gradient(Chain_out_avg2(:,50))));
            end
        end
        y_sim1=cell2mat(y_sim1);%tumor mass ss
        y_sim2=cell2mat(y_sim2);%wave speed ss
        y_sim3=cell2mat(y_sim3);%wave slope ss
        y_sim2_1=cell2mat(y_sim2_1);%tumor mass gog
        y_sim2_2=cell2mat(y_sim2_2);%wave speed gog
        y_sim2_3=cell2mat(y_sim2_3);%wave slope gog
%         [T Y_obs y_sim y_sim2]
        d1_ss=abs(y_sim1-Y_obs1);%SS model
        d2_ss=abs(y_sim2-Y_obs2);
        d3_ss=abs(y_sim3-Y_obs3);
        d1_gog=abs(y_sim2_1-Y_obs1);%GoG model
        d2_gog=abs(y_sim2_2-Y_obs2);
        d3_gog=abs(y_sim2_3-Y_obs3);
        [T Y_obs d1_ss d1_gog d2_ss d2_gog d3_ss d3_gog]
        D1_ss=[D1_ss d1_ss];
        D2_ss=[D2_ss d2_ss];
        D3_ss=[D3_ss d3_ss];
        D1_gog=[D1_gog d1_gog];
        D2_gog=[D2_gog d2_gog];
        D3_gog=[D3_gog d3_gog];
        count=0;%to count rejection rate
        %WINNER FOR TUMOR SIZE
        if norm(d1_ss,Inf)<norm(d1_gog,Inf) %SS wins in tumor size
            SS_win1=SS_win1+1;
        elseif norm(d1_ss,Inf)>norm(d1_gog,Inf)%GOG wins in tumor size
            GOG_win1=GOG_win1+1;
        else
            draws1=draws1+1;
        end
        %WINNER FOR WAVE SPEED
        if norm(d2_ss,Inf)<norm(d2_gog,Inf) %SS wins in wave speed
            SS_win2=SS_win2+1;
        elseif norm(d2_ss,Inf)>norm(d2_gog,Inf)%GOG wins in wave speed
            GOG_win2=GOG_win2+1;
        else
            draws2=draws2+1;
        end
        %WINNER FOR WAVE SLOPE
        if norm(d3_ss,Inf)<norm(d3_gog,Inf) %SS wins in wave slope
            SS_win3=SS_win3+1;
        elseif norm(d3_ss,Inf)>norm(d3_gog,Inf)%GOG wins in wave slope
            GOG_win3=GOG_win3+1;
        else
            draws3=draws3+1;
        end
        
end
%Average of distances for this Nr of cell cycles
%specifically
D1_ss=mean(D1_ss);
D2_ss=mean(D2_ss);
D3_ss=mean(D3_ss);
D1_gog=mean(D1_gog);
D2_gog=mean(D2_gog);
D3_gog=mean(D3_gog);

if model==1
    ss_D1_ss=[ss_D1_ss D1_ss];%tumor size
    ss_D1_gog=[ss_D1_gog D1_gog];
    ss_D2_ss=[ss_D2_ss D2_ss];%Wave speed
    ss_D2_gog=[ss_D2_gog D2_gog];
    ss_D3_ss=[ss_D3_ss D3_ss];%Wave slope
    ss_D3_gog=[ss_D3_gog D3_gog];
else
    gog_D1_ss=[gog_D1_ss D1_ss];
    gog_D1_gog=[gog_D1_gog D1_gog];
    gog_D2_ss=[gog_D2_ss D2_ss];
    gog_D2_gog=[gog_D2_gog D2_gog];
    gog_D3_ss=[gog_D3_ss D3_ss];
    gog_D3_gog=[gog_D3_gog D3_gog];
end

%Alternative evaluation
if model==1
    All_SS_Results1=[All_SS_Results1 SS_win1];%tumor size
    All_SS_Results2=[All_SS_Results2 SS_win2];%Wave speed
    All_SS_Results3=[All_SS_Results3 SS_win3];%Wvae slope
else
    All_GOG_Results1=[All_GOG_Results1 GOG_win1];%tumor size
    All_GOG_Results2=[All_GOG_Results2 GOG_win2];%Wave speed
    All_GOG_Results3=[All_GOG_Results3 GOG_win3];%Wave slope
end
end

end
toc