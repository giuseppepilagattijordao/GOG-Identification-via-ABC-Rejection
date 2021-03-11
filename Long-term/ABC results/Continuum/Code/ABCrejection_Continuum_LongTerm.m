
%ABC rejection for comparing Continuum Approximations
clear all
clc

tic
winner=[];%register winner for all trials
LR=[];
AIC=[];
BIC=[];
%-----------------------------------------------------------
%GO OR GROW
%lists to store accepted values of each parameter
alpha_vector=[];
v_vector=[];
mu_vector=[];
%parameter values
v=1; 
alpha=0.01;
mu=0.001;
q_p=1;
q_m=1;
%INITIAL CONDITIONS
dt=0.5;
dx=0.01;
x_inf=50;%25 for short term (50x50 lattice) / 50 for long-term (100x100 lattice)
x=0:dx:x_inf;
%SS
%Short-term
% a=exp(-1*(x));
%Long-term
a=exp(-0.5*(x));
%GOG
%Short-term
% m=exp(-1*(x))*0.5;
% p=exp(-1*(x))*0.5;
%Long-term
m=exp(-0.5*(x))*0.5;
p=exp(-0.5*(x))*0.5;
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

all_D1_SS=[];
all_D2_SS=[];
all_D3_SS=[];
all_D1_GOG=[];
all_D2_GOG=[];
all_D3_GOG=[];
tic
for model=1:2
for CellCycles=0:2:30
display('Generating data, please wait...')
if model==1
    [C,a_out]=SingleSpeciesFunction(x,dx,dt,a,CellCycles,alpha,v,mu);
    Y_obs1=a_out;%probability of occupancy 
    Y_obs2=C; %WaveSpeed
    Y_obs3=max(abs(gradient(a_out)));%max slope of wave 
elseif model==2
    [C,C_analytical,p_out,m_out] = GerleeContinuumFunction(x,dx,dt,p,m,CellCycles,alpha,v,mu,q_p,q_m);
    Y_obs1=p_out'+m_out';
    Y_obs2=C; %Wave Speed
    Y_obs3=max(abs(gradient(p_out'+m_out')));
end
display('Data generated.')
epsilon=10;
Y_sim=[];
Y_sim2=[];
y_sim=0;
Count_avg=0;
ALPHA_SS=[];
NU_SS=[];
ALPHA_GOG=[];
NU_GOG=[];
N=100;
SS_win1=0;
SS_win2=0;
SS_win3=0;
GOG_win1=0;
GOG_win2=0;
GOG_win3=0;
draws1=0;
draws2=0;
draws3=0;

draws=0;
D1_ss=[];
D2_ss=[];
D3_ss=[];
D1_gog=[];
D2_gog=[];
D3_gog=[];
for trial=1:N
    ['Time: ' num2str(CellCycles) ', ' num2str(trial)]
    ['Tumor Size: ' num2str(SS_win1) ', ' num2str(GOG_win1) ', ' num2str(draws1) newline 'Wave Speed: ' num2str(SS_win2) ', ' num2str(GOG_win2) ', ' num2str(draws2) newline 'Wave Slope: ' num2str(SS_win3) ', ' num2str(GOG_win3) ', ' num2str(draws3)]
        %Long-term
        alpha_hat1=0.02*rand()+7e-4; %sample from prior for SS model
        v_hat1=2*rand();
        alpha_hat2=0.02*rand()+7e-4; %sample from prior for GoG model
        v_hat2=2*rand();
        %parallel computations:
        y_sim1=cell(1,1); %sliced variable
        y_sim2=cell(1,1);
        y_sim3=cell(1,1);
        y_sim2_1=cell(1,1); %sliced variable
        y_sim2_2=cell(1,1);
        y_sim2_3=cell(1,1);
        parfor i = 1:2
            if i == 1
              [C,a_out]=SingleSpeciesFunction(x,dx,dt,a,CellCycles,alpha_hat1,v_hat1,mu);
              y_sim1{i}=a_out;%probability of occupancy
              y_sim2{i}=C;%wave speed
              y_sim3{i}=max(abs(gradient(a_out)));%max slope of wave 
            else
              [C,C_analytical,p_out,m_out] = GerleeContinuumFunction(x,dx,dt,p,m,CellCycles,alpha_hat2,v_hat2,mu,q_p,q_m);
              y_sim2_1{i}=p_out'+m_out';%probability of occupancy
              y_sim2_2{i}=C;%wave speed
               y_sim2_3{i}=max(abs(gradient(p_out'+m_out')));%max slope of wave 
            end
        end
        y_sim1=cell2mat(y_sim1);
        y_sim2=cell2mat(y_sim2);
        y_sim3=cell2mat(y_sim3);
        y_sim2_1=cell2mat(y_sim2_1);
        y_sim2_2=cell2mat(y_sim2_2);
        y_sim2_3=cell2mat(y_sim2_3);
        d1_ss=abs(y_sim1-Y_obs1);%SS model
        d2_ss=abs(y_sim2-Y_obs2);%SS model
        d3_ss=abs(y_sim3-Y_obs3);%SS model
        d1_gog=abs(y_sim2_1-Y_obs1);%GoG model
        d2_gog=abs(y_sim2_2-Y_obs2);%GoG model
        d3_gog=abs(y_sim2_3-Y_obs3);%GoG model
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