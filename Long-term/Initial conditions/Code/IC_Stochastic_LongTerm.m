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

%PARAMETER
alpha=0.01;
v=1;
mu=0.001;
q_p=1;
q_m=1;
%INITIAL CONDITIONS FOR 20X400------------------
% N1=20;
% N2=400;
% %for SS
%     Chain1=zeros(N1,N2);
%     for i=1:size(Chain1,1)
%         for j=180:220
%             Chain1(i,j)=1;
%         end
%     end
% %for GOG
%     Chain2=zeros(N1,N2);
%     r=rand();
%     for o=1:size(Chain2,1)
%         for j=180:220
%             if r<0.5
%                 Chain2(o,j)=1;
%             else
%                 Chain2(o,j)=2;
%             end
%         end
%     end
%-------------------------
%Initial condition filled circle
N1=100;
N2=100;
Chain=zeros(N1,N2);
Chain2=Chain;
C=[round(N1/2),round(N1/2)];
R=ones(size(C,1),1)*10;%radius of circles
Chain1=MakeCircles(Chain,N1,N2,C,R);
for i=1:numel(Chain1)
    if Chain1(i)==1
        Chain2(i)=1+round(rand()); %randomly assign motile/proliferative to each position of initial condition
    end
end

figure
imagesc(Chain1)
colormap(gca,[1,1,1;0,0,0;0,0,0])
title('SS')
xlabel('j')
ylabel('i')
figure
imagesc(Chain2)
colormap(gca,[1,1,1;0,0,1;1,0,0])
title('GOG')
xlabel('j')
ylabel('i')
