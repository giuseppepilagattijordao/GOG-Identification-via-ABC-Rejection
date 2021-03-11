%ABC rejection
clear all
N=1000;
epsilon=0.025;
%Initialize theta_1
theta_1=0;
theta_i=theta_1;
Theta=[theta_i];
for i=1:N
    i
    y_obs=(1/2)*normrnd(0,1/100)+(1/2)*normrnd(0,1);
    theta_hat=normrnd(theta_i,1); 
     y_hat=normrnd(theta_hat,1,[1 100]);
    while distance(y_obs,y_hat,theta_hat)>=epsilon
        %data from target distribution
        y_obs=(1/2)*normrnd(0,1/100)+(1/2)*normrnd(0,1);
            %Generate candidate value 
        theta_hat=normrnd(theta_i,1); % q()
        %Generate data set
        y_hat=normrnd(theta_hat,1,[1 100]);
        %Set next theta equal to the candidate value 
    end  
        theta_iplus1=theta_hat;%accept
        Theta=[Theta theta_iplus1];
    end
%     
figure
[f,xi]=ksdensity(Theta);
bar(xi,f)
hold on
%actual distribution
plot(xi,(normpdf(xi,0,sqrt(1))+normpdf(xi,0,sqrt(1/100)))/2,'black','LineWidth',1.5)
xlabel('\theta')
ylabel('Posterior distribution')