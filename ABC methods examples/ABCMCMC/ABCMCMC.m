%ABC MCMC
clear all
N=10000;
epsilon=0.025;
%Initialize theta_1
theta_1=0;
theta_i=theta_1;
Theta=[theta_i];
i=1;
    while i<N
        %data from target distribution
        y_obs=(1/2)*normrnd(0,1/100)+(1/2)*normrnd(0,1);
            %Generate candidate value 
        theta_hat=normrnd(theta_i,1); % q()
        %Generate data set
        y_hat=normrnd(theta_hat,1,[1 100]);
        %Set next theta equal to the candidate value 
        %with Metropolis-Hastings probability:
        if distance(y_obs,y_hat,theta_hat)<epsilon
            alpha=min(1,(20*rand()-10)*normrnd(theta_hat,1)/((20*rand()-10)*normrnd(theta_i,1))*1);
        else
            alpha=0;
        end
        r=rand();
        if r<alpha %accepted
            theta_iplus1=theta_hat;
            Theta=[Theta theta_iplus1];
        else
            theta_iplus1=theta_i;
        end
        if i<N
            i=i+1;
        end
    end
%     
figure
[f,xi]=ksdensity(Theta);
bar(xi,f)
hold on
%actual distribution
plot(xi,(normpdf(xi,0,sqrt(1))+normpdf(xi,0,sqrt(1/100)))/2,'black','LineWidth',1.5)
xlabel('\theta')
ylabel('Posterior density')     