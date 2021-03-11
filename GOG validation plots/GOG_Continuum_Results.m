alpha=1;
v=5;
mu=0.001;
T=[40 50];
beta=10;
%__________________________________________
%Figure with q_p > q_m
%numerical parameters
dx=0.25;
x=0:dx:100;
dt=0.19;
%Initial conditions
p0=exp(-beta*x);
m0=zeros*x;
%transition parameters
q_m=10;
q_p=20;

figure
for t=T
    %define color string
    if t==40
        cl='green'
    elseif t==50
        cl='blue'
    end
    %run model
    [c,C_analytical,p_out,m_out]=GerleeContinuumFunction(x,dx,dt,p0,m0,t,alpha,v,mu,q_p,q_m);
    %plots
    plot(x,p_out,cl)
    hold on
    plot(x,m_out,[cl '--'])
    xlim([0 100])
    ylim([0 1])
    xlabel('x')
    ylabel('Probability of occupancy')
end
%__________________________________________
%Figure with q_p < q_m
%numerical parameters
dx=0.1;
x=0:dx:100;
dt=0.1;
%Initial conditions
p0=exp(-beta*x);
m0=zeros*x;
%transition parameters
q_m=20;
q_p=10;

figure
for t=T
    %define color string
    if t==40
        cl='green'
    elseif t==50
        cl='blue'
    end
    %run model
    [c,C_analytical,p_out2,m_out2]=GerleeContinuumFunction(x,dx,dt,p0,m0,t,alpha,v,mu,q_p,q_m);
    %plots
    plot(x,p_out2,cl)
    hold on
    plot(x,m_out2,[cl '--'])
    xlim([0 100])
    ylim([0 1])
    xlabel('x')
    ylabel('Probability of occupancy')
end

%Figure with q_p > q_m
%numerical parameters
dx=0.25;
x=0:dx:100;
dt=0.09;
%Initial conditions
p0=exp(-beta*x);
m0=zeros*x;
%transition parameters
q_m=0;
q_p=10;

figure
for t=T
    %define color string
    if t==40
        cl='green'
    elseif t==50
        cl='blue'
    end
    %run model
    [c,C_analytical,p_out3,m_out3]=GerleeContinuumFunction(x,dx,dt,p0,m0,t,alpha,v,mu,q_p,q_m);
    %plots
    plot(x,p_out3,cl)
    hold on
    plot(x,m_out3,[cl '--'])
    xlim([0 100])
    ylim([0 1])
    xlabel('x')
    ylabel('Probability of occupancy')
end