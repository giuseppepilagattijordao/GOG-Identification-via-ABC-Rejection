%Test Analytical Speed 
C_analytical=[];
alpha=1;
v=5;
mu=0.001;
q_p=22;
%TO COMPARE WITH CONTINUOUS MODEL
display('First')
%q_m
for q_m=0:30
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];
end
figure(1)
plot(0:30,C_analytical,'black')
xlabel('q_m (cell cycle^{-1})')
ylabel('wave speed c')
%q_p
display('Second')
C_analytical=[];
q_m=6.3;
for q_p=0:30
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];
end
figure(2)
plot(0:30,C_analytical,'black')
xlabel('q_p (cell cycle^{-1})')
ylabel('wave speed c')
%TO COMPARE WITH STOCHASTIC MODEL
%q_m
display('Third')
C_analytical=[];
v=5/2;%rescale motility rate to take into account the 2-Dimensions of motion
q_p=15;
for q_m=0:30
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];
end
figure(3)
plot(0:30,C_analytical,'black')
xlabel('q_m (cell cycle^{-1})')
ylabel('wave speed c')
%q_p
display('Fourth')
C_analytical=[];
q_m=15;
for q_p=0:30
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];
end
figure(4)
plot(0:30,C_analytical,'black')
xlabel('q_p (cell cycle^{-1})')
ylabel('wave speed c')
%OTHER PARAMETERS
%alpha
display('Fifth')
C_analytical=[];
v=5;
q_p=20;
q_m=10;
for alpha=0.0001:0.1:3
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];
end
figure(5)
plot(0.0001:0.1:3,C_analytical,'black')
xlabel('\alpha (proliferation rate)')
ylabel('wave speed c')

%nu
C_analytical=[];
display('sixth')
alpha=1;
for v=0.0001:0.1:100
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];
end
figure(6)
plot(0.0001:0.1:100,C_analytical,'black')
xlabel('\nu (motility rate)')
ylabel('wave speed c')
%mu
display('seventh')
C_analytical=[];
v=5;
for mu=0:0.01:1
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];
end
k=0:0.01:1;
figure(7)
plot(k,C_analytical,'black')
xlabel('\mu (apoptosis rate)')
ylabel('wave speed c')
%for mu=0:0.01:1, find(abs(C_analytical-0.1349)<0.0001)=68
%and k=0:0.01:1 means k(68)=0.67