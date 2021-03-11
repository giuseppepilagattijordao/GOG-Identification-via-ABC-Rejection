function [c,C_analytical,p,m] = GerleeContinuumFunction(x,dx,dt,p,m,CellCycles,alpha,v,mu,q_p,q_m)
%GERLEE CONTINUUM MODEL

%NON DIMNESIONALIZE (TIME AND SPACE)
C=[]; %for all wave speeds
c=0;
C_analytical=[];
Xf=[];
Xi=[];
beta=10;

if CellCycles==0
    p=p';
    m=m';
    return
end

%CALCULATE ANALYTICAL TRAVELLING WAVE SPEEDS
c_analytical=AnalyticalSpeed(alpha,v,mu,q_p,q_m);
C_analytical=[C_analytical c_analytical];

%NUMERICAL PARAMETERS
h=1;
D_alpha=(h^2)*alpha/2;
D_v=(h^2)*v/2;

p0=p;
m0=m;
epsilon=1e-8;

if size(p,2)>1
    p=p';
    m=m';
end
% m=zeros(length(x),1);
p_previous=p;
m_previous=m;
j40_check=1;%for the snapshots to be at the right time point
j50_check=1;
for j=0:dt:CellCycles %time iterations

    m_u_up=0;
    dp=0.0001*ones(length(x),1);
    dm=0.0001*ones(length(x),1);

    while norm(dp,inf)>epsilon || norm(dm,inf)>epsilon
       %---------------------------
       %Start Boundary conditions for i=1
        P(1)=p(2)-p(1);
        M(1)=m(2)-m(1);
        d1_p(1)=0;%dF/du_iminus1=0
        d2_p(1)=-1;%dF/du_i=-1
        d3_p(1)=1;%dF/du_iplus1=1
        d1_m(1)=0;%dF/du_iminus1=0
        d2_m(1)=-1;%dF/du_i=-1
        d3_m(1)=1;%dF/du_iplus1=1

        for i=2:length(x)-1
            
            
            %Crank-Nicolson (Average between implicit and explicit
%             IMPLICIT
            d2pdx2=(p(i-1)+p(i+1)-2*p(i))/dx^2;
            d2mdx2=(m(i-1)+m(i+1)-2*m(i))/dx^2;
            %Algebraic equation
            P(i)=-(p(i)-p_previous(i))/dt+(D_alpha)*(1-p(i)-m(i))*(d2pdx2)+alpha*p(i)*(1-p(i)-m(i))-(q_m+mu)*p(i)+q_p*m(i);
            M(i)=-(m(i)-m_previous(i))/dt+(D_v)*((1-p(i))*(d2mdx2)+m(i)*(d2pdx2))-(q_p+mu)*m(i)+q_m*p(i);
            if isnan(P(i))||isnan(M(i))
                display('stop')
                %return
            end
                  
            d1_p(i)=(D_alpha/(dx^2))*(1-p(i)-m(i));
            d2_p(i)=-1/dt+(D_alpha/(dx^2))*(-p(i-1)+4*p(i)-p(i+1)+2*m(i)-2)+alpha*(1-2*p(i)-m(i))-(q_m+mu);
            d3_p(i)=(D_alpha/(dx^2))*(1-p(i)-m(i));
            d1_m(i)=(D_v/(dx^2))*(1-p(i));
            d2_m(i)=-1/dt+(D_v/(dx^2))*(-2+p(i-1)+p(i+1))-(q_p+mu);
            d3_m(i)=(D_v/(dx^2))*(1-p(i));
                      
        end
        % End Boundary Conditions

        d1_p(length(x))=0;
        d2_p(length(x))=-1;
        d3_p(length(x))=1;
        d1_m(length(x))=0;
        d2_m(length(x))=-1;
        d3_m(length(x))=1;
        P(length(x))=p(length(x))-p(length(x)-1);
        M(length(x))=m(length(x))-m(length(x)-1);

        %-----------------------------------------
        d1_p=[0 d1_p];
        d3_p=[d3_p 0];
        d1_m=[0 d1_m];
        d3_m=[d3_m 0];
        dp=tridiag(d1_p',d2_p',d3_p',-P);
        dm=tridiag(d1_m',d2_m',d3_m',-M);

        if isnan(dp(1))||isnan(dm(1))
            display('stop')
        end
        %-------------------------------------------
        p=p+dp;
        m=m+dm;
        
        
    end   
      %keep record of previous u
     p_previous=p;
     m_previous=m;
%%%---------------------------------------------
end
%WAVE SPEED---------------------------------
p_i=p0;
m_i=m0;
pm_i=p_i+m_i; %sum of m and p

p_f=p;
m_f=m;
pm_f=p_f+m_f; %sum of m and p
%get point close to 1/2
x_i=1800;
x_f=min(find(abs(pm_f-max(pm_f)/2)<=(0.01)));
c=abs((dx*x_f-dx*x_i))/CellCycles;%wave speed
if isempty(c)==1
   c=0; 
end
end
