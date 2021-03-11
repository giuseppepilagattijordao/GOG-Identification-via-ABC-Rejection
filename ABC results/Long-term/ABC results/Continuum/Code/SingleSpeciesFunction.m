function [c,a]=SingleSpeciesFunction(x,dx,dt,a_previous,CellCycles,alpha,v,mu)
%STOCHASTIC VS CONTINUUM - SINGLE-SPECIES

%CONTINUUM MODEL

%NON DIMNESIONALIZE (TIME AND SPACE)
C=[]; %for all wave speeds
C_analytical=[];
C_analytical2=[];
Xf=[];
Xi=[];
c=0;
if CellCycles==0
    a=a_previous;
    return
end

%PARAMETERS
beta=10;
%Initial condiitons for stochastic
N1=20;
N2=400;
Chain=zeros(N1,N2);
% Chain(round(N1/2),round(N2/2))=1;
C_0=0.6;

InitialField=Chain;
InitialField=Chain(:,181:220);
while nnz(InitialField)/numel(InitialField)<C_0
    Site=randi(numel(InitialField));
    InitialField(Site)=1;%get position in the area of initial cells
    start=sub2ind(size(Chain),1,180); %row should be 1
    Chain(Site+start)=1;%populate original matrix properly
end
%CALCULATE ANALYTICAL TRAVELLING WAVE SPEEDS
c_analytical=sqrt(2*(alpha-mu)*(alpha+v));%from point (0,0)
C_analytical=[C_analytical c_analytical];
%NUMERICAL PARAMETERS
h=1;
D_alpha=(h^2)*alpha/4;
D_v=(h^2)*v/4;
%Continuum initial conditions
D=1/4;
x0=40;
a=a_previous;
epsilon=1e-6;
a0=a;
a_previous=a;
A=zeros(1,length(x));
d1_a=zeros(1,length(x)-1);
d2_a=zeros(1,length(x));
d3_a=zeros(1,length(x)-1);
for j=0:dt:CellCycles %time iterations
    if j-dt<1000 && j+dt>1000 && abs(floor(j-dt)-floor(j))==1
        a1000=a;
    end
    if j-dt<2000 && j+dt>2000 && abs(floor(j-dt)-floor(j))==1
        a2000=a;
    end
    if j==10/dt+1 %evaluate wave's speed using two snap shots
        a_i=a;
    elseif j==20/dt+1
        a_f=a;
        %get point close to 1/2
        x_i=min(find(abs(a_i-max(a_i)/2)<=(0.005)));
        x_f=min(find(abs(a_f-max(a_f)/2)<=(0.005)));
        c=dx*(x_f-x_i)/(20-10);%wave speed
        Xf=[Xf x_f];
        Xi=[Xi x_i];
        C=[C c];
    end

%     j
    da=0.0001*ones(length(x),1);

    while norm(da,Inf)>epsilon
%            norm(da)
        
       %---------------------------
       %Start Boundary conditions for i=1
        A(1)=(a(2)-a(1));
        d1_a(1)=0;%dF/da_iminus1=0
        d2_a(1)=-1;%dF/da_i=-1
        d3_a(1)=1;%dF/da_iplus1=1
        
        for i=2:length(x)-1
            
            %define algebraic nonlinear equation to solve
            A(i)=-(a(i)-a_previous(i))/dt+(D_v/dx^2)*(a(i-1)-2*a(i)+a(i+1))+alpha*a(i)*(1-a(i))-mu*a(i);
            if isnan(A(1))
                display('stop')
                return
            end
           
            d1_a(i)=(D_v/(dx^2));
            d2_a(i)=-1/dt-2*((D_v/dx^2))+alpha*(1-a(i))-alpha*a(i)-mu;
            d3_a(i)=(D_v/(dx^2));
            
        end
        
        % End Boundary Conditions
        d1_a(end)=0;
        d2_a(end)=-1;
        d3_a(end)=1;
        A(end)=(a(end)-a(end-1));
        %-----------------------------------------
        d1_a=[0 d1_a];
        d3_a=[d3_a 0];
        da=tridiag(d1_a',d2_a',d3_a',-A);
%        dm=Thomas(d1_a',d2_a',d3_a',-A);
        if isnan(da(1))
            display('stop')
            return
        end
        %-------------------------------------------
        a=a+da';       
    end   
      %keep record of previous u
     a_previous=a;
end
%Wave Speed
a_i=a0;
a_f=a;
%get point close to 1/2
x_i=1800;
x_f=min(find(abs(a_f-max(a_f)/2)<=(0.01)));
c=abs(dx*x_f-dx*x_i)/CellCycles;%wave speed
if isempty(c)==1
   c=0; 
end
end