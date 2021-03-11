
 function [M,P,Chain] = Gerlee_Main_loop3D(Chain,alpha,v,q_p,q_m,mu,T,t)
% profile on
[P_index_vect, M_index_vect]=get_index_vectors(Chain);%getting indexes of each species in the matrix
% M_vector=[];
% P_vector=[];
nrows=size(Chain,1);
VonNeumann=[-nrows*nrows,-nrows,-1,1,nrows,nrows*nrows];
P=sum(Chain(:)==1);%count number of P cells
M=sum(Chain(:)==2);%count number of M cells
while t<T
%     t
    %Initial check
    if M+P==0
        break
    end
    %CALCULATE PROPENSITIES
    a1=alpha*P; %Proliferation occurrences
    a2=v*M; %motility occurrences
    a3=q_m*P; %P-M transition occurrences
    a4=q_p*M; %M-P transition occurrences
    a5=mu*P; %death of P cell occurrences
    a6=mu*M; %death of M cell occurrences
    
    a=[a1,a2,a3,a4,a5,a6];
%    a=sort(a);
    %----------------------------------------------------------
    %RANDOM NUMBERS
    r1=rand();
    r2=rand();
    %TIME INCREMENT
    a0=sum(a);
    dt=(1/a0)*log(1/r1);
    %CHOOSE (POSSIBLE) EVENT TO PERFORM  
    %if r2 falls between 0 and a1/a0, choose a1, if it falls between (a1+a2)/a0 choose a2, etc...
    if r2>0 && r2<=a1/a0
        pos=1;
    elseif r2>a1/a0 && r2<=(a1+a2)/a0
        pos=2;
    elseif r2>(a1+a2)/a0 && r2<=(a1+a2+a3)/a0
        pos=3;
    elseif r2>(a1+a2+a3)/a0 && r2<=(a1+a2+a3+a4)/a0
        pos=4;
    elseif r2>(a1+a2+a3+a4)/a0 && r2<=(a1+a2+a3+a4+a5)/a0
        pos=5;
    elseif r2>(a1+a2+a3+a4+a5)/a0 && r2<=1
        pos=6;
    end
%     pos=GillespieReactionChoice(a,r2);
%     pos = min(find(r2*a0 <= cumsum(a)));
    %UPDATE TIME
    t_previous=t;%use this time to count integer times for Q(i,i,t) function
    t=t+dt;
    %-----------------------------------------------------------
%     %PICK RANDOM ADEQUATE CELL FROM MATRIX TO PERFORM CHOSEN EVENT
%     if pos==1 || pos==3 || pos==5
% %         Site=P_index_vect(randperm(size(P_index_vect,1),1),:);
%         Site=P_index_vect(randi(numel(P_index_vect)));
%     elseif pos==2 || pos==4 || pos==6
% %         Site=M_index_vect(randperm(size(M_index_vect,1),1),:);
%         Site=M_index_vect(randi(numel(M_index_vect)));
%     end
    
    %PICK RANDOM ADEQUATE CELL FROM MATRIX TO PERFORM CHOSEN EVENT
    if pos==1 || pos==3 || pos==5
%         Site=P_index_vect(randperm(size(P_index_vect,1),1),:);
        Site=P_index_vect(randi(numel(P_index_vect)));
            if pos==1
                %Checking neighborhood
                VacantSites=zeros(1,length(VonNeumann));%preallocating for speed, instead of setting to []
                for i=1:length(VonNeumann)
                    ng=VonNeumann(i);
                    if Site+ng<=numel(Chain) && Site+ng>=1 && Chain(Site+ng)==0
                        VacantSites(i)=ng;
                    end
                end
                VacantSites=VacantSites(VacantSites~=0);%exclude any zero in VacantSites
                if isempty(VacantSites)==1
                    randomSite=0; %division fails = do nothing
                    pos=7; %do nothing event
                else
                    randomSite=VacantSites(randi(numel(VacantSites))); %pick random site to move daughter cell to
                end
            end
        elseif pos==2 || pos==4 || pos==6
        %         Site=M_index_vect(randperm(size(M_index_vect,1),1),:);
                Site=M_index_vect(randi(numel(M_index_vect)));
                %Checking neighborhood
            if pos==2
                VacantSites=zeros(1,length(VonNeumann));
                for i=1:length(VonNeumann)
                    ng=VonNeumann(i);
                    if Site+ng<=numel(Chain) && Site+ng>=1 && Chain(Site+ng)==0
                        VacantSites(i)=ng;
                    end
                end
                VacantSites=VacantSites(VacantSites~=0);%exclude any zero in VacantSites
                if isempty(VacantSites)==1
                    randomSite=0;%motion fails
                    pos=7;%do nothing event
                else
                    randomSite=VacantSites(randi(numel(VacantSites))); %pick random site to move to
                end
            end
    end
    
%     Neighbourhood=Neighbours3D(Chain,Site);%checking neighbourhood
    if pos==1 && randomSite+Site>=1 && randomSite+Site<=numel(Chain) && Chain(randomSite+Site)==0%isempty(Neighbourhood)==0 %event = P cell division
%         daughterCell=Neighbourhood(randi(numel(Neighbourhood)));
%         P_index_vect=[P_index_vect;daughterCell]; %update P vector with daughter cell
        daughterCell=randomSite+Site;
        P_index_vect(end+1)=daughterCell;
        Chain(daughterCell)=1; %update matrix with daughter cell
        P=P+1;

    elseif pos == 2 && randomSite+Site>=1 && randomSite+Site<=numel(Chain) && Chain(randomSite+Site)==0%isempty(Neighbourhood)==0 %event = M cell motion
%             newPlace=Neighbourhood(randi(numel(Neighbourhood)));
%             M_index_vect=[M_index_vect;newPlace]; %update vector
%             M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
%             Chain(newPlace)=2;
%             Chain(Site)=0;
            newPlace=randomSite+Site;
            M_index_vect(end+1)=newPlace;
            M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
            Chain(newPlace)=2;
            Chain(Site)=0;
    elseif pos==3 %event = P-M transition
%         P_index_vect=P_index_vect(P_index_vect~=Site);%delete previous site from vector
%         M_index_vect = [M_index_vect;Site]; %add Site to M vector
%         Chain(Site)=2; %update matrix
%         P=P-1;
%         M=M+1;
        
        P_index_vect=P_index_vect(P_index_vect~=Site);%delete previous site from vector
        M_index_vect(end+1)=Site;
        Chain(Site)=2; %update matrix
        P=P-1;
        M=M+1;
    elseif pos==4 %event = M-P transition
%         M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
%         P_index_vect = [P_index_vect;Site]; %add Site to P vector
%         Chain(Site)=1; %update matrix
%         M=M-1;
%         P=P+1;
        M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
        P_index_vect(end+1)=Site;
        Chain(Site)=1; %update matrix
        M=M-1;
        P=P+1;
    elseif pos==5 %event = P cell dies
        P_index_vect=P_index_vect(P_index_vect~=Site);%delete previous site from vector
        Chain(Site)=0;
        P=P-1;
        
    elseif pos==6 %event = M cell dies
        M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
        Chain(Site)=0;
        M=M-1;
        
    end
%     [t, q_m, q_p, M, P, M+P]
% %TRANSFORM MULTIDIMENSIONAL MATRIX INTO COORDINATES X,Y,Z FOR PLOTTING
% [X_p,Y_p,Z_p,X_m,Y_m,Z_m]=multi2coord(Chain);
% %PLOT
% [x,y,z] = meshgrid(1:200,1:200,1:200);
% scatter3(x(:),y(:),z(:),5,Chain(:))
% pause(0.0001)
%     X = Chain(:,:,1);
%     Y = Chain(:,:,2);
%     Z = Chain(:,:,3); 
%     scatter3(X,Y,Z,25,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
%     scatter3(X_p,Y_p,Z_p,25,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
%     hold on
%      scatter3(X_m,Y_m,Z_m,25, ...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[1 1 0]);
%     legend('Proliferating','Moving')
%     set(gca,'XLim',[0 200],'YLim',[0 200],'ZLim',[0 200])
% %     plot3(X_p,Y_p,Z_p,'o');hold on
% %     plot3(X_m,Y_m,Z_m,'o')
%     pause(1)
%     hold off
%     x = Chain(:,:,1);
%     y = Chain(:,:,2);
%     z = Chain(:,:,3); 
%     plot3(x,y,z,'.-','MarkerSize',10)
%     pause(0.00001)
end

% %TRANSFORM MULTIDIMENSIONAL MATRIX INTO COORDINATES X,Y,Z FOR PLOTTING
% [X_p,Y_p,Z_p,X_m,Y_m,Z_m]=multi2coord(Chain);
% %PLOT
%     scatter3(X_p,Y_p,Z_p,25,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
%     hold on
%      scatter3(X_m,Y_m,Z_m,25, ...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[1 1 0]);
%     legend('Proliferating','Moving')
%     set(gca,'XLim',[0 200],'YLim',[0 200],'ZLim',[0 200])
% %     plot3(X_p,Y_p,Z_p,'o');hold on
% %     plot3(X_m,Y_m,Z_m,'o')
%     pause(0.1)
%     hold off
% profile viewer
end