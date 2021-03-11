
 function [M,P,Chain,Chain50,Chain100,Chain200,Q_ij,WaveSpeed,TumorSize_cellcycles,WaveSpeed_vect] = Gerlee_Main_loop(Chain,alpha,v,q_p,q_m,mu,T,t)
% profile on
[P_index_vect, M_index_vect]=get_index_vectors(Chain);%getting indexes of each species in the matrix
% M_vector=[];
% P_vector=[];
P=sum(Chain(:)==1);%count number of P cells
M=sum(Chain(:)==2);%count number of M cells
Q_ij=Matrix2Binary(Chain);%matrix to count number of occupancies of each site
cnt=0; %to count number of reactions
WaveSpeed_avg=0;
TumorSize_cellcycles=nnz(Chain)/numel(Chain);%for data collection
WaveSlope_vect=[];
WaveSpeed_vect=[];
if T==0
    WaveSpeed_vect=[0];
end
nrows=size(Chain,1);
VonNeumann=[-1,1,nrows,-nrows];
randomWalkNeighbours=[-nrows,-nrows+1,1,nrows+1,nrows,nrows-1,-1,-nrows-1];
eventcount=0;
c1=0;
c2=0;
c3=0;
c4=0;
c5=0;
c6=0;
c7=0;
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
    r=rand(1,2);
    r1=r(1);
    r2=r(2);
    %TIME INCREMENT
    a0=sum(a);
%     a0=a1+a2+a3+a4+a5+a6;
    dt=(1/(a0))*log(1/r1);
    
    %CHOOSE (POSSIBLE) EVENT TO PERFORM  
    %if r2 falls between 0 and a1/a0, choose a1, if it falls between (a1+a2)/a0 choose a2, etc...
    if r2>=0 && r2<=a(1)/a0
        pos=1;
    elseif r2>a(1)/a0 && r2<=(a(1)+a(2))/a0
        pos=2;
    elseif r2>(a(1)+a(2))/a0 && r2<=(a(1)+a(2)+a(3))/a0
        pos=3;
    elseif r2>(a(1)+a(2)+a(3))/a0 && r2<=(a(1)+a(2)+a(3)+a(4))/a0
        pos=4;
    elseif r2>(a(1)+a(2)+a(3)+a(4))/a0 && r2<=(a(1)+a(2)+a(3)+a(4)+a(5))/a0
        pos=5;
    elseif r2>(a(1)+a(2)+a(3)+a(4)+a(5))/a0 && r2<=1
        pos=6;
    else
        display('here!')
        break
    end
%     pos=GillespieReactionChoice(a,r2);
%     pos = min(find(r2*a0 <= cumsum(a)));
    %UPDATE TIME
    t_previous=t;%use this time to count integer times for Q(i,i,t) function
    t=t+dt;
    %-----------------------------------------------------------
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
%         elseif pos==5
%                 All_Cells_vect=[P_index_vect M_index_vect];
%                 Site=All_Cells_vect(randi(numel(All_Cells_vect)));
%             if P>0 && M>0 %choose any if there are both types of cells
%                 r3=rand();
%             elseif M==0 %choose only P cells
%                 r3=0;
%             elseif P==0 %choose only M cells
%                 r3=1;
%             end
% 
%             if r3<0.5 && P>0
%                 Site=P_index_vect(randi(numel(P_index_vect)));
%             elseif r3>=0.5 && M>0
%                 Site=M_index_vect(randi(numel(M_index_vect)));
%             end
    else
        display('here2!')
        break
    end

%     Neighbourhood=Neighbours(Chain,Site);%checking neighbourhood
    if pos==1 && randomSite+Site>=1 && randomSite+Site<=numel(Chain) && Chain(randomSite+Site)==0%isempty(Neighbourhood)==0 %event = P cell division
%         daughterCell=Neighbourhood(randi(numel(Neighbourhood)));
        daughterCell=randomSite+Site;
%         P_index_vect=[P_index_vect;daughterCell]; %update P vector with daughter cell
        P_index_vect(end+1)=daughterCell;
        Chain(daughterCell)=1; %update matrix with daughter cell
        P=P+1;
%         c1=c1+1;
%         dt=exprnd(alpha)/a0;
%         t=t+dt;
        eventcount=eventcount+1;
    elseif pos == 2 && randomSite+Site>=1 && randomSite+Site<=numel(Chain) && Chain(randomSite+Site)==0%isempty(Neighbourhood)==0 %event = M cell motion
%             newPlace=Neighbourhood(randi(numel(Neighbourhood)));
            newPlace=randomSite+Site;
%             M_index_vect=[M_index_vect;newPlace]; %update vector
            M_index_vect(end+1)=newPlace;
            M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
            Chain(newPlace)=2;
            Chain(Site)=0;
%             c2=c2+1;
%             dt=exprnd(v)/a0;
%             t=t+dt;
            eventcount=eventcount+1;
    elseif pos==3 %event = P-M transition
        P_index_vect=P_index_vect(P_index_vect~=Site);%delete previous site from vector
        
%         M_index_vect = [M_index_vect;Site]; %add Site to M vector
        M_index_vect(end+1)=Site;
        Chain(Site)=2; %update matrix
        P=P-1;
        M=M+1;
%         c3=c3+1;
%         dt=exprnd(q_m)/a0;
%         t=t+dt;
        eventcount=eventcount+1;
    elseif pos==4 %event = M-P transition
        M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
%         P_index_vect = [P_index_vect;Site]; %add Site to P vector
        P_index_vect(end+1)=Site;
        Chain(Site)=1; %update matrix
        M=M-1;
        P=P+1;
%         c4=c4+1;
%         dt=exprnd(q_p)/a0;
%         t=t+dt;
        eventcount=eventcount+1;
%     elseif pos==5 %event = cell dies
%         if Chain(Site)==1
%             P=P-1;
%             P_index_vect=P_index_vect(P_index_vect~=Site);%delete previous site from vector
%         elseif Chain(Site)==2
%             M=M-1;
%             M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
%         end
%         Chain(Site)=0;
    elseif pos==5 %event = P cell dies
        P_index_vect=P_index_vect(P_index_vect~=Site);%delete previous site from vector
          Chain(Site)=0;
          P=P-1;
%           if ismember(Site,P_index_vect)==1 %Death of P-cell
%               P=P-1;
%               P_index_vect=P_index_vect(P_index_vect~=Site);%delete previous site from vector
%           elseif ismember(Site,M_index_vect)==1
%               M=M-1;
%               M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
%           end
        
% %         c5=c5+1;
% %         dt=exprnd(mu);
% %         t=t+dt;
    elseif pos==6 %event = M cell dies
        M_index_vect=M_index_vect(M_index_vect~=Site);%delete previous site from vector
        Chain(Site)=0;
        M=M-1;
%         c6=c6+1;
        eventcount=eventcount+1;
%     elseif pos==7 %do nothing
%         t=t-dt; %get the time back and try again
%         c7=c7+1;
    end
%     [t, q_m, q_p, M, P, M+P]
    %-------------------------------------------------------------
    %PLOT
%     colormap hot
%     imagesc(Chain)
%     title(['Total nr of cells: ' num2str(P+M) ', T = ' num2str(t)])
%     %colormap(gca,[0,1,0;0,1,1;1,1,1])
%     colorbar('Ticks',[0,1,2],...
%          'TickLabels',{'Empty Space','Proliferating','Moving'})
%     caxis([0,2]) %to refrain the colobar scale values from moving in the plot
% %     pause(1)
% % plot(Chain(100,:))
% pause(0.00000001)

%CALCULATE Q(i,j,t)
%FOR SHORT-TERM
%     if t_previous<1/6 && t>1/6
%         Chain50=Chain;
%     end
%     if t_previous<1 && t>1
%         Chain100=Chain;
%     end
%     if t_previous<48/36 && t>48/36
%         Chain200=Chain;
%     end
%FOR LONG-TERM
%     if t_previous<14 && t>14
%         Chain50=Chain;
%     end
%     if t_previous<20 && t>20
%         Chain100=Chain;
%     end
%     if t_previous<26 && t>26
%         Chain200=Chain;
%     end
%SHORT-TERM 2
%     if t_previous<1/3 && t>1/3
%         Chain50=Chain;
%     end
%     if t_previous<1 && t>1
%         Chain100=Chain;
%     end
%     if t_previous<2 && t>2
%         Chain200=Chain;
%     end
        
if abs(floor(10*t_previous)-floor(10*t))==1%whenever there's a new cycle
    
    
%     [t eventcount P+M]
eventcount=0;
Q_ij=Q_ij+Matrix2Binary(Chain);
Q_ij=Q_ij/max(max(Q_ij));
% WaveSpeed_avg=WaveSpeed_avg+AveragePropagationSpeed(Q_ij/cnt,cnt);
% [t WaveSpeed_avg/cnt AveragePropagationSpeed(Q_ij,cnt)]

%     imagesc(Q_ij)
%     colormap(gca,[0,1,0;0,1,1;1,1,1])
%     colorbar('Ticks',[0,1,2],...
%          'TickLabels',{'Empty Space','Proliferating','Moving'})
%     caxis([0,2]) %to refrain the colobar scale values from moving in the plot
%     pause(.1)
TumorSize_cellcycles=[TumorSize_cellcycles nnz(Chain)/numel(Chain)];%for data collection
WaveSpeed_vect=[WaveSpeed_vect AveragePropagationSpeed(Q_ij,T)];
% WaveSpeed_vect=[WaveSpeed_vect AveragePropagationSpeed2(Q_ij,T)];
% WaveSlope_vect=[WaveSlope_vect max(gradient(TumorSize_cellcycles))];
cnt=cnt+1;
end
end
% c=[c1,c2,c3,c4,c5,c6];
%occupancy probability
% Q_ij=Q_ij/cnt;
% WaveSpeed=1;
WaveSpeed=AveragePropagationSpeed(Q_ij,T);%Gerlee IC
% WaveSpeed=AveragePropagationSpeed2(Q_ij,T);%Simpson IC
% TumorSize=M_vector(end)+P_vector(end);
%     figure
%     imagesc(Chain)
%     %colormap(gca,[0,1,0;0,1,1;1,1,1])
%     colorbar('Ticks',[0,1,2],...
%          'TickLabels',{'Empty Space','Proliferating','Moving'})
%     caxis([0,2]) %to refrain the colobar scale values from moving in the plot
%     pause(.1)
% profile viewer
Chain50=1;
Chain100=1;
Chain200=1;
end