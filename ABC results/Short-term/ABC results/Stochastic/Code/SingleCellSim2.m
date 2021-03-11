function [Chain,Chain0,Chain500,Chain1000,Chain2000,Cells_index_vect,WaveSpeed,Density_history,c2,c,WaveSpeed_vect] = SingleCellSim2(T,alpha,v,mu,Chain,tau)
%     profile on
    
Cells_index_vect=get_index_vectors_singlecell(Chain);%getting linear indexes of each species in the matrix
t=0;
Cells=nnz(Chain); %number of cells
Chain0=Chain;
nrows=size(Chain,1);
VonNeumann=[-1,1,-nrows,nrows];
%for wave speed
WaveSpeed=0;
Q_ij=zeros(size(Chain,1),size(Chain,2));
Density_history=nnz(Chain)/numel(Chain);
WaveSpeed_vect=[];
c=0;
c2=0;
if t==T
    WaveSpeed_vect=[0];
end
    while t<T
%          t
        if Cells==0
            break
        end
        
        %MOTILITY: choose each cell on lattice randomnly
        dummie_vect=Cells_index_vect;
        n=0;
        while n<length(Cells_index_vect)%isempty(dummie_vect)==0
            c=c+1;
%             Cell=dummie_vect(randi(numel(dummie_vect)));%pick random cell
%             dummie_vect=dummie_vect(dummie_vect~=Cell);%delete previous cell from vector to avoid repetition 
            Cell=Cells_index_vect(randi(numel(Cells_index_vect)));%pick random cell
            randomSite=VonNeumann(randi(numel(VonNeumann))); %pick random site to move to
%             Neighborhood=Neighbours(Chain,Cell); %check neighborhood
            %probability to move, v
            r=rand();
            [row col]=ind2sub([size(Chain,1) size(Chain,2)],Cell);%sub position of current cell
%             [row2 col2]=ind2sub([size(A,1) size(A,2)],Cell+randomSite);%sub position of site candidate
            if r<v && not(row==size(Chain,1) && randomSite==1) && not(row==1 && randomSite==-1) && not(col==size(Chain,2) && randomSite==nrows) && not(col==1 && randomSite==-nrows) && randomSite+Cell>=1 && randomSite+Cell<=numel(Chain) && Chain(randomSite+Cell)==0 %isempty(Neighborhood)==0 
%                 newPlace=Neighborhood(randi(numel(Neighborhood)));%pick random place to move to
                newPlace=randomSite+Cell;
                Chain(newPlace)=1;
                Chain(Cell)=0;
                Cells_index_vect=[Cells_index_vect; newPlace];%add new location
                Cells_index_vect=Cells_index_vect(Cells_index_vect~=Cell);%remove previous location
                c2=c2+1;
%                 imagesc(Chain)
%                 pause(1)
            end
            n=n+1;
        end
        %PROLIFERATION: choose each cell on lattice randomly
        dummie_vect2=Cells_index_vect;
        n=0;
        while n<length(Cells_index_vect)%isempty(dummie_vect2)==0
%             Cell=dummie_vect2(randi(numel(dummie_vect2)));%pick random cell
%             dummie_vect2=dummie_vect2(dummie_vect2~=Cell);%delete previous cell from vector to avoid repetition 
            Cell=Cells_index_vect(randi(numel(Cells_index_vect)));
            randomSite=VonNeumann(randi(numel(VonNeumann))); %pick random site to move to
%             Neighborhood=Neighbours(Chain,Cell); %check neighborhood
            %probability to proliferate
            r=rand();
            if r<alpha && randomSite+Cell>=1 && randomSite+Cell<=numel(Chain) && Chain(randomSite+Cell)~=1 %isempty(Neighborhood)==0 
%                 daughterCell=Neighborhood(randi(numel(Neighborhood)));%pick random place to move to
                daughterCell=randomSite+Cell;
                Chain(daughterCell)=1;
                Cells_index_vect=[Cells_index_vect; daughterCell];%add new cell location
            end
            n=n+1;
        end
        
        %DEATH
        r=rand();
        if r<mu
            Cell=Cells_index_vect(randi(numel(Cells_index_vect)));%pick random cell
            Chain(Cell)=0;%cell dies
        end
        
        %update t
        t_previous=t;
        
        t=t+tau;
        %FOR SHORT-TERM
%         if t_previous<1/6 && t>1/6
%             Chain500=Chain;
%         end
%         if t_previous<1 && t>1
%             Chain1000=Chain;
%         end
%         if t_previous<1.5 && t>1.5
%             Chain2000=Chain;
%         end
% FOR LONG-TERM
        if t_previous<14 && t>14
            Chain500=Chain;
        end
        if t_previous<20 && t>20
            Chain1000=Chain;
        end
        if t_previous<26 && t>26
            Chain2000=Chain;
        end
%------------------------------------------------------------------
%SHORT-TERM 2
%         if t_previous<1/3 && t>1/3
%             Chain500=Chain;
%         end
%         if t_previous<1 && t>1
%             Chain1000=Chain;
%         end
%         if t_previous<2 && t>2
%             Chain2000=Chain;
%         end
        
        
%         if floor(t_previous)==499 && floor(t)==500
%             Chain500=Chain;
%         elseif floor(t_previous)==999 && floor(t)==1000
%             Chain1000=Chain;
%         elseif floor(t_previous)==1999 && floor(t)==2000
%             Chain2000=Chain;
%         end
               
        %WaveSpeed:
        Q_ij=Q_ij+Matrix2Binary(Chain);
        Q_ij=Q_ij/max(max(Q_ij));
%         WaceSpeed=1;
        WaveSpeed=AveragePropagationSpeed(Q_ij,T);%Gerlee IC
%         WaveSpeed=AveragePropagationSpeed2(Q_ij,T);%Simpson IC
        if abs(floor(10*t_previous)-floor(10*t))==1%whenever there's a new cycle
        %Density history
        Density_history=[Density_history nnz(Chain)/numel(Chain)];
            WaveSpeed_vect=[WaveSpeed_vect AveragePropagationSpeed(Q_ij,T)];
%         WaveSpeed_vect=[WaveSpeed_vect AveragePropagationSpeed2(Q_ij,T)];
        end
%         imagesc(Chain)
%         colorbar()
% %         plot(Chain(round(size(Chain,1)),:))
%         pause(0.001)
    end
    
    %In case T doesnt go to these values
        Chain500=1;
        Chain1000=1;
        Chain2000=1;
end