function [c] = AveragePropagationSpeed(Chain,T)

    if T==0
        c=0;
        return
    end
    N1=size(Chain,1); %number of rows
    N2=size(Chain,2); %number of columns
    %ATTEMPT 6
    %Calculate diagonal distance's average (4 diagonal distances, using 4
    %points that belong to 2 diagonals- a cross on the map)
    Diagonal1=[];
    Diagonal2=[];
    Chain_sites=[];
    for i=1:N1
        Diagonal1=[Diagonal1;i i];
        Diagonal2=[Diagonal2;i N1-i+1];
    end
    for i=1:N1
        for j=1:N2
            if Chain(i,j)~=0
                Chain_sites=[Chain_sites;i j];
            end
        end
    end
    Diagonal1=intersect(Chain_sites,Diagonal1,'rows');
    Diagonal2=intersect(Chain_sites,Diagonal2,'rows');
    if isempty(Diagonal1) || isempty(Diagonal2)
        c=0;
        return
    end
    if length(Diagonal1(:,1))<2 || length(Diagonal2(:,1))<2 %check number of points in each diagonal
        c=0;%there are not enough cells to form a wave front
        return
    end
    
    prev_tlpoint=Diagonal1(1,:);
    prev_brpoint=Diagonal1(1,:);
    prev_blpoint=Diagonal2(1,:);
    prev_trpoint=Diagonal2(1,:);
    for i=2:length(Diagonal1) %get thetop left point and bottom right point
        tlpoint=Diagonal1(i,:);
        brpoint=Diagonal1(i,:);
        if tlpoint(1)<prev_tlpoint(1)%get top left point
            prev_tlpoint=tlpoint;
        end
        if brpoint(1)>prev_brpoint(1)%get bottom right point
            prev_brpoint=brpoint;
        end
    end
    tlpoint=prev_tlpoint;
    brpoint=prev_brpoint;
    for i=2:length(Diagonal2) %get the top right point and bottom left point
        blpoint=Diagonal2(i,:);
        trpoint=Diagonal2(i,:);
        if blpoint(1)>prev_blpoint(1)%get top left point
            prev_blpoint=blpoint;
        end
        if trpoint(1)<prev_trpoint(1)%get bottom right point
            prev_trpoint=trpoint;
        end
    end
    trpoint=prev_trpoint;
    blpoint=prev_blpoint;
    
%     F=[tlpoint; trpoint; blpoint; brpoint];
    c1=sqrt((100-tlpoint(1))^2+(100-tlpoint(2))^2);
    c2=sqrt((100-blpoint(1))^2+(100-blpoint(2))^2);
    c3=sqrt((100-trpoint(1))^2+(100-trpoint(2))^2);
    c4=sqrt((100-brpoint(1))^2+(100-brpoint(2))^2);
    c=mean([c1 c2 c3 c4]);
    c=c/T;
    
    %Attempt 5
    %for each line, get first cell's distance from (100,100) in i j
%     F=zeros(N1,N2);
%     WaveFront=[];
%     %FIRST SIDE SCAN (LEFT TO RIGHT)
%     for i=1:N1
%         j=1;%reset columns
%         while Chain(i,j)==0 && j<N2
%             j=j+1;
%         end
%         if j~=N2
%             WaveFront=[WaveFront;i j];%points of the wave front 
%         end
%     end
%     WaveFront;
%     for k=1:length(WaveFront)
%         point=WaveFront(k,:);
%         F(point(1),point(2))=1;
%         imagesc(F);hold on
%         pause(0.1)
%     end
%     %Scan from other side (RIGHT TO LEFT)
%     for i=1:N1
%         j=N2;%reset columns
%         while Chain(i,j)==0 && j>1
%             j=j-1;
%         end
%         if j~=1
%             WaveFront=[WaveFront;i j];%points of the wave front 
%         end
%     end
%     WaveFront
%     
%     for k=1:length(WaveFront)
%         point=WaveFront(k,:);
%         F(point(1),point(2))=1;
% %         imagesc(F);hold on
% %         pause(0.1)
%     end
%     %CALCULATE THE WAVE FRONT SPEED BASED ON FRONT CELLS
%     D_i=0;
%     D_j=0;
%     for h=1:length(WaveFront)
%         point=WaveFront(h,:);
%        d_j=abs(100-point(2)); 
%        d_i=abs(100-point(1));
%        D_i=D_i+d_i;
%        D_j=D_i+d_j;
%     end
%     D_i=D_i/(length(WaveFront)/2)
%     D_j=D_j/(length(WaveFront)/2)
%     D=(D_i+D_j)/2
%     c=D/T;
    %Attempt 4
    %calculate the middle row and middle column displacement
%     
%     for j=1:N2 
%        if Chain(round(N1/2),j)>0 %if probability of occupancy > 0 from the left
%            leftMostCell=[round(N1/2) j];
%            break
%        end
%     end
%     for j=N2:-1:1
%        if Chain(round(N1/2),j)>0%from the right
%            rightMostCell=[round(N1/2) j];
%            break
%        end
%     end
%      for i=1:N1 
%        if Chain(i,round(N1/2))>0 %if probability of occupancy > 0 from the left
%            topMostCell=[i round(N1/2)];
%            break
%        end
%     end
%     for i=N1:-1:1
%        if Chain(i,round(N1/2))>0%from the right
%            bottomMostCell=[i round(N1/2)]
%            break
%        end
%     end
%     i_displacement=(abs(leftMostCell-round(N1)/2)+abs(rightMostCell-round(N1)/2))/2; %take average of propagation
%     i_speed=i_displacement/T;
%     j_displacement=(abs(topMostCell-round(N2)/2)+abs(bottomMostCell-round(N2)/2))/2; %take average of propagation
%     j_speed=j_displacement/T;
%     c=(i_speed+j_speed)/2;
%     
    %ATTEMPT 3
%Take the average propagation speed in i- an j- direction
%for each cell, get its displacement in i and j direction
% and make an average of the displacement
% and divide by T
%     count=0;
%     D_i=0;
%     D_j=0;
%     for i=1:N1
%         for j=1:N2
%             if Chain(i,j)~=0
%                 d_i=abs(100-i);
%                 d_j=abs(100-j);
%                 D_i=D_i+d_i;
%                 D_j=D_j+d_j;
%                 count=count+1;
%             end
%         end
%     end
%    %average displacement for each direction
%     D_i=D_i/count;
%     D_j=D_j/count;
%     %average displacement
%     D=(D_i+D_j)/2;
%     c=D/T;
%     
    
    
%    % ATTEMPT 2 --- TAKE WAVEFRONT-MOST SITE I AN J
%     WaveFront_i=0;
%     WaveFront_j=0;
%    for j=1:N2
%        for i=1:N1
%            if and(Chain(i,j)~=0,i>WaveFront_i) %if a cell is found even more to the right, consider the front of the wave in j-direction
%                WaveFront_i=i %we have a new Wave front (further to the right)
%            elseif and(Chain(i,j)~=0,j>WaveFront_j)%if a cell is found even more to the right, consider the front of the wave in i-direction
%                 WaveFront_j=j %we have a new Wave front (further down)
%            end
%        end
%    end
%       WaveFront_i=WaveFront_i-100;%get displacement from site (100,100)
%       WaveFront_j=WaveFront_j-100;
%       c=(WaveFront_i+WaveFront_j)/2; %average of the 2 directions
%       c=c/T;
    
    
    %ATTEMPT 1 - average number of occupied cell's per row divided by 2
%     row_avg=0;
%     column_avg=0;
%     
%     for i=1:N1
% %         row_avg=row_avg+n/nnz(Chain_avg(:,round(length(Chain_avg)/2)));
%         row_avg=row_avg+nnz(Chain(i,:))/2;
%     end
%     for j=1:N2
%         n=nnz(Chain(:,j))/2;%number of cells in row i
% %         column_avg=column_avg+n/nnz(Chain_avg(round(length(Chain_avg)/2),:));
%         column_avg=column_avg+n;
%     end
%     row_avg;
%     column_avg;
%     row_avg=row_avg/column_avg
%     column_avg=column_avg/row_avg
%     c=(row_avg+column_avg)/2
%     c=c/T
end