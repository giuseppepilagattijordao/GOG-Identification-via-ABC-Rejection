function [P_index_vect, M_index_vect] = get_index_vectors(A)
%get index vectors from populations P and M in matrix A
%with linear indexes
Cells_index_vect=[];
P_index_vect=[];
M_index_vect=[];
for i=1:numel(A)
    if A(i)==1
        P_index_vect=[P_index_vect;i];
    end
    if A(i)==2
        M_index_vect=[M_index_vect;i];
    end 
end
% %find P cells in A (entries equal to 1)
% [rowp,colp]=ind2sub([size(A) size(A)],find(A==1));
% %find M cells in A (entries equal to 2)
% [rowm,colm]=ind2sub([size(A) size(A)],find(A==2));
% P_index_vect=[];
% M_index_vect=[];
%     for i=1:size(rowp)
%         P_index_vect=[P_index_vect ; [rowp(i) colp(i)]];
%     end
% for i=1:size(rowm)
%     M_index_vect=[M_index_vect ; [rowm(i) colm(i)]];
% end
end