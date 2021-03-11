function [Cells_index_vect] = get_index_vectors_singlecell(A)
%get linear index vectors from populations P and M in matrix A

Cells_index_vect=[];
for i=1:numel(A)
    if A(i)==1
        Cells_index_vect=[Cells_index_vect;i];
    end
end

%find cells in A (entries equal to 1)
% [row,col]=ind2sub([size(A) size(A)],find(A==1));
% Cells_index_vect=[];
%     for i=1:size(row)
%         Cells_index_vect=[Cells_index_vect ; [row(i) col(i)]];
%     end
% for i=1:size(row)
%     Cells_index_vect=[Cells_index_vect ; [row(i) col(i)]];
% end
end