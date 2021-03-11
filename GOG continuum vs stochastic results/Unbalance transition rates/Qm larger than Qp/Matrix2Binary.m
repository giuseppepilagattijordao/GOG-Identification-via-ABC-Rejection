function [binaryMatrix]=Matrix2Binary(matrix)

N1=size(matrix,1);
N2=size(matrix,2);
binaryMatrix=zeros(N1,N2);
    for u=1:N1
        for l=1:N2
            if matrix(u,l)~=0 %if it's not empty
                binaryMatrix(u,l)=1;
            end
        end
    end

end