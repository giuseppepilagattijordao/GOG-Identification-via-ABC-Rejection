function [Chain]=MakeCircles(Chain,N1,N2,C,R)
%C is a matrix of the coordinates of all the centers of the circles
%R is the vector of each circle

for i=1:length(R)
    r=R(i);%get specific radius
    c=C(i,:);%get specific center
    Chain=Full_circle(Chain,N1,N2,c,r); %draw circle and get Chain
end

end