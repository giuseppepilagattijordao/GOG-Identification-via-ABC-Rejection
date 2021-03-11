%-------------------------
%Initial condition filled circle
N1=50;
N2=50;
Chain=zeros(N1,N2);
Chain2=Chain;
C=[round(N1/2),round(N1/2)];
R=ones(size(C,1),1)*5;%radius of circles
Chain1=MakeCircles(Chain,N1,N2,C,R);
for i=1:numel(Chain1)
    if Chain1(i)==1
        Chain2(i)=1+round(rand()); %randomly assign motile/proliferative to each position of initial condition
    end
end
% 
figure
imagesc(Chain1)
colormap(gca,[1,1,1;0,0,0;0,0,0])
title('SS')
xlabel('j')
ylabel('i')
figure
imagesc(Chain2)
colormap(gca,[1,1,1;0,0,1;1,0,0])
title('GOG')
xlabel('j')
ylabel('i')
