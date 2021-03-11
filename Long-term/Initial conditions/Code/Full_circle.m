function [Chain] = Full_circle(Chain,N1,N2,c,r)

% r=N/4.162;%radius to get closest to 1800 cells
% c=[round(N/2),round(N/2)];%center of circle
P=0;

for y=1:N1
%     y;
    for x=1:N2
%         x;
%         r^2-(x-c(1))^2;
        if r^2-(x-c(1))^2>0 && y>=-sqrt(r^2-(x-c(1))^2)+c(2) && y<=sqrt(r^2-(x-c(1))^2)+c(2)
            Chain(y,x)=1;
            P=P+1;
        end
    end
end
% imagesc(Chain)
% xlim([0 N])
% ylim([0 N])
% colorbar('Ticks',[0,1,2],...
%      'TickLabels',{'Empty Space','Proliferating','Moving'})
% pause(.1)

end