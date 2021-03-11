function [d]=distance(Y_obs,Y_sim,theta)
    r=rand();
    if r<0.5
        d=abs(Y_sim(1));
        return
    end
    if r>=0.5
        d=abs(mean(Y_sim));
        return
    end
end