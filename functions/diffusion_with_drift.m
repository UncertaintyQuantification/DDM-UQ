function pos = diffusion_with_drift(len,step,pos0,N,mu)
% Simulates the brownian motion with drift of particles

    %pos = zeros(N,len,2);
    pos(:,1,1) = pos0(:,1); % x position at time 0
    pos(:,1,2) = pos0(:,2); % y position at time 0
    
    for i=2:len
        pos(:,i,1) = pos(:,i-1,1)+normrnd(mu,step,N,1);
        pos(:,i,2) = pos(:,i-1,2)+normrnd(mu,step,N,1);
    end   
end