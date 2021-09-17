function pos = OU_diff_drift_every_2_time_points(len,N, pos0,drift_dir,step,rho,mu)
% skip a step ever other step
    pos1_x = pos0(:,1)+step.*randn(N,1);
    pos1_y = pos0(:,2)+step.*randn(N,1);
    poscurr(:,1,1) = pos1_x;
    poscurr(:,1,2) = pos1_y;
    if drift_dir==1
        theta = 2 * pi * rand(N,1);   % mean velocity orientation
        for i=2:(2*len)
            poscurr(:,i,1) = rho.*(poscurr(:,i-1,1)-pos0(:,1)-(i-2)*mu*cos(theta))+pos0(:,1)+(i-1)*mu*cos(theta)+step*sqrt((1-rho^2)).*randn(N,1);
            poscurr(:,i,2) = rho.*(poscurr(:,i-1,2)-pos0(:,2)-(i-2)*mu*sin(theta))+pos0(:,2)+(i-1)*mu*sin(theta)+step*sqrt((1-rho^2)).*randn(N,1);
%             poscurr(:,i,1) = rho.*(poscurr(:,i-1,1)-pos1_x-(i-2)*mu.*cos(theta))+pos1_x+(i-1)*mu.*cos(theta)+step*sqrt((1-rho^2)).*randn(N,1);
%             poscurr(:,i,2) = rho.*(poscurr(:,i-1,2)-pos1_y-(i-2)*mu.*sin(theta))+pos1_y+(i-1)*mu.*sin(theta)+step*sqrt((1-rho^2)).*randn(N,1);
        end
    else
        for i=2:(2*len)
            poscurr(:,i,1) = rho.*(poscurr(:,i-1,1)-pos0(:,1)-(i-2)*mu)+pos0(:,1)+(i-1)*mu+step*sqrt((1-rho^2)).*randn(N,1);
            poscurr(:,i,2) = rho.*(poscurr(:,i-1,2)-pos0(:,2)-(i-2)*mu)+pos0(:,2)+(i-1)*mu+step*sqrt((1-rho^2)).*randn(N,1);
%             poscurr(:,i,1) = rho.*(poscurr(:,i-1,1)-pos1_x-(i-2)*mu)+pos1_x+(i-1)*mu+step*sqrt((1-rho^2)).*randn(N,1);
%             poscurr(:,i,2) = rho.*(poscurr(:,i-1,2)-pos1_y-(i-2)*mu)+pos1_y+(i-1)*mu+step*sqrt((1-rho^2)).*randn(N,1);
        end 
    end
    pos = poscurr(:,(1:len)*2,:);
end

