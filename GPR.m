function [mean_pred,sigma_pred] = GPR(x_sample,y_sample,new_x,theta,sigma,alpha,num_samp,inverse_R_R)

    % Gaussian Process Regression
    Ath = diag(theta);  % make diagonal matrix
    if size(new_x,2) > size(new_x,1)   % if something wrong, make correct
        new_x = new_x';
    end
    
    Mxy = new_x'*Ath*new_x;
    Mxx = diag(Mxy);
    
    Kxy = x_sample'*Ath*x_sample;
    Kxx = diag(Kxy);
    
    Nxy = new_x'*Ath*x_sample;
    
    corr = Mxx + Kxx - 2*Nxy';
    r_r = exp(-corr);          % make correlation matrix r (between unknown point x and known point)
    
    % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
    a = ones(num_samp,1);
    mean_pred = alpha + r_r'*inverse_R_R*(y_sample-alpha*a);
    mu = r_r'*inverse_R_R*a - 1;
    sigma_pred = sigma*(1 - r_r'*inverse_R_R*r_r + mu'*(a'*inverse_R_R*a)\mu);
    % sigma_pred = abs(sigma*(1 - r_r'*inverse_R_R*r_r + mu'*(a'*inverse_R_R*a)\mu));

%     % For simple kriging (remove the last part of the Ordinary kriging)
%     % It can be used only when the mean function is a known deterministic function
    % mean_pred = alpha+r_r'*inverse_R_R*(y_sample-alpha*a);
    % sigma_pred = sigma*(1-r_r'*inverse_R_R*r_r); 
    
%     % unit_test for r_r
    % for i = 1 : num_iter_Samp
    %         rrr(i,1) = exp(-(theta'*(x_sample(:,i)-new_x(:)).^2));
    % end

end