function EI_p = EI_acq(initial_Dom,initial_Domy,new_x,theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight)  
    % EI evaulation
    r_r = exp(-sum(bsxfun(@times, Correlation(initial_Dom,new_x), permute(theta, [3,1,2])), 3));

    % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
    a = ones(num_samp,1);
    pred_mean = alpha + r_r'*inv_R*(initial_Domy-alpha);
    mu = r_r'*inv_R*a - 1;
    pred_sig = sigma*(1 - r_r'*inv_R*r_r + mu'*(a'*inv_R*a)\mu);
    
    % sigma_pred = abs(sigma*(1 - r_r'*inverse_R_R*r_r + mu'*(a'*inverse_R_R*a)\mu));

    % % For simple kriging (remove the last part of the Ordinary kriging)
    % % It can be used only when the mean function is a known deterministic function
    % mean_pred = alpha+r_r'*inv_R*(y_sample-alpha*a);
    % sigma_pred = sigma*(1-r_r'*inv_R*r_r); 
    
    % % unit_test for r_r
    % for i = 1 : num_iter_Samp
    %         rrr(i,1) = exp(-(theta'*(initial_Dom(:,i)-new_x(:)).^2));
    % end

    % Based on the GPR model, EI process is progressed.
    u = (min_obj - pred_mean);     
    ss = sqrt(pred_sig);             
        
    p0 = 10^(floor(1+log10(abs(u))));

    switch EI_acq_mode
        case 'normal'
            if ss > 0
                EI_p =  ((u-p0*ratio_or_weight)*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
            else
                EI_p = 0;
            end
        case 'weighted'
            % ref. On the Design of Optimization Strategies Based on...
            % Global Response Surface Approximation Models (2005)
            % Journal of Global Optimization
            if ss > 0
                EI_p =  (ratio_or_weight*u*normcdf(u/ss) + ss*(1-ratio_or_weight)*normpdf(u/ss));   % EI calculation process
            else
                EI_p = 0;
            end
        otherwise
            error('Please select one')
    end
end