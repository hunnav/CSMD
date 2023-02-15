options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
[hyp] = fmincon(@(x) -acq_PI(x,S),Initial_theta,[],[],[],[],ones(size(x_sample,1),1)*0,ones(size(x_sample,1),1)*100,[],options);

function [pred_mean] = acq_PI(x,S)

   % EI evaulation

    % correlation r (for speed)
    nSample = size(initial_Dom,2);   % # of the Samples 
    dim = size(initial_Dom,1);       % dim of inputs    
    % Gaussian Process Regression
    r = zeros(nSample,1,dim);
    new_x = new_x';

    for i = 1:dim
        Mxx = new_x(i)^2;
        Kxx = (initial_Dom(i,:).^ 2);
        Nxy = new_x(i)'*initial_Dom(i,:);
        
        r(:,:,i) = Mxx + Kxx' - 2*Nxy';
    end

    r_r = exp(-sum(bsxfun(@times, r, permute(theta, [3,1,2])), 3));

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
        

    if ss > 0
        PI_p =  normcdf(u/ss);    % PI calculation process
    else
        PI_p = 0;
    end

end