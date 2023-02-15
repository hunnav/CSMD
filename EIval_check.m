function [bayesopt_x,r] = EIval_check(Domain,Domain_y,theta,sigma,alpha_kriging,inv_R,low_Range,upper_Range,ratio,beta,num_initial_value, Minimum_Value,EI_mode)
    num_samp = size(Domain,2); % number of the sample
    dim = size(Domain,1);      % dimension of the sample

    if  Minimum_Value(end,2) == inf
        min_obj = min(Domain_y);
    else
        min_obj = Domain_y(num_initial_value + Minimum_Value(end,1));
    end

    while 1
        while 1
            try
                switch EI_mode
                    case 'ga'
                        options = optimoptions('ga','InitialPopulationRange',[low_Range;upper_Range],'PopulationSize',200,'UseParallel',true);
                        bayesopt_x = ga(@(x) -EI_acq(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio) ...
                            ,dim,[],[],[],[],low_Range*ones(dim,1),upper_Range*ones(dim,1),[],options);       % x point to make EI maximum

                    case 'pso'
                        options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
                        bayesopt_x = particleswarm(@(x) -EI_acq(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio) ...
                            ,dim,low_Range*ones(dim,1),upper_Range*ones(dim,1),options);

                    case 'fmincon'   % only difference is method to find mininum value
                        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
                        bayesopt_x = fmincon(@(x) -EI_acq(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio) ...
                            ,zeros(dim,1),[],[],[],[],low_Range*ones(1,dim),upper_Range*ones(1,dim),[],options)';   % different part

                    case 'ga+fmincon'  % only difference is method to find mininum value
                        hybridopts = optimoptions('fmincon','Display', 'off', 'algorithm', 'interior-point',...
                            'HessianApproximation','bfgs','FiniteDifferenceType', 'central','OptimalityTolerance',1e-10,'UseParallel',true);   % different part
                        options = optimoptions('ga','InitialPopulationRange',[low_Range;upper_Range],'PopulationSize',200,'UseParallel',true,'HybridFcn',{'fmincon',hybridopts});
                        bayesopt_x = ga(@(x) -EI_acq(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio) ...
                            ,dim,[],[],[],[],low_Range*ones(1,dim),upper_Range*ones(1,dim),[],options);

                    case 'multi_start'   % only difference is method to find mininum value
                        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
                        problem = createOptimProblem('fmincon','objective',@(x)-EI_acq(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio) ...
                            ,'x0',zeros(dim,1),'lb',low_Range,'ub',upper_Range,'options',options);
                        ms = MultiStart('UseParallel',true);
                        bayesopt_x = run(ms,problem,500);
                end
                bayesopt_x = bayesopt_x';
                break
            catch
                disp('There is some error, repeat again.')
            end
        end
        r = Correlation(Domain,bayesopt_x);
        if min(sum(r, 3)) > beta
            break
        else
            if abs(ratio) >= 5
                break
            else
                ratio = ratio+1;
            end
        end
    end
end


function EI_p = EI_acq(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio)  
    % EI evaulation

    r_r = exp(-sum(bsxfun(@times, Correlation(Domain,x), permute(theta, [3,1,2])), 3));

    % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
    a = ones(num_samp,1);
    pred_mean = alpha_kriging + r_r'*inv_R*(Domain_y-alpha_kriging);
    mu = r_r'*inv_R*a - 1;
    pred_sig = sigma*(1 - r_r'*inv_R*r_r + mu'*(a'*inv_R*a)\mu);
    
    % % For simple kriging (remove the last part of the Ordinary kriging)
    % % It can be used only when the mean function is a known deterministic function
    % mean_pred = alpha+r_r'*inv_R*(y_sample-alpha*a);
    % sigma_pred = sigma*(1-r_r'*inv_R*r_r); 
    
    % % unit_test for r_r
    % for i = 1 : num_iter_Samp
    %         rrr(i,1) = exp(-(theta'*(initial_Dom(:,i)-new_x(:)).^2));
    % end

    u = (min_obj - pred_mean);     
    ss = sqrt(pred_sig);              
    p0 = 10^(floor(1+log10(abs(u))));

    if ss > 0
        EI_p =  ((u-p0*ratio)*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
    else
        EI_p = 0;
    end

end