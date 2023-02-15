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
                switch S.acqui.solver
                    case 'ga'
                        options = optimoptions('ga','InitialPopulationRange',[low_Range;upper_Range],'PopulationSize',200,'UseParallel',true);
                        bayesopt_x = ga(AA.acqui.mode.equation,dim,[],[],[],[],low_Range*ones(dim,1),upper_Range*ones(dim,1),[],options);       % x point to make EI maximum

                    case 'pso'
                        options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
                        bayesopt_x = particleswarm(AA.acqui.mode.equation,dim,low_Range*ones(dim,1),upper_Range*ones(dim,1),options);

                    case 'fmincon'   % only difference is method to find mininum value
                        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
                        bayesopt_x = fmincon(AA.acqui.mode.equation,zeros(dim,1),[],[],[],[],low_Range*ones(1,dim),upper_Range*ones(1,dim),[],options)';   % different part

                    case 'ga+fmincon'  % only difference is method to find mininum value
                        hybridopts = optimoptions('fmincon','Display', 'off', 'algorithm', 'interior-point', 'HessianApproximation','bfgs','FiniteDifferenceType', 'central','OptimalityTolerance',1e-10,'UseParallel',true);   % different part
                        options = optimoptions('ga','InitialPopulationRange',[low_Range;upper_Range],'PopulationSize',200,'UseParallel',true,'HybridFcn',{'fmincon',hybridopts});
                        bayesopt_x = ga(AA.acqui.mode.equation,dim,[],[],[],[],low_Range*ones(1,dim),upper_Range*ones(1,dim),[],options);

                    case 'multi_start'   % only difference is method to find mininum value
                        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
                        problem = createOptimProblem('fmincon','objective',@(x)-acq_EI(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio),'x0',zeros(dim,1),'lb',low_Range,'ub',upper_Range,'options',options);
                        bayesopt_x = run(MultiStart('UseParallel',true),problem,500);
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



% @(x) -acq_EI(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio)


function EI_p = acq_EI(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio)  
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


function PI_p = acq_PI(Domain,Domain_y,x,theta,sigma,alpha_kriging,num_samp,inv_R,min_obj,ratio)  
    % PI evaulation

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

    if ss > 0
        PI_p =  normcdf(u/ss);    % PI calculation process
    else
        PI_p = 0;
    end

end

function [pred_mean] = acq_MP(x,S)
    % MP evaulation

    r_r = exp(-sum(bsxfun(@times, Correlation(Domain,x), permute(theta, [3,1,2])), 3));
    
    % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
    a = ones(num_samp,1);
    pred_mean = alpha + r_r'*inv_R*(initial_Domy-alpha);

    const_group = [g1,g2,g3,g4]; % Group the constraints from the design varaible
                                 % Make sure that g value should be minus
                                 % if design variable violate the given constraints
    vio_ind = find(find(const_group < 0)>0); % Find which constraints are violated with the design variable
    
    if ismepty(vio_ind) ~= 1
        pred_mean = pred_mean + sum(abs([g1,g2,g3,g4]).*vio_ind);
    else
    end
    % 현재 constraints 가 위반 할 시, constraint ( g1,g2,g3,g4) 의 값이 음수값을 갖도록 통일하면
    % 코드가 되게 편해질 듯. ( 각 constraints 마다 매번 처리를 해줘야 하지만,
    % 프리 스로세스로 해줄 수 있다고 개인적으로 생각이 듬.
    % 그리고 마지막으로 MP+EI process는 직접 code를 짜서 진행해보면 좋을 듯.
    % 일반적으로 EI+MP를 많이 쓰는 것으로 인지하고 있음.
    % 짜기전에 나머지 code 다 정리하고(S로), EI+MP 어떻게 짜야하는지 와서 물어보면 될 듯.
end