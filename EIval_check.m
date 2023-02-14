function bayesopt_x = EIval_check(initial_Dom,initial_Domy,theta,sigma,alpha,inv_R,low_Range,upper_Range,min_obj,EI_mode,EI_acq_mode,ratio_or_weight)
    num_samp = size(initial_Dom,2); % number of the sample
    dim = size(initial_Dom,1);      % dimension of the sample
    
    switch EI_mode
        case 'ga'
           options = optimoptions('ga','InitialPopulationRange',[low_Range;upper_Range],'PopulationSize',200,'UseParallel',true);
           bayesopt_x = ga(@(x) -EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
               ,dim,[],[],[],[],low_Range*ones(dim,1),upper_Range*ones(dim,1),[],options);       % x point to make EI maximum  

        case 'pso'
           options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
           bayesopt_x = particleswarm(@(x) -EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
               ,dim,low_Range*ones(dim,1),upper_Range*ones(dim,1),options);   

        case 'fmincon'   % only difference is method to find mininum value  
           options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
           bayesopt_x = fmincon(@(x) -EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
               ,zeros(dim,1),[],[],[],[],low_Range*ones(1,dim),upper_Range*ones(1,dim),[],options)';   % different part
          
        case 'ga+fmincon'  % only difference is method to find mininum value 
           hybridopts = optimoptions('fmincon','Display', 'off', 'algorithm', 'interior-point',...
                                     'HessianApproximation','bfgs','FiniteDifferenceType', 'central','OptimalityTolerance',1e-10,'UseParallel',true);   % different part
           options = optimoptions('ga','InitialPopulationRange',[low_Range;upper_Range],'PopulationSize',200,'UseParallel',true,'HybridFcn',{'fmincon',hybridopts});
           bayesopt_x = ga(@(x) -EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
               ,dim,[],[],[],[],low_Range*ones(1,dim),upper_Range*ones(1,dim),[],options);

        case 'multi_start'   % only difference is method to find mininum value
           options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
           problem = createOptimProblem('fmincon','objective',@(x)-EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
               ,'x0',zeros(dim,1),'lb',low_Range,'ub',upper_Range,'options',options);
           ms = MultiStart('UseParallel',true);
           bayesopt_x = run(ms,problem,500);
    
%         case 'random_sampling'    % only difference is method to find mininum value
%             Need to be fixed
%             % Setting for the EI evaluation
%             real_x = lhsdesign(10000,1)*(upper_Range-low_Range)+low_Range;
%             % EI evaulation
%             for i = 1: size(real_x,2)
%                 real_y(i,1) = f(real_x(1,i));
%                 % With new data Real_funcy
%                 [pred_mean,pred_sig] = GPR(initial_Dom,initial_Domy,real_x(1,i),theta,sigma,alpha,num_samp,inv_R);
%                 krig_funcy(i,1) = pred_mean;
%                 krig_devia(i,1) = sqrt(pred_sig);
%                 % Based on the GPR model, EI process is progressed.
%                 u = (min_obj - krig_funcy(i,1));
%                 ss = krig_devia(i,1);
% 
%                 if ss > 0
%                    EI_p(i,1) =  (u*normcdf(u/ss) + ss*normpdf(u/ss));
%                 else
%                    EI_p(i,1) = 0;
%                 end
%             end
%             EIval_entire = EI_p;
%             krig_devia(krig_devia<1e-6) = 0;
     end

end