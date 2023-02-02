function bayesopt_x = EIval_check(EI_mode,initial_Dom,initial_Domy,theta,sigma,alpha,inv_R,low_Range,upper_Range,min_obj)
    num_samp = size(initial_Dom,2); % number of the sample
    dim = size(initial_Dom,1);      % dimension of the sample
    
    switch EI_mode
        case 'ga'
           EI_p = @(x) -EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj);        
           options = optimoptions('ga','InitialPopulationRange',[low_Range;upper_Range],'PopulationSize',150,'UseParallel',true); % what....????????
           ga_x = ga(EI_p,dim,[],[],[],[],low_Range,upper_Range,[],options);                                % x point to make EI maximum
           bayesopt_x = ga_x;                                      % maximum value of function EI

        case 'optimoptions'
           EI_p = @(x) -EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj);
           options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
           particleswarm_x = particleswarm(EI_p,dim,low_Range*ones(dim,1),upper_Range*ones(dim,1),options);
           bayesopt_x = particleswarm_x;     


        case 'fmincon'   % only difference is method to find mininum value
           EI_p = @(x) -EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj);
           options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
           fmin_x = fmincon(EI_p,zeros(dim,1),[],[],[],[],low_Range,upper_Range,[],options);   % different part
           bayesopt_x = fmin_x;
          
        case 'ga+fmincon'  % only difference is method to find mininum value
           EI_p = @(x) -EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj);
           hybridopts = optimoptions('fmincon','Display', 'off', 'algorithm', 'interior-point',...
                                     'HessianApproximation','bfgs','FiniteDifferenceType', 'central','OptimalityTolerance',1e-10);   % different part
           options = optimoptions('ga','InitialPopulationRange',[low_Range,upper_Range],'PopulationSize',500,'UseParallel',true,'HybridFcn',{'fmincon',hybridopts});
           [hybrid_x,EI_p] = ga(EI_p,dim,[],[],[],[],low_Range,upper_Range,[],options);
           bayesopt_x = hybrid_x;

        case 'multi_start'   % only difference is method to find mininum value
           EI_p = @(x) -EI_acq(initial_Dom,initial_Domy,x,theta,sigma,alpha,num_samp,inv_R,min_obj);
           options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
           problem = createOptimProblem('fmincon','objective',EI_p,'x0',zeros(dim,1),'lb',low_Range,'ub',upper_Range,'options',options);
           ms = MultiStart('UseParallel',true);
           [multi_x,EI_p] = run(ms,problem,500);
           bayesopt_x = multi_x;
    
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