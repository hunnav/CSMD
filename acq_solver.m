function [S] = acq_solver(S)
% Sigma, alpha ... 는 따로 바꾸지 않음. 사실 EI_acq 또한 S에 있는 정보들로만 구성되게 하면 편할 듯.
% S 안에 있는 식을 코드가 읽고 불러와서 solver안에 들어갈 수 있는지 모르겠음.
% 가령, EI_acq를 S 안에 넣고 bayesopt_x = ga(@(x) -S.acqui.mode.equation 이라고 했을 때,
% 불러와 질 수 있나? 그거는 확인을 해봐야할 듯.
% AA.acqui.mode.equation = @(x) x^2 + x^3 +1;
% bayesopt_x = ga(AA.acqui.mode.equation,1,[],[],[],[],[],[],[],[]);
% 위의 간단한 코드로 실행해본 결과 가능한 것으로 나타남.

switch S.acqui.solver
    case 'ga'
        options = optimoptions('ga','InitialPopulationRange',[low_Range;upper_Range],'PopulationSize',200,'UseParallel',true);
        bayesopt_x = ga(@(x) -EI_acq(S.prob.domain,S.prob.domainy,x,S.Hypopt.theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
            ,dim,[],[],[],[],low_Range*ones(dim,1),upper_Range*ones(dim,1),[],options);       % x point to make EI maximum

    case 'pso'
        options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
        bayesopt_x = particleswarm(@(x) -EI_acq(S.prob.domain,S.prob.domainy,x,S.Hypopt.theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
            ,dim,low_Range*ones(dim,1),upper_Range*ones(dim,1),options);

    case 'fmincon'   % only difference is method to find mininum value
        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
        bayesopt_x = fmincon(@(x) -EI_acq(S.prob.domain,S.prob.domainy,x,S.Hypopt.theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
            ,zeros(dim,1),[],[],[],[],low_Range*ones(1,dim),upper_Range*ones(1,dim),[],options)';   % different part

    case 'ga+fmincon'  % only difference is method to find mininum value
        hybridopts = optimoptions('fmincon','Display', 'off', 'algorithm', 'interior-point',...
            'HessianApproximation','bfgs','FiniteDifferenceType', 'central','OptimalityTolerance',1e-10,'UseParallel',true);   % different part
        options = optimoptions('ga','InitialPopulationRange',[low_Range;upper_Range],'PopulationSize',200,'UseParallel',true,'HybridFcn',{'fmincon',hybridopts});
        bayesopt_x = ga(@(x) -EI_acq(S.prob.domain,S.prob.domainy,x,S.Hypopt.theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
            ,dim,[],[],[],[],low_Range*ones(1,dim),upper_Range*ones(1,dim),[],options);

    case 'multi_start'   % only difference is method to find mininum value
        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
        problem = createOptimProblem('fmincon','objective',@(x)-EI_acq(S.prob.domain,S.prob.domainy,x,S.Hypopt.theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight) ...
            ,'x0',zeros(dim,1),'lb',low_Range,'ub',upper_Range,'options',options);
        ms = MultiStart('UseParallel',true);
        bayesopt_x = run(ms,problem,500);

end

end