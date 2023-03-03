function S = Acq_solver(S)

switch S.acqui.solver
    case 'fmincon'
        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
        S.acqui.x = fmincon(@(x) -Acqfn(x,S), zeros(1,S.prob.dim),[],[],[],[],S.prob.min*ones(S.prob.dim,1),S.prob.max*ones(S.prob.dim,1),[],options)';

    case 'ga'
        options = optimoptions('ga','InitialPopulationRange',[S.prob.min;S.prob.max],'PopulationSize',300,'UseParallel',true);
        S.acqui.x = ga(@(x) -Acqfn(x,S), S.prob.dim,[],[],[],[],S.prob.min*ones(1,S.prob.dim),S.prob.max*ones(1,S.prob.dim),[],options)'; 

    case 'ga+fmincon' 
        hybridopts = optimoptions('fmincon','Display', 'off', 'algorithm', 'interior-point', 'HessianApproximation','bfgs','FiniteDifferenceType', 'central','OptimalityTolerance',1e-10,'UseParallel',true);
        options = optimoptions('ga','InitialPopulationRange',[S.prob.min;S.prob.max],'PopulationSize',300,'UseParallel',true,'HybridFcn',{'fmincon',hybridopts});
        S.acqui.x = ga(@(x) -Acqfn(x,S), S.prob.dim,[],[],[],[],S.prob.min*ones(1,S.prob.dim),S.prob.max*ones(1,S.prob.dim),[],options)';

    case 'pso'
        options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',500,'UseParallel',true,'SelfAdjustmentWeight',2,'SocialAdjustmentWeight',2,'InertiaRange',[0.1,2]);
        S.acqui.x = particleswarm(@(x) -Acqfn(x,S), S.prob.dim,S.prob.min*ones(S.prob.dim,1),S.prob.max*ones(S.prob.dim,1),options)';

    case 'multi_start'
        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
        problem = createOptimProblem('fmincon','objective',@(x) -Acqfn(x,S), 'x0',zeros(S.prob.dim,1),'lb',S.prob.min,'ub',S.prob.max,'options',options)';
        S.acqui.x = run(MultiStart('UseParallel',true),problem,500);
end

end