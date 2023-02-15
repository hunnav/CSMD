function S = Acq_solver(S)

switch S.acqui.solver
    case 'ga'
        options = optimoptions('ga','InitialPopulationRange',[S.prob.min;S.prob.max],'PopulationSize',200,'UseParallel',true);
        S.acqui.x = ga(@(x) -Acqfn(S),S.prob.dim,[],[],[],[],S.prob.min*ones(S.prob.dim,1),S.prob.max*ones(S.prob.dim,1),[],options);       % x point to make EI maximum

    case 'pso'
        options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
        S.acqui.x = particleswarm(@(x) -Acqfn(S),S.prob.dim,S.prob.min*ones(S.prob.dim,1),S.prob.max*ones(S.prob.dim,1),options);

    case 'fmincon'   % only difference is method to find mininum value
        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
        S.acqui.x = fmincon(@(x) -Acqfn(S),zeros(S.prob.dim,1),[],[],[],[],S.prob.min*ones(1,S.prob.dim),S.prob.max*ones(1,dim),[],options)';   % different part

    case 'ga+fmincon'  % only difference is method to find mininum value
        hybridopts = optimoptions('fmincon','Display', 'off', 'algorithm', 'interior-point', 'HessianApproximation','bfgs','FiniteDifferenceType', 'central','OptimalityTolerance',1e-10,'UseParallel',true);   % different part
        options = optimoptions('ga','InitialPopulationRange',[S.prob.min;S.prob.max ],'PopulationSize',200,'UseParallel',true,'HybridFcn',{'fmincon',hybridopts});
        S.acqui.x = ga(@(x) -Acqfn(S),S.prob.S.prob.dim,[],[],[],[],S.prob.min*ones(1,S.prob.dim),S.prob.max*ones(1,S.prob.dim),[],options);

    case 'multi_start'   % only difference is method to find mininum value
        options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central');
        problem = createOptimProblem('fmincon','objective',@(x) -Acqfn(S),'x0',zeros(S.prob.dim,1),'lb',S.prob.min,'ub',S.prob.max,'options',options);
        S.acqui.x = run(MultiStart('UseParallel',true),problem,500);
end

end