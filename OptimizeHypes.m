function S = OptimizeHypes(S) % optimize hyperparameters based on MLE

    nSample = size(S.add.domain,2);   % # of the Samples 
    if size(S.Hypopt.initheta) ~= S.prob.dim
        error("dims of hyperparameter and sample input do not match");
    end
    
    if S.add.cnt == 0
        for i = 1:S.prob.dim
            Kxy = S.add.domain(i,:)'*S.add.domain(i,:);
            Kxx = repmat(diag(Kxy), 1, nSample);
            S.Hypopt.R(:,:,i) = Kxx + Kxx' - 2*Kxy;
        end
    else
        S.Hypopt.R(:,end+1,:) = S.Hypopt.r;
        S.Hypopt.R(end+1,:,:) = [permute(S.Hypopt.r,[2,1,3]) zeros(1,1,S.prob.dim)];
    end
    
    % 2. SOLVER SELECTIONS (NOT DETERMINED YET)
    if or(mod(S.add.cnt,S.Hypopt.frequency)==0,S.add.cnt<100*S.prob.dim)
        switch (lower(S.Hypopt.solver))  % all solvertypes are to find minimum of constrained nonlinear multivariable function
            case 'fmincon'
                options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
                [S.Hypopt.theta] = fmincon(@(x) -MLE(x), S.Hypopt.initheta,[],[],[],[],ones(size(S.add.domain,1),1)*0,ones(size(S.add.domain,1),1)*100,[],options);

            case 'fminunc'
                options = optimoptions(@fminunc,'Display','off','algorithm','quasi-newton');
                [S.Hypopt.theta] = fminunc(@(x) -MLE(x),S.Hypopt.initheta,options);

            case 'ga'
                options = optimoptions('ga','PopulationSize',300,'UseParallel',true);
                [S.Hypopt.theta] = ga(@(x) -MLE(x), length(S.Hypopt.initheta),[],[],[],[],ones(size(S.add.domain,1),1)*0,ones(size(S.add.domain,1),1)*100,[],[],options);

            case 'particleswarm'
                options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
                [S.Hypopt.theta] = particleswarm(@(x) -MLE(x),length(S.Hypopt.initheta),ones(size(S.add.domain,1),1)*0,ones(size(S.add.domain,1),1)*100,options);

            otherwise
                error("solver type is not specific");
        end
    else
        MLE(S.Hypopt.theta);
    end


    function MLE_result = MLE(theta) % make correlation matrix R with exponential kernel(KernelExponential)
        C_xx = exp(-sum(bsxfun(@times, S.Hypopt.R, permute(theta, [3,1,2])), 3));
        
        % Add the nugget term only when the eigenvalue of correlation matrix is too small to inverse the correlation matrix
        ew = eig(C_xx);
        nugget = 1e-6*eye(size(C_xx,1));       
        if sum(abs(ew)<1e-10)>= 1
            C_xx = C_xx + nugget;
        end

        % C_xx is now reversible
        R_chol = chol(C_xx);
        S.Hypopt.invR = R_chol\(R_chol'\eye(size(R_chol,1)));
        Xi = ones(nSample,1);

        S.Hypopt.alpha = (Xi'*S.Hypopt.invR*S.add.domainy)/(Xi'*S.Hypopt.invR*Xi);
        S.Hypopt.sigma = 1/nSample*(S.add.domainy-S.Hypopt.alpha*Xi)'*S.Hypopt.invR*(S.add.domainy-S.Hypopt.alpha*Xi);  

        % MLE
        MLE_result = -0.5*(nSample*log(S.Hypopt.sigma) + log(det(C_xx))); 
    end
end
