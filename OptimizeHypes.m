function S = OptimizeHypes(S) % Optimize hyperparameters based on MLE

nSample = size(S.add.domain,2);
if size(S.Hypopt.initheta) ~= S.prob.dim
    error("Dimension of hyperparameter and sample do not match");
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

% SOLVER SELECTIONS
if strcmp(S.acqui.mode,'MP')
    j = S.prob.numconstraint+1;
else
    j = 1;
end
for i = 1:j
    if i == 1
        Value = S.add.domainy;
    else
        Value = S.add.constdomainy(:,i-1);
    end
    if or(mod(S.add.cnt,S.Hypopt.frequency)==0,S.add.cnt<100*S.prob.dim)
        switch (lower(S.Hypopt.solver))
            case 'fmincon'
                options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
                S.Hypopt.theta(i,:) = fmincon(@(x) -MLE(x), S.Hypopt.initheta,[],[],[],[],ones(S.prob.dim,1)*0,ones(S.prob.dim,1)*100,[],options);

            case 'fminunc'
                options = optimoptions(@fminunc,'Display','off','algorithm','quasi-newton','UseParallel',true);
                S.Hypopt.theta(i,:) = fminunc(@(x) -MLE(x), S.Hypopt.initheta,options);

            case 'ga'
                options = optimoptions('ga','PopulationSize',300,'UseParallel',true);
                S.Hypopt.theta(i,:) = ga(@(x) -MLE(x), length(S.Hypopt.initheta),[],[],[],[],ones(S.prob.dim,1)*0,ones(S.prob.dim,1)*100,[],[],options);

            case 'pso'
                options = optimoptions('particleswarm','SwarmSize',200,'UseParallel',true);
                S.Hypopt.theta(i,:) = particleswarm(@(x) -MLE(x), S.prob.dim,ones(S.prob.dim,1)*0,ones(S.prob.dim,1)*100,options);

            otherwise
                error("solver type is not specific");
        end
    else
        MLE(S.Hypopt.theta);
    end
end

    function MLE_result = MLE(theta) % Exponential kernel
        C_xx = exp(-sum(bsxfun(@times, S.Hypopt.R, permute(theta, [3,1,2])), 3));

        % Add the nugget term only when the eigenvalue of correlation matrix is too small to inverse the correlation matrix
        ew = eig(C_xx);
        nugget = 1e-6*eye(size(C_xx,1));
        if sum(abs(ew)<1e-10)>= 1
            C_xx = C_xx + nugget;
        end

        R_chol = chol(C_xx);
        S.Hypopt.invR = R_chol\(R_chol'\eye(size(R_chol,1)));
        Xi = ones(nSample,1);

        S.Hypopt.alpha(i,1) = (Xi'*S.Hypopt.invR*Value)/(Xi'*S.Hypopt.invR*Xi);
        S.Hypopt.sigma(i,1) = 1/nSample*(Value-S.Hypopt.alpha(i,1)*Xi)'*S.Hypopt.invR*(Value-S.Hypopt.alpha(i,1)*Xi);

        % MLE
        MLE_result = -0.5*(nSample*log(S.Hypopt.sigma(i,1)) + log(det(C_xx)));
    end
end
