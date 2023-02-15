function [theta,alpha_kriging,sigma,inv_R,R] = optimizeHypes(Initial_theta, theta, Domain, Domain_y, R_before, r, Iteration, divider, solvertype) % optimize hyperparameters based on MLE

    nSample = size(Domain,2);   % # of the Samples 
    dim = size(Domain,1);       % dim of inputs  
    if size(Initial_theta) ~= dim
        error("dims of hyperparameter and sample input do not match");
    end
    
    R = zeros(nSample,nSample,dim);
    if Iteration == 0
        for i = 1:dim
            Kxy = Domain(i,:)'*Domain(i,:);
            Kxx = repmat(diag(Kxy), 1, nSample);
            R(:,:,i) = Kxx + Kxx' - 2*Kxy;
        end
    else
        for j = 1:dim
            R(:,:,j) = [R_before(:,:,j) r(:,:,j); r(:,:,j)' 0];
        end
    end
    
    % 2. SOLVER SELECTIONS (NOT DETERMINED YET)
    if or(mod(Iteration,divider)==0,Iteration<100*dim)
        switch (lower(solvertype))  % all solvertypes are to find minimum of constrained nonlinear multivariable function
            case 'fmincon'
                options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
                [theta] = fmincon(@(x) -subG_MLE(x),Initial_theta,[],[],[],[],ones(size(Domain,1),1)*0,ones(size(Domain,1),1)*100,[],options);

            case 'fminunc'
                options = optimoptions(@fminunc,'Display','off','algorithm','quasi-newton');
                [theta] = fminunc(@(x) -subG_MLE(x),Initial_theta,options);

            case 'ga'
                options = optimoptions('ga','PopulationSize',300,'UseParallel',true);
                [theta] = ga(@(x) -subG_MLE(x), length(Initial_theta),[],[],[],[],0,10,[],[],options);

            case 'particleswarm'
                options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
                [theta] = particleswarm(@(x) -subG_MLE(x),length(Initial_theta),0,10,options);

            otherwise
                error("solver type is not specific");
        end
    else
        subG_MLE(theta);
    end

    function out_GA = subG_MLE(theta) % make correlation matrix R with exponential kernel(KernelExponential)
        C_xx = exp(-sum(bsxfun(@times, R, permute(theta, [3,1,2])), 3));
        
        % Add the nugget term only when the eigenvalue of correlation matrix is too small to inverse the correlation matrix
        ew = eig(C_xx);
        nugget = 1e-6*eye(size(C_xx,1));       
        if sum(abs(ew)<1e-10)>= 1
            C_xx = C_xx + nugget;
        end

        % C_xx is now reversible
        R_chol = chol(C_xx);
        inv_R = R_chol\(R_chol'\eye(size(R_chol,1)));
        Xi = ones(nSample,1);

        alpha_kriging = (Xi'*inv_R*Domain_y)/(Xi'*inv_R*Xi);
        sigma = 1/nSample*(Domain_y-alpha_kriging*Xi)'*inv_R*(Domain_y-alpha_kriging*Xi);  

        % MLE
        out_GA = -0.5*(nSample*log(sigma) + log(det(C_xx))); 
    end
end
