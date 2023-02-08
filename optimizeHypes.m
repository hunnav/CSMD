function [hyp,alpha,sigmasq,invC,R] = optimizeHypes(Initial_theta, theta, x_sample, y_sample, R_before, r, miniter, divider, solvertype) % optimize hyperparameters based on MLE

    nSample = size(x_sample,2);   % # of the Samples 
    dim = size(x_sample,1);       % dim of inputs  
    if size(Initial_theta) ~= dim
        error("dims of hyperparameter and sample input do not match");
    end
    
    R = zeros(nSample,nSample,dim);
    if miniter == 0
        for i = 1:dim
            Kxy = x_sample(i,:)'*x_sample(i,:);
            Kxx = repmat(diag(Kxy), 1, nSample);
            R(:,:,i) = Kxx + Kxx' - 2*Kxy;
        end
    else
        for j = 1:dim
            R(:,:,j) = [R_before(:,:,j) r(:,:,j); r(:,:,j)' 0];
        end
    end
    
    % 2. SOLVER SELECTIONS (NOT DETERMINED YET)
    if or(rem(miniter,divider)==0,miniter<100*dim)
        switch (lower(solvertype))  % all solvertypes are to find minimum of constrained nonlinear multivariable function
            case 'fmincon'
                options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
                [hyp] = fmincon(@(x) -subG_MLE(x),Initial_theta,[],[],[],[],ones(size(x_sample,1),1)*0,ones(size(x_sample,1),1)*100,[],options);

            case 'fminunc'
                options = optimoptions(@fminunc,'Display','off','algorithm','quasi-newton');
                [hyp] = fminunc(@(x) -subG_MLE(x),Initial_theta,options);

            case 'ga'
                options = optimoptions('ga','PopulationSize',300,'UseParallel',true);
                [hyp] = ga(@(x) -subG_MLE(x), length(Initial_theta),[],[],[],[],0,10,[],[],options);

            case 'particleswarm'
                options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
                [hyp] = particleswarm(@(x) -subG_MLE(x),length(Initial_theta),0,10,options);

            otherwise
                error("solver type is not specific");
        end
    else
        hyp = theta;
        subG_MLE(theta);
    end

    function out_GA = subG_MLE(theta) % make correlation matrix R with exponential kernel(KernelExponential)
        C_xx = exp(-sum(bsxfun(@times, R, permute(theta, [3,1,2])), 3));
        
        % Add the nugget term only when the eigenvalue of correlation matrix is too small to inverse the correlation matrix
        ew = eig(C_xx);
        nugget = 1e-6*eye(size(C_xx,1));       
        C_xx = C_xx + nugget * sum(abs(ew)<1e-10);

        % C_xx is now reversible
        R_chol = chol(C_xx);
        invC = R_chol\(R_chol'\eye(size(R_chol,1)));
        Xi = ones(nSample,1);

        alpha = (Xi'*invC*y_sample)/(Xi'*invC*Xi);
        sigmasq = 1/nSample*(y_sample-alpha*Xi)'*invC*(y_sample-alpha*Xi);  

        % MLE
        out_GA = -0.5*(nSample*log(sigmasq) + log(det(C_xx))); 
    end
end
