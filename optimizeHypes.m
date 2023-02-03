function [hyp,alpha,sigmasq,invC] = optimizeHypes(theta, x_sample, y_sample, solvertype) % optimize hyperparameters based on MLE

    if nargin < 4   % default solvertype
       solvertype = 'fmincon';
    end

    function out_GA = subG_MLE(theta) % make correlation matrix R with exponential kernel(KernelExponential)
        nSample = size(x_sample,2);   % # of the Samples 
        dim = size(x_sample,1);       % dim of inputs
        if size(theta) ~= dim    
            error("dims of hyperparameter and sample input do not match");
        end
        
        Ath = diag(theta);  
        Kxy = x_sample'*Ath*x_sample;
        Kxx = repmat(diag(Kxy), 1, nSample);
        corr = Kxx + Kxx' - 2*Kxy;
        C_xx = exp(-corr);

        % Add the nugget term only when the eigenvalue of correlation matrix is too small to inverse the correlation matrix
        ew = eig(C_xx);
        nugget = 1e-4*eye(size(C_xx,1));
        for i = 1: length(ew) 
            if(abs(ew(i))<1e-10)
                C_xx = C_xx + nugget;
            end
        end
        % C_xx is now reversible
        R_chol = chol(C_xx);
        invC = R_chol\(R_chol'\eye(size(R_chol,1)));
        Xi = ones(nSample,1);

        alpha = (Xi'*invC*y_sample)/(Xi'*invC*Xi);
        sigmasq = 1/nSample*(y_sample-alpha*Xi)'*invC*(y_sample-alpha*Xi);  

        % MLE
        out_GA = -0.5*(nSample*log(sigmasq) + log(det(C_xx))); 
    end

    f = @(x) -subG_MLE(x);   % to find maximum, we use '-'sign

    % 2. SOLVER SELECTIONS (NOT DETERMINED YET)
    switch (lower(solvertype))  % all solvertypes are to find minimum of constrained nonlinear multivariable function
        case 'fmincon'
            options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
            [hyp] = fmincon(f,theta,[],[],[],[],ones(size(x_sample,1),1)*0,ones(size(x_sample,1),1)*100,[],options); 

        case 'fminunc'
            options = optimoptions(@fminunc,'Display','off','algorithm','quasi-newton');
            [hyp] = fminunc(f,theta,options);
        
        case 'ga'
            options = optimoptions('ga','PopulationSize',300,'UseParallel',true);
            [hyp] = ga(f, length(theta),[],[],[],[],0,1e5,[],[],options);

        case 'particleswarm'
            options = optimoptions('particleswarm','UseParallel',true,'SwarmSize',200);
            [hyp] = particleswarm(f,length(theta),0,10,options);
                  
        otherwise
            error("solver type is not specifie");
    end
end
