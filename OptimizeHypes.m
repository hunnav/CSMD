function S = OptimizeHypes(S) % Optimize hyperparameters based on MLE

if size(S.Hypopt.initheta) ~= S.prob.dim
    error("Dimension of hyperparameter and sample do not match");
end

if S.add.cnt == 0
    for i = 1:S.prob.dim
        Kxy = S.add.domain(i,:)'*S.add.domain(i,:);
        Kxx = repmat(diag(Kxy), 1, size(S.add.domain,2));
        S.Hypopt.R(:,:,i) = Kxx + Kxx' - 2*Kxy;
    end
else
    S.Hypopt.R(:,end+1,:) = S.Hypopt.r;
    S.Hypopt.R(end+1,:,:) = [permute(S.Hypopt.r,[2,1,3]) zeros(1,1,S.prob.dim)];
end

% SOLVER SELECTIONS
if S.prob.surconst == true
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
        if S.Hypopt.validation == true
            initheta = zeros(S.Hypopt.validation_trynum,S.prob.dim);
            initheta(1,:) = zeros(1,S.prob.dim);
            initheta(2,:) = S.Hypopt.theta(i,:);
            initheta(3:end-1,:) = (S.Hypopt.max-S.Hypopt.min).*lhsdesign(S.Hypopt.validation_trynum-3,S.prob.dim) + (S.Hypopt.min);
            initheta(end,:) = zeros(1,S.prob.dim);
            num_domainy = size(Value,1);
            training_num = round(num_domainy*S.Hypopt.validation_ratio);
            test_num = num_domainy-training_num;
            except_num = sort(randperm(num_domainy,test_num));

            idx = true(size(S.add.domain));
            idx(:,except_num) = false;
            training_domain = reshape(S.add.domain(idx), S.prob.dim, training_num);
            idx = ~idx;
            validation_domain = reshape(S.add.domain(idx), S.prob.dim, test_num);

            idx = true(size(Value));
            idx(except_num) = false;
            training_domainy = Value(idx);
            idx = ~idx;
            validation_domainy = Value(idx);

            idx = true(size(S.Hypopt.R));
            idx(except_num,:,:) = false;
            idx(:,except_num,:) = false;
            training_R = reshape(S.Hypopt.R(idx), training_num, training_num, S.prob.dim);

            domainy = training_domainy;
            R = training_R;
            nSample = training_num;
            k = 0;
        else
            domainy = Value;
            R = S.Hypopt.R;
            nSample = size(S.add.domain,2);
        end

        while 1
            if S.Hypopt.validation == true
                k = k+1;
                S.Hypopt.initheta = initheta(k,:);
            end
            while 1
                try
                    switch (lower(S.Hypopt.solver))
                        case 'fmincon'
                            options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
                            S.Hypopt.theta(i,:) = fmincon(@(x) -MLE(x), S.Hypopt.initheta,[],[],[],[],ones(S.prob.dim,1)*S.Hypopt.min,ones(S.prob.dim,1)*S.Hypopt.max,[],options);

                        case 'fminunc'
                            options = optimoptions(@fminunc,'Display','off','algorithm','quasi-newton','UseParallel',true);
                            S.Hypopt.theta(i,:) = fminunc(@(x) -MLE(x), S.Hypopt.initheta,options);

                        case 'ga'
                            options = optimoptions('ga','InitialPopulationRange',[S.Hypopt.min;S.Hypopt.max],'PopulationSize',1000,'UseParallel',true);
                            S.Hypopt.theta(i,:) = ga(@(x) -MLE(x), S.prob.dim,[],[],[],[],ones(S.prob.dim,1)*S.Hypopt.min,ones(S.prob.dim,1)*S.Hypopt.max,[],options);

                        case 'ga+fmincon'
                            hybridopts = optimoptions('fmincon','Display', 'off', 'algorithm', 'interior-point', 'HessianApproximation','bfgs','FiniteDifferenceType', 'central','OptimalityTolerance',1e-10,'UseParallel',true);
                            options = optimoptions('ga','InitialPopulationRange',[S.Hypopt.min;S.Hypopt.max],'PopulationSize',1000,'UseParallel',true,'HybridFcn',{'fmincon',hybridopts});
                            S.Hypopt.theta(i,:) = ga(@(x) -MLE(x), S.prob.dim,[],[],[],[],ones(S.prob.dim,1)*S.Hypopt.min,ones(S.prob.dim,1)*S.Hypopt.max,[],options);

                        case 'pso'
                            options = optimoptions('particleswarm','SwarmSize',200,'UseParallel',true);
                            S.Hypopt.theta(i,:) = particleswarm(@(x) -MLE(x), S.prob.dim,ones(S.prob.dim,1)*S.Hypopt.min,ones(S.prob.dim,1)*S.Hypopt.max,options);

                        case 'dace'
                            if strcmp(S.Hypopt.dace_reg, 'regpolyauto')
                                if nSample >= (S.prob.dim+1)*(S.prob.dim+2)*(S.prob.dim+3)/6
                                    regmode = 'regpoly3';
                                elseif nSample >= (S.prob.dim+1)*(S.prob.dim+2)/2
                                    regmode = 'regpoly2';
                                elseif nSample >= 1+S.prob.dim
                                    regmode = 'regpoly1';
                                else
                                    regmode = 'regpoly0';
                                end
                            else
                                regmode = S.Hypopt.dace_reg;
                            end
                            if S.prob.surconst == true
                                if strcmp(S.acqui.const_mode, 'same')
                                    dace_domain = S.add.domain;
                                else
                                    if i == 1
                                        dace_domain = S.add.domain;
                                    else
                                        if strcmp(S.acqui.const_mode, 'different')
                                            dace_domain(:,1:S.prob.numinitsam) = S.add.iniconstdomain(:,1:S.prob.numinitsam);
                                        else
                                            dace_domain(:,1:S.prob.numinitsam) = S.add.iniconstdomain(:,1:S.prob.numinitsam,i-1);
                                        end
                                    end
                                end
                            else
                                dace_domain = S.add.domain;
                            end
                            if S.Hypopt.validation == true
                                idx = true(size(dace_domain));
                                idx(:,except_num) = false;
                                dace_domain = reshape(dace_domain(idx), S.prob.dim, training_num);
                            end
                            dmodel = dacefit(dace_domain', domainy, regmode, S.Hypopt.dace_cor, S.Hypopt.initheta, (S.Hypopt.min+1e-10)*ones(1,S.prob.dim), S.Hypopt.max*ones(1,S.prob.dim));
                            S.Hypopt.theta(i,:) = dmodel.theta;
                            MLE(S.Hypopt.theta(i,:));

                        otherwise
                            error("solver type is not specific");
                    end
                    break
                catch
                    disp('There is some error, repeat again.')
                end
            end

            if S.Hypopt.validation == true
                S.Hypopt.r = zeros(training_num,1,S.prob.dim);
                result = zeros(test_num,1);
                for test_sequence = 1:test_num
                    x = validation_domain(:,test_sequence);
                    x = x';
                    for d = 1:S.prob.dim
                        Mxx = x(d)^2;
                        Kxx = (training_domain(d,:).^ 2);
                        Nxy = x(d)'*training_domain(d,:);

                        S.Hypopt.r(:,:,d) = Mxx + Kxx' - 2*Nxy';
                    end
                    r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(i,:), [3,1,2])), 3));
                    a = ones(size(training_domain,2),1);
                    pred_mean = S.Hypopt.alpha(i) + r_r'*S.Hypopt.invR*(training_domainy-S.Hypopt.alpha(i));
                    mu = r_r'*S.Hypopt.invR*a - 1;
                    pred_sig = sqrt(S.Hypopt.sigma(i)*(1 - r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu));
                    if and(validation_domainy(test_sequence)>=pred_mean-pred_sig*1.96, validation_domainy(test_sequence)<=pred_mean+pred_sig*1.96)
                        result(test_sequence,1) = 0;
                    else
                        fprintf('Use new initial theta (%g)\n',k);
                        result(test_sequence,1) = 1;
                        break
                    end
                end
                if or(sum(result) == 0, k == S.Hypopt.validation_trynum)
                    domainy = Value;
                    R = S.Hypopt.R;
                    nSample = size(S.add.domain,2);
                    MLE(S.Hypopt.theta(i,:));
                    break
                end
            else
                break
            end
        end
    else
        domainy = Value;
        R = S.Hypopt.R;
        nSample = size(S.add.domain,2);
        MLE(S.Hypopt.theta(i,:));
    end
end

    function MLE_result = MLE(theta) % Exponential kernel
        C_xx = exp(-sum(bsxfun(@times, R, permute(theta, [3,1,2])), 3));

        % Add the nugget term only when the eigenvalue of correlation matrix is too small to inverse the correlation matrix
        nugget = 1e-6*eye(size(C_xx,1));
        if sum(abs(eig(C_xx))<1e-6)>= 1
            C_xx = C_xx + nugget;
        end

        R_chol = chol(C_xx);
        S.Hypopt.invR = R_chol\(R_chol'\eye(size(R_chol,1)));
        Xi = ones(nSample,1);
        S.Hypopt.alpha(i,1) = (Xi'*S.Hypopt.invR*domainy)/(Xi'*S.Hypopt.invR*Xi);
        S.Hypopt.sigma(i,1) = 1/nSample*(domainy-S.Hypopt.alpha(i,1)*Xi)'*S.Hypopt.invR*(domainy-S.Hypopt.alpha(i,1)*Xi);

        % MLE
        if det(C_xx) == 0
            MLE_result = -0.5*(nSample*log(S.Hypopt.sigma(i,1)));
        else
            MLE_result = -0.5*(nSample*log(S.Hypopt.sigma(i,1)) + log(det(C_xx)));
        end
    end
end
