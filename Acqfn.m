function Result = Acqfn(x,S)

S.Hypopt.r = zeros(size(S.add.domain,2),1,S.prob.dim);
x = x';

for i = 1:S.prob.dim
    Mxx = x(i)^2;
    Kxx = (S.add.domain(i,:).^ 2);
    Nxy = x(i)'*S.add.domain(i,:);

    S.Hypopt.r(:,:,i) = Mxx + Kxx' - 2*Nxy';
end

switch S.acqui.mode
    case "PI"
        % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
        a = ones(size(S.add.domain,2),1);

        r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(1,:), [3,1,2])), 3));
        pred_mean = S.Hypopt.alpha(1) + r_r'*S.Hypopt.invR*(S.add.domainy-S.Hypopt.alpha(1));
        mu = r_r'*S.Hypopt.invR*a - 1;
        pred_sig = S.Hypopt.sigma(1)*(1 - r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);

        if S.prob.surconst == true
            pred_constraint_mean = zeros(S.prob.numconstraint,1);
            pred_constraint_sig = zeros(S.prob.numconstraint,1);
            for i = 1 : S.prob.numconstraint
                r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(i+1,:), [3,1,2])), 3));
                pred_constraint_mean(i) = S.Hypopt.alpha(i+1,1) + r_r'*S.Hypopt.invR*(S.add.constdomainy(:,i)-S.Hypopt.alpha(i+1,1));
                mu = r_r'*S.Hypopt.invR*a - 1;
                pred_constraint_sig(i) = S.Hypopt.sigma(i+1,1)*(1-r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);
            end
        end

        u = (S.acqui.minobj - pred_mean);
        ss = sqrt(pred_sig);

        if ss > 0
            Result = normcdf(u/ss);   % PI calculation process
            if S.prob.surconst == 1
                for i = 1:S.prob.numconstraint
                    Result = Result*normcdf((0-pred_constraint_mean(i))/pred_constraint_sig(end,i+1));
                end
            end
        else
            Result = 0;
        end

    case "EI"
        % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
        a = ones(size(S.add.domain,2),1);

        r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(1,:), [3,1,2])), 3));
        pred_mean = S.Hypopt.alpha(1) + r_r'*S.Hypopt.invR*(S.add.domainy-S.Hypopt.alpha(1));
        mu = r_r'*S.Hypopt.invR*a - 1;
        pred_sig = S.Hypopt.sigma(1)*(1 - r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);

        if S.prob.surconst == true
            pred_constraint_mean = zeros(S.prob.numconstraint,1);
            pred_constraint_sig = zeros(S.prob.numconstraint,1);
            for i = 1 : S.prob.numconstraint
                r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(i+1,:), [3,1,2])), 3));
                pred_constraint_mean(i) = S.Hypopt.alpha(i+1,1) + r_r'*S.Hypopt.invR*(S.add.constdomainy(:,i)-S.Hypopt.alpha(i+1,1));
                mu = r_r'*S.Hypopt.invR*a - 1;
                pred_constraint_sig(i) = S.Hypopt.sigma(i+1,1)*(1-r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);
            end
        end

        u = (S.acqui.minobj - pred_mean);
        ss = sqrt(pred_sig);

        p0 = 10^(floor(1+log10(abs(u))));
        if ss > 0
            Result = ((u-p0*S.acqui.exploratio)*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
            if S.prob.surconst ==1
                for i = 1:S.prob.numconstraint
                    Result = Result*normcdf((S.add.surconst_standard(i)-pred_constraint_mean(i))/pred_constraint_sig(i));
                end
            end
        else
            Result = 0;
        end

    case "LCB"

    case "UCB"

    case "MP"  
        if S.prob.surconst == 1
            j = S.prob.numconstraint+1;
        else
            j = 1;
        end
        pred_mean = zeros(j,1);
        for i = 1 : j
            r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(i,:), [3,1,2])), 3));
            if i == 1
                Value = S.add.domainy;
            else
                Value = S.add.constdomainy(:,i-1);
            end
            pred_mean(i,1) = S.Hypopt.alpha(i,1) + r_r'*S.Hypopt.invR*(Value-S.Hypopt.alpha(i,1));
        end
        if S.prob.surconst == 1
            vio_ind = pred_mean(2:S.prob.numconstraint+1)>0; % Find which constraints are violated with the design variable
            Result = pred_mean(1,1) + sum(abs(pred_mean(2:S.prob.numconstraint+1)).*vio_ind);
        else
            Result = pred_mean(1,1);
        end

    case "MP+EI"
        if mod(S.add.cnt,2)==0
            a = ones(size(S.add.domain,2),1);

            r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(1,:), [3,1,2])), 3));
            pred_mean = S.Hypopt.alpha(1) + r_r'*S.Hypopt.invR*(S.add.domainy-S.Hypopt.alpha(1));
            mu = r_r'*S.Hypopt.invR*a - 1;
            pred_sig = S.Hypopt.sigma(1)*(1 - r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);

            if S.prob.surconst == true
                pred_constraint_mean = zeros(S.prob.numconstraint,1);
                pred_constraint_sig = zeros(S.prob.numconstraint,1);
                for i = 1 : S.prob.numconstraint
                    r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(i+1,:), [3,1,2])), 3));
                    pred_constraint_mean(i) = S.Hypopt.alpha(i+1,1) + r_r'*S.Hypopt.invR*(S.add.constdomainy(:,i)-S.Hypopt.alpha(i+1,1));
                    mu = r_r'*S.Hypopt.invR*a - 1;
                    pred_constraint_sig(i) = S.Hypopt.sigma(i+1,1)*(1-r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);
                end
            end

            u = (S.acqui.minobj - pred_mean);
            ss = sqrt(pred_sig);

            p0 = 10^(floor(1+log10(abs(u))));
            if ss > 0
                Result = ((u-p0*S.acqui.exploratio)*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
                if S.prob.surconst == 1
                    for i = 1:S.prob.numconstraint
                        Result = Result*normcdf((0-pred_constraint_mean(i))/pred_constraint_sig(i));
                    end
                end
            else
                Result = 0;
            end
        else
            if S.prob.surconst == 1
                j = S.prob.numconstraint+1;
            else
                j = 1;
            end
            pred_mean = zeros(j,1);
            for i = 1 : j
                r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(i,:), [3,1,2])), 3));
                if i == 1
                    Value = S.add.domainy;
                else
                    Value = S.add.constdomainy(:,i-1);
                end
                pred_mean(i,1) = S.Hypopt.alpha(i,1) + r_r'*S.Hypopt.invR*(Value-S.Hypopt.alpha(i,1));
            end
            if S.prob.surconst == 1
                vio_ind = pred_mean(2:S.prob.numconstraint+1)>0; % Find which constraints are violated with the design variable
                Result = pred_mean(1,1) + sum(abs(pred_mean(2:S.prob.numconstraint+1)).*vio_ind);
            else
                Result = pred_mean(1,1);
            end
        end
end