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
        r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(1,:), [3,1,2])), 3));
        % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
        a = ones(size(S.add.domain,2),1);
        pred_mean = S.Hypopt.alpha(1) + r_r'*S.Hypopt.invR*(S.add.domainy-S.Hypopt.alpha(1));
        mu = r_r'*S.Hypopt.invR*a - 1;
        pred_sig = S.Hypopt.sigma(1)*(1 - r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);
        
        u = (S.acqui.minobj - pred_mean);
        ss = sqrt(pred_sig);

        if ss > 0
            Result =  normcdf(u/ss);    % PI calculation process
        else
            Result = 0;
        end

    case "EI"
        r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(1,:), [3,1,2])), 3));
        % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
        a = ones(size(S.add.domain,2),1);
        pred_mean = S.Hypopt.alpha(1) + r_r'*S.Hypopt.invR*(S.add.domainy-S.Hypopt.alpha(1));
        mu = r_r'*S.Hypopt.invR*a - 1;
        pred_sig = S.Hypopt.sigma(1)*(1 - r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);

        u = (S.acqui.minobj - pred_mean);
        ss = sqrt(pred_sig);

        p0 = 10^(floor(1+log10(abs(u))));

        if ss > 0
            Result =  ((u-p0*S.acqui.exploratio)*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
        else
            Result = 0;
        end

    case "LCB"

    case "UCB"

    case "MP"
        pred_mean = zeros(S.prob.numconstraint+1,1);
        for i = 1 : S.prob.numconstraint+1
            r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(i,:), [3,1,2])), 3));
            if i == 1
                Value = S.add.domainy;
            else
                Value = S.add.constdomainy(:,i-1);
            end
            pred_mean(i,1) = S.Hypopt.alpha(i,1) + r_r'*S.Hypopt.invR*(Value-S.Hypopt.alpha(i,1));
        end
        vio_ind = pred_mean(2:S.prob.numconstraint+1)<0; % Find which constraints are violated with the design variable
       
        if isempty(vio_ind) ~= 1
            Result = pred_mean(1,1) + sum(abs(pred_mean(2:S.prob.numconstraint+1)).*vio_ind);
        else
            Result = pred_mean(1,1);
        end
        
    case "MP+EI"
        if mod(S.add.cnt,2)==0
            r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(1,:), [3,1,2])), 3));
            % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
            a = ones(size(S.add.domain,2),1);
            pred_mean = S.Hypopt.alpha(1) + r_r'*S.Hypopt.invR*(S.add.domainy-S.Hypopt.alpha(1));
            mu = r_r'*S.Hypopt.invR*a - 1;
            pred_sig = S.Hypopt.sigma(1)*(1 - r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);

            u = (S.acqui.minobj - pred_mean);
            ss = sqrt(pred_sig);

            p0 = 10^(floor(1+log10(abs(u))));

            if ss > 0
                Result =  ((u-p0*S.acqui.exploratio)*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
            else
                Result = 0;
            end
        else
            pred_mean = zeros(S.prob.numconstraint+1,1);
            for i = 1 : S.prob.numconstraint+1
                r_r = exp(-sum(bsxfun(@times, S.Hypopt.r, permute(S.Hypopt.theta(i,:), [3,1,2])), 3));
                if i == 1
                    Value = S.add.domainy;
                else
                    Value = S.add.constdomainy(:,i-1);
                end
                pred_mean(i,1) = S.Hypopt.alpha(i,1) + r_r'*S.Hypopt.invR*(Value-S.Hypopt.alpha(i,1));
            end
            vio_ind = pred_mean(2:S.prob.numconstraint+1)<0; % Find which constraints are violated with the design variable

            if isempty(vio_ind) ~= 1
                Result = pred_mean(1,1) + sum(abs(pred_mean(2:S.prob.numconstraint+1)).*vio_ind);
            else
                Result = pred_mean(1,1);
            end
        end

        % 그리고 마지막으로 MP+EI process는 직접 code를 짜서 진행해보면 좋을 듯.
        % 일반적으로 EI+MP를 많이 쓰는 것으로 인지하고 있음.
        % 짜기전에 나머지 code 다 정리하고(S로), EI+MP 어떻게 짜야하는지 와서 물어보면 될 듯.

    otherwise
end