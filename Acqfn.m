function Result = Acqfn(S)

r_r = exp(-sum(bsxfun(@times, Correlation(S), permute(S.Hypopt.theta, [3,1,2])), 3));

% For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
a = ones(size(S.add.domain,2),1);
pred_mean = S.Hypopt.alpha + r_r'*S.Hypopt.invR*(S.add.domainy-S.Hypopt.alpha);


switch S.acqui.mode
    case "PI"
        % PI evaulation
        mu = r_r'*S.Hypopt.invR*a - 1;
        pred_sig = sigma*(1 - r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);

        % % For simple kriging (remove the last part of the Ordinary kriging)
        % % It can be used only when the mean function is a known deterministic function
        % mean_pred = alpha+r_r'*inv_R*(y_sample-alpha*a);
        % sigma_pred = sigma*(1-r_r'*inv_R*r_r);

        % % unit_test for r_r
        % for i = 1 : num_iter_Samp
        %         rrr(i,1) = exp(-(theta'*(initial_Dom(:,i)-new_x(:)).^2));
        % end

        u = (S.acqui.minobj - pred_mean);
        ss = sqrt(pred_sig);

        if ss > 0
            Result =  normcdf(u/ss);    % PI calculation process
        else
            Result = 0;
        end

    case "EI"
        % EI evaulation
        mu = r_r'*S.Hypopt.invR*a - 1;
        pred_sig = sigma*(1 - r_r'*S.Hypopt.invR*r_r + mu'*(a'*S.Hypopt.invR*a)\mu);
        u = (S.acqui.minobj - pred_mean);
        ss = sqrt(pred_sig);

        p0 = 10^(floor(1+log10(abs(u))));

        if ss > 0
            Result =  ((u-p0*ratio)*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
        else
            Result = 0;
        end

    case "LCB"

    case "UCB"

    case "MP"
        % MP evaulation
        const_group = [g1,g2,g3,g4]; % Group the constraints from the design varaible
        % Make sure that g value should be minus
        % if design variable violate the given constraints
        vio_ind = find(find(const_group < 0)>0); % Find which constraints are violated with the design variable

        if ismepty(vio_ind) ~= 1
            pred_mean = pred_mean + sum(abs([g1,g2,g3,g4]).*vio_ind);
        else
        end
        % 현재 constraints 가 위반 할 시, constraint ( g1,g2,g3,g4) 의 값이 음수값을 갖도록 통일하면
        % 코드가 되게 편해질 듯. ( 각 constraints 마다 매번 처리를 해줘야 하지만,
        % 프리 스로세스로 해줄 수 있다고 개인적으로 생각이 듬.
        % 그리고 마지막으로 MP+EI process는 직접 code를 짜서 진행해보면 좋을 듯.
        % 일반적으로 EI+MP를 많이 쓰는 것으로 인지하고 있음.
        % 짜기전에 나머지 code 다 정리하고(S로), EI+MP 어떻게 짜야하는지 와서 물어보면 될 듯.

    case "EIMP"

    otherwise
end