function [S] = acquifcn(S)
% 각 parts 별로 PI,EI,LCB,UCB,... function을 생성하여 만들고
% 그 function을 이용해서 최대한 짧게 정리하면 좋을 듯.
% 각 항목별로 너무 길이가 길어지면 시각적으로 한번에 어떤 methodology를 이용할 수 있는지 확인이 불가.
% 또한 acq_solver라는 code를 통해 코드를 더 짧게 만들면 좋을 듯.
% (파일 내 내용은 수정하지 않았으나,
% acq_solver 생성)

switch S.acqui.mode
    case "PI"

    case "EI"

        nSample = size(initial_Dom,2);   % # of the Samples
        dim = size(initial_Dom,1);       % dim of inputs

        % Gaussian Process Regression
        r = zeros(nSample,1,dim);
        new_x = new_x';

        for i = 1:dim
            Mxx = new_x(i)^2;
            Kxx = (initial_Dom(i,:).^ 2);
            Nxy = new_x(i)'*initial_Dom(i,:);

            r(:,:,i) = Mxx + Kxx' - 2*Nxy';
        end

        r_r = exp(-sum(bsxfun(@times, r, permute(theta, [3,1,2])), 3));

        % For ordinary kriging (if the value of 1 consist of function f(x), it is universal kriging)
        a = ones(num_samp,1);
        pred_mean = alpha + r_r'*inv_R*(initial_Domy-alpha);
        mu = r_r'*inv_R*a - 1;
        pred_sig = sigma*(1 - r_r'*inv_R*r_r + mu'*(a'*inv_R*a)\mu);

        % Based on the GPR model, EI process is progressed.
        u = (min_obj - pred_mean);
        ss = sqrt(pred_sig);

        p0 = 10^(floor(1+log10(abs(u))));

        if ss > 0
            EI_p =  ((u-p0*ratio_or_weight)*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
        else
            EI_p = 0;
        end




    case "LCB"

    case "UCB"

    case "MP"
        
    case "EIMP"

    otherwise
end