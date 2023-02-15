options = optimoptions(@fmincon,'Display', 'off', 'algorithm', 'interior-point','HessianApproximation','bfgs','FiniteDifferenceType', 'central','UseParallel',true);
[hyp] = fmincon(@(x) -acq_MP(x,S),Initial_theta,[],[],[],[],ones(size(x_sample,1),1)*0,ones(size(x_sample,1),1)*100,[],options);

function [pred_mean] = acq_MP(x,S)

    % correlation r (for speed)
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
end
