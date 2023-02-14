clear;
clc;
format("shortG")
rng("shuffle")
delete(gcp('nocreate'))        % returns the current pool if one exists, otherwise pool will be empty & delete it
parpool('threads')             % creates and returns a thread-based pool.

%% 1.Setting for BayesOpt

[S] = Input_struc;

Constraint = zeros(max_iter,1);
Add = zeros(num_initial_value,1);
Standard = zeros(max_iter,6);
A = -min(Objective)+1; 
D = exp(log(K+1)/(K/2)); 
Domain_y = Objective+A;
Minimum_Value = [0,inf];

%% 2.Iteration for BayesOpt

% Infilling criterion : EI process
while Iteration < max_iter       % until to be maximum iterations
    tic

    % Hyperparameter optimization with new samples based on MLE (maximum likelyhood estimation)
    [theta,alpha_kriging,sigma,inv_R,R] = optimizeHypes(Initial_theta, theta, Domain, Domain_y, R, r, Iteration, divider, MLE_mode);
    
    % EI process to extract the new point(Dom_EI)
    if  Minimum_Value(end,2) == inf
        min_obj = min(Domain_y);
    else
        min_obj = Domain_y(num_initial_value + Minimum_Value(end,1));
    end
    ratio_or_weight_modify = ratio_or_weight;
    while 1
        while 1
            try
                new_x = EIval_check(Domain, Domain_y, theta, sigma, alpha_kriging, inv_R, low_Range, upper_Range, min_obj, EI_mode, EI_acq_mode, ratio_or_weight_modify);
                break
            catch
                disp('There is some error, repeat again.')
            end
        end
        r = Correlation(Domain,new_x);
        if min(sum(r, 3)) > beta
            break
        else
            if abs(ratio_or_weight_modify) >= 5
                break
            else
                if strcmp(EI_acq_mode,'normal')
                    ratio_or_weight_modify = ratio_or_weight_modify+1;
                else
                    ratio_or_weight_modify = ratio_or_weight_modify-1;
                end
            end
        end
    end
    
    Iteration = Iteration + 1    % add the number of iteration
    n = size(Domain,2)+1;
    Domain(:,n) = new_x; 
    Objective(n,1) = f(new_x(1),new_x(2),new_x(3),new_x(4),new_x(5),new_x(6),new_x(7));    % we need to change it according to the conditon

    gg1 = 0-g1(new_x(1),new_x(2),new_x(3),new_x(4),new_x(5),new_x(6),new_x(7));            % we need to change it according to the conditon
    gg2 = 0-g2(new_x(1),new_x(2),new_x(3),new_x(4),new_x(5),new_x(6),new_x(7));            % we need to change it according to the conditon
    gg3 = 0-g3(new_x(1),new_x(2),new_x(3),new_x(4),new_x(5),new_x(6),new_x(7));            % we need to change it according to the conditon
    gg4 = 0-g4(new_x(1),new_x(2),new_x(3),new_x(4),new_x(5),new_x(6),new_x(7));            % we need to change it according to the conditon
    Constraint(n,1) = max([gg1,gg2,gg3,gg4,0]);
    
    [Standard,Domain_y,Modified_Objective,Add] = Scaling(Standard,Iteration,Objective,Constraint,K);

    if (Add(end) == 0 && Objective(end) < Minimum_Value(end,2))
        Minimum_Value(end+1,:) = [Iteration, Objective(end)];
    end
    Minimum_Value(end,1:2)

    toc
end

%% 3.Result

Minimum_Value_x = Domain(:,Objective==Current_Minimum_Value)'
Current_Minimum_Value