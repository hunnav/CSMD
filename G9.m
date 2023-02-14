clear;
clc;
format("shortG")
rng("shuffle")
delete(gcp('nocreate'))        % returns the current pool if one exists, otherwise pool will be empty & delete it
parpool('threads')             % creates and returns a thread-based pool.

%% 1.Setting for BayesOpt

[S] = input_struc;

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
                Dom_EI = EIval_check(Domain, Domain_y, theta, sigma, alpha_kriging, inv_R, low_Range, upper_Range, min_obj, EI_mode, EI_acq_mode, ratio_or_weight_modify);
                break
            catch
                disp('There is some error, repeat again.')
            end
        end
        r = Correlation(Domain,Dom_EI);
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
    x = Dom_EI';
    gg1 = 0-g1(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    gg2 = 0-g2(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    gg3 = 0-g3(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    gg4 = 0-g4(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    
    % Add sample from the EI process
    n = size(Domain,2)+1;
    Domain(:,n) = x; 
    Objective(n,1) = f(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
    Constraint(n,1) = max([gg1,gg2,gg3,gg4,0]);
    
    Iteration = Iteration + 1    % add the number of iteration
    Standard(Iteration,1) = Iteration;
    if Objective(n,1)+A < 1
        A = -Objective(n,1)+1;    
    end
    Standard(Iteration,2) = median(Objective+A);
    B = exp(log(Standard(Iteration,2))/K);
    Reverse_Modified_Objective = K - log(Objective+A)/log(B);
    Reverse_Modified_Objective(Reverse_Modified_Objective < 0) = 0;
    Modified_Objective = K-log(Reverse_Modified_Objective+1)/log(D);
    non_zero_numbers = Constraint(Constraint ~= 0);
    if isempty(non_zero_numbers)
        Standard(Iteration,3) = 0;
        C = 5;
    else
        Standard(Iteration,3) = median(non_zero_numbers);
        C = exp(log(Standard(Iteration,3))/(K/2));
        Add = log(Constraint+C)/log(C)-1;
    end
    Domain_y = Modified_Objective + Add;
    Standard(Iteration,4:7) = [A,B,C,D];

    if (Add(end) == 0 && Objective(end) < Minimum_Value(end,2))
        Minimum_Value(end+1,:) = [Iteration, Objective(end)];
    end
    Minimum_Value(end,1:2)

    toc
end

%% 3.Result

Minimum_Value_x = Domain(:,Objective==Current_Minimum_Value)'
Current_Minimum_Value