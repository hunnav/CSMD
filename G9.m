%% 0.Setting for Matlab

clear; clc;
rng("shuffle")
delete(gcp('nocreate'))        % returns the current pool if one exists, otherwise pool will be empty & delete it
parpool('threads')             % creates and returns a thread-based pool.

%% 1.Setting for Bayesopt

[S] = Input_struc;

%% 2.Iteration for Bayesopt

while Iteration < max_iter       % until to be maximum iterations
    tic

    % Hyperparameter optimization with new samples based on MLE (maximum likelyhood estimation)
    [theta,alpha_kriging,sigma,inv_R,R] = optimizeHypes(Initial_theta, theta, Domain, Domain_y, R, r, Iteration, divider, MLE_mode);
    
    % EI process to extract the new point(Dom_EI)
    [x,r] = EIval_check(Domain, Domain_y, theta, sigma, alpha_kriging, inv_R, low_Range, upper_Range, ratio, beta, num_initial_value, Minimum_Value, EI_mode);

    % Scaling and add a new point
    Minimum_Value = Add_new_point(num_initial_value, Iteration, Domain, Objective, f, g1, g2, g3, g4, )

    % Print
    Iteration
    Minimum_Value(end,1:2)
    
    toc
end

%% 3.Result

Minimum_Value_x = Domain(:,Objective==Minimum_Value(end,2))'
Minimum_Value(end,2)