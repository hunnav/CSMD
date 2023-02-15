%% 0.Setting for Matlab

clear; clc;
format("shortG")
rng("shuffle")
delete(gcp('nocreate'))        % returns the current pool if one exists, otherwise pool will be empty & delete it
parpool('threads')             % creates and returns a thread-based pool.

%% 1.Setting for Bayesopt

[S] = Input_struc;

%% 2.Iteration for Bayesopt

while S.add.cnt < S.prob.maxiter       % until to be maximum iterations
    tic

    % Hyperparameter optimization with new samples based on MLE (maximum likelyhood estimation)
    S = OptimizeHypes(S);

    % Acquisition function for the new point
    S = New_point(S);

    % Scaling and add a new point
    S = Add_point(S);

    % Print
    S.add.cnt
    S.add.minimum_Value(end,1:2)

    toc
end

%% 3.Result

disp(['Minimum_Value_x', num2str(S.add.domain(:,S.add.objective==S.add.minimum_Value(end,2))')])
disp(['Minimum_Value', S.add.minimum_Value(end,2)])