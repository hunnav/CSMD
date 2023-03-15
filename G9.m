%% 0.Setting for Matlab

clear; close; clc;
rng("shuffle")
delete(gcp('nocreate'))
parpool('threads')

for iter = 1:10
    %% 1.Setting for Bayesopt

    [S] = Input_struc;

    %% 2.Iteration for Bayesopt

    while S.add.cnt < S.prob.maxiter
        tic

        % Hyperparameter optimization based on MLE
        S = OptimizeHypes(S);

        % Acquisition function for the new point
        S = New_point(S);

        % Scaling and add a new point
        S = Add_point(S);

        % Print
        S = Print(S);

        toc
    end

    %% 3.Result

    fprintf('\n Minimum_Value_x : %s',num2str(S.add.domain(:,S.add.objective==S.add.minimum_Value(end,2))'));
    fprintf('\n Minimum_Value : %g\n\n',S.add.minimum_Value(end,2));
    save(sprintf('%s_%d', S.prob.filename, iter), 'S');

end

% Quantile Transformation
% bayesian optimization with mixed integer