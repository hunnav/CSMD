%% 0.Setting for Matlab

clear; clc;
rng("shuffle")
delete(gcp('nocreate'))
parpool('threads')             

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
    fprintf('\n Iteration : %g',S.add.cnt);
    fprintf('\n Current minimum value : %g (%g)\n\n',S.add.minimum_Value(end,2),S.add.minimum_Value(end,1));

    toc
end

%% 3.Result

fprintf('\n Minimum_Value_x : %s',num2str(S.add.domain(:,S.add.objective==S.add.minimum_Value(end,2))'));
fprintf('\n Minimum_Value : %g\n\n',S.add.minimum_Value(end,2));