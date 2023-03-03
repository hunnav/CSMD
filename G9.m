
%% 0.Setting for Matlab

clear; close; clc;
rng("shuffle")
delete(gcp('nocreate'))
parpool('threads')

for iter = 1:15
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

        S.print.x(S.add.cnt) = S.add.cnt;
        S.print.y(S.add.cnt) = norm([2.3331998;1.942837;-0.479546;4.387811;-0.632617;1.039775;1.600996]-[S.add.domain(:,S.add.objective==S.add.minimum_Value(end,2))],'inf');   % Also can use pdist2
        addpoints(animatedline,S.print.x,S.print.y); drawnow;

        if S.print.y(S.add.cnt) < S.print.near(end,end) && S.add.constraint(end) == 0
            S.print.near(end+1,:) = [S.add.domain(:,end)', S.add.objective(end),S.print.y(S.add.cnt)];
        end

        toc
    end

    %% 3.Result

    fprintf('\n Minimum_Value_x : %s',num2str(S.add.domain(:,S.add.objective==S.add.minimum_Value(end,2))'));
    fprintf('\n Minimum_Value : %g\n\n',S.add.minimum_Value(end,2));

    filename = sprintf('Dace_%d.mat', iter);
    save(filename, '-struct', 'S');
end

% MP+EI  일단 나중에~
% bayesian optimization with mixed integer