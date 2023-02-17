%% 0.Setting for Matlab

clear; close; clc;
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

    x(S.add.cnt) = S.add.cnt;
    y(S.add.cnt) = norm([2.3331998;1.942837;-0.479546;4.387811;-0.632617;1.039775;1.600996]-[S.add.domain(:,S.add.objective==S.add.minimum_Value(end,2))],'inf');   % Also can use pdist2
    addpoints(animatedline,x,y);
    drawnow;

    toc
end

%% 3.Result

fprintf('\n Minimum_Value_x : %s',num2str(S.add.domain(:,S.add.objective==S.add.minimum_Value(end,2))'));
fprintf('\n Minimum_Value : %g\n\n',S.add.minimum_Value(end,2));

% 새로 만든 MP 제대로 되는지 확인하기
% fmincon 왜 안되는지 확인