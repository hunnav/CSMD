function [S] = Input_struc

%% Variable

% Problem
S.prob.f_name = 'G9';                                                                                     % Function name to save data
S.prob.f = @(a,b,c,d,e,f,g) (a-10).^2+5*(b-12).^2+c.^4+3*(d-11).^2+10*e.^6+7*f.^2+g.^4-4*f.*g-10*f-8*g;   % Objective Function
S.prob.c1 = @(a,b,c,d,e,f,g) 127-2*a.^2-3*b.^4-c-4*d.^2-5*e;                                              % Constraint1 (plz defrom the equation to be "constraint >=0" form)
S.prob.c2 = @(a,b,c,d,e,f,g) 282-7*a-3*b-10*c.^2-d+e;                                                     % Constraint2 (plz defrom the equation to be "constraint >=0" form)
S.prob.c3 = @(a,b,c,d,e,f,g) 196-23*a-b.^2-6*f.^2+8*g;                                                    % Constraint3 (plz defrom the equation to be "constraint >=0" form)
S.prob.c4 = @(a,b,c,d,e,f,g) -4*a.^2-b.^2+3*a.*b-2*c.^2-5*f+11*g;                                         % Constraint4 (plz defrom the equation to be "constraint >=0" form)
S.prob.dim = 7;                                                                                           % Dimension of the problem
S.prob.numconstraint = 4;                                                                                 % Number of constraints
S.prob.maxiter = 300;                                                                                     % Maximum iterations for Bayesopt
S.prob.min = -10;                                                                                         % Minimum value of design variables
S.prob.max = 10;                                                                                          % Maximum value of design variables
S.prob.numinitsam = 10;                                                                                   % # of initial samples    
S.prob.surconst = false;                                                                                  % Using surrogate constraint mode, different with MP (true,false)

% Hyperparameter optimization
S.Hypopt.solver = 'dace';                                                                                 % Maximum likelihood estimation solver ('fmincon','fminunc','ga','ga+fmincon','pso','dace')
S.Hypopt.min = 0;                                                                                         % Minimum value of hyperparameter
S.Hypopt.max = 10;                                                                                        % Maximum value of hyperparameter
S.Hypopt.validation = true;                                                                               % Using training-validation mode (true,false)
S.Hypopt.validation_ratio = 0.9;                                                                          % Only for validation mode, Validation ratio (0~1)
S.Hypopt.validation_trynum = 10;                                                                          % Only for validation mode, number of initial theta trials.
S.Hypopt.frequency = 10;                                                                                  % How often to calculate hyperparameters, only when >= S.prob.dim*100
S.Hypopt.initheta = zeros(1,S.prob.dim);                                                                  % Initial theta
S.Hypopt.dace_reg = 'regpoly1';                                                                           % Only for Dace, Regression function in dace method ('regpoly0','regpoly1','regpoly2','regpoly3','regpolyauto') (Mininum # : 1, (1+dim), (dim+1)*(dim+2)/2, (dim+1)*(dim+2)*(dim+3)/6)
S.Hypopt.dace_cor = 'correxp';                                                                            % Only for Dace, Correlation function in dace method ('corrcubic','correxp','correxpg','corrgauss','corrlin','corrspherical','corrspline')

% Acquisition function
S.acqui.mode = 'EI';                                                                                      % Acquisition function for bayesian optimization ('PI', 'EI', 'LCB', 'UCB', 'MP', 'MP+EI')
S.acqui.solver = 'ga+fmincon';                                                                            % Solver for the acquisition function ('fmincon','ga','ga+fmincon','pso','multi_start')
S.acqui.constmode = 'same';                                                                               % Only for surrogate constraint mode, Determine initial constraint domain ('same','different','seperate')
S.acqui.mindis = 0.01;                                                                                    % Minimum distance between samples

% Add_point
S.add.scaling = 'yeojohnson';                                                                             % Choose the scailing mode ('boxcox', 'yeojohnson', 'log', 'log2')                                                             
S.add.log_domscale = 20;                                                                                  % Only for log and log2, Scale for the domain


%% Default

% Problem
S.prob.filename = sprintf('%s_%d_%s_%d_%s_%s_%s', ...                                                     % File name (Hyperparameter solver_validation_acqui_acqui solver)
    S.prob.f_name,S.prob.maxiter,S.Hypopt.solver,S.Hypopt.validation, ...
    S.acqui.mode,S.acqui.solver,S.add.scaling);

% Hyperparameter optimization
if S.prob.surconst == true
    S.Hypopt.theta = zeros(1+S.prob.numconstraint,S.prob.dim);                                            % Optimize hyperparameters from MLE process
else
    S.Hypopt.theta = zeros(1,S.prob.dim);                                                                 % Optimize hyperparameters from MLE process
end
S.Hypopt.R = zeros(S.prob.numinitsam,S.prob.numinitsam,S.prob.dim);                                       % Correlation matrix between observations
S.Hypopt.r = zeros(S.prob.numinitsam,1,S.prob.dim);                                                       % Correlation matrix among observations and new data point
S.Hypopt.invR = zeros(S.prob.numinitsam,S.prob.numinitsam);                                               % Inverse matrix of R matrix
if S.prob.surconst == true
    S.Hypopt.alpha = zeros(1+S.prob.numconstraint,1);                                                     % Alpha from kriging
    S.Hypopt.sigma = zeros(1+S.prob.numconstraint,1);                                                     % Sigma from kriging
else
    S.Hypopt.alpha = 0;                                                                                   % Alpha from kriging
    S.Hypopt.sigma = 0;                                                                                   % Sigma from kriging
end

% Acquisition function
S.acqui.exploratio = 0;                                                                                   % Exploration ratio between exploitation and exploration
S.acqui.minobj = 0;                                                                                       % Current feasible minimum value
S.acqui.x = zeros(S.prob.dim,1);                                                                          % New sample

% Add_point
% S.add.domain = (S.prob.max-S.prob.min).*lhsdesign(S.prob.dim,S.prob.numinitsam) + (S.prob.min);           % Initial samples for Bayesopt (Latin hyper Cube sampling)
S.add.domain = [-8.58675620812388 -6.39015938271466 -4.13290068298605 -8.73922222746539	6.61028657019135 -1.94974747676626 8.73838320483492 -8.67671267639097 1.62402507112903 5.28067216773188;
    6.05957938481184	-9.86934724613246	6.74704959884142	8.12841436909128	-6.83479777590614	7.71458410823048	-7.80309737241254	-1.23446710354779	-7.79462233412538	-1.25125215933072;
    7.16985146038418	-0.761081186790541	-7.23817699685566	2.55630915902762	8.91438625222899	0.489652209552931	4.10391368291101	-1.58348934760120	-0.156536098023999	-5.06114936076976;
    -6.41206713345401	3.01054723572729	7.42173375063370	-1.65867095466551	-2.00137555513768	3.38111795234074	4.32289772533690	3.23609574600468	-5.62313244560910	7.67316499603713;
    0.100318077678697	5.30867496807755	3.67512002850975	5.74295197185217	-8.08957241047017	6.33992611175073	-3.60327245668628	7.31727961324541	4.43385629602521	-3.68228004144281;
    2.30730136372282	9.75272579305742	-5.78968200309360	-1.10863638836131	2.91684053468354	-8.41445758851710	-7.01511778678660	-6.60769827329952	-3.19063207616620	3.49626158617308;
    -1.72518756874041	-3.11788473310928	0.0719233000777759	-6.77227488931732	1.09679744149467	-6.23498846417084	1.34553116895889	5.76750509838933	8.36756025697434	-8.79144278654358];
S.add.objective = transpose(S.prob.f(S.add.domain(1,:), S.add.domain(2,:),...                             % Objective function value of initial samples for Bayesopt
                                     S.add.domain(3,:), S.add.domain(4,:),...
                                     S.add.domain(5,:), S.add.domain(6,:),...
                                     S.add.domain(7,:)));
if S.prob.surconst == true
    switch S.acqui.constmode
        case 'same'
            for i = 1:S.prob.numconstraint
                S.add.original_constdomainy(:,i) = transpose(S.prob.(['c', num2str(i)])( ...              % Domainy of constraint surrogate model
                    S.add.domain(1,:), S.add.domain(2,:), S.add.domain(3,:), S.add.domain(4,:), ...
                    S.add.domain(5,:), S.add.domain(6,:), S.add.domain(7,:)));
            end
        case 'different'
            S.add.iniconstdomain = (S.prob.max-S.prob.min).*lhsdesign( ...                                % Initial samples for constraint (Latin hyper Cube sampling)
                S.prob.dim,S.prob.numinitsam)+(S.prob.min);
            for i = 1:S.prob.numconstraint
                S.add.original_constdomainy(:,i) = transpose(S.prob.(['c', num2str(i)])( ...              % Domainy of constraint surrogate model
                    S.add.iniconstdomain(1,:), S.add.iniconstdomain(2,:), ...
                    S.add.iniconstdomain(3,:), S.add.iniconstdomain(4,:), ...
                    S.add.iniconstdomain(5,:), S.add.iniconstdomain(6,:), ...
                    S.add.iniconstdomain(7,:)));
            end
        case 'seperate'
            for i = 1:S.prob.numconstraint
                S.add.iniconstdomain(:,:,i) = (S.prob.max-S.prob.min ...                                  % Initial samples for constraint (Latin hyper Cube sampling)
                    ).*lhsdesign(S.prob.dim,S.prob.numinitsam)+(S.prob.min);                              % Domainy of constraint surrogate model
                S.add.original_constdomainy(:,i) = transpose(S.prob.(['c', num2str(i)])( ...
                    S.add.iniconstdomain(1,:,i), S.add.iniconstdomain(2,:,i),...
                    S.add.iniconstdomain(3,:,i), S.add.iniconstdomain(4,:,i),...
                    S.add.iniconstdomain(5,:,i), S.add.iniconstdomain(6,:,i),...
                    S.add.iniconstdomain(7,:,i)));
            end
    end
    S.add.constdomainy = S.add.original_constdomainy;                                                     % For modified constdomainy
end
if  S.prob.surconst == true                                                                               % Original value of constraint function
    S.add.constraint = max([S.add.original_constdomainy(:,:), ...
        zeros(size(S.add.original_constdomainy, 1), 1)], [], 2);
else
    S.add.constraint = zeros(S.prob.numinitsam,1);
end
S.add.obPLUScon = S.add.objective;                                                                        % Objective plus constraint (Only use in PI and EI)
S.add.standard = zeros(S.prob.maxiter,5);                                                                 % Storage of the data for scaling
S.add.A = -min(S.add.objective)+1;                                                                        % Initial sum term of log function for objective function
S.add.B = 0;                                                                                              % Initial base of log function for objective function
S.add.C = exp(log(S.add.log_domscale+1)/(S.add.log_domscale/2));                                          % Initial base of log function for objective function 2
S.add.domainy = S.add.objective+S.add.A;                                                                  % y domain value for Bayesopt
S.add.minimum_Value = [0,inf];                                                                            % Record of minimum value
S.add.cnt = 0;                                                                                            % number of iteration

% Print
S.print.x = nan(S.prob.maxiter,1);                                                                        % x value
S.print.y_min = nan(S.prob.maxiter,1);                                                                    % y value (min norm)
S.print.y_current = nan(S.prob.maxiter,1);                                                                % y value (current distance)
S.print.near = inf(1,S.prob.dim+2);                                                                       % Save the nearest point

end