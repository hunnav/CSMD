function [S] = Input_struc

%% Problem

% Variable
S.prob.f = @(a,b,c,d,e,f,g) (a-10).^2+5*(b-12).^2+c.^4+3*(d-11).^2+10*e.^6+7*f.^2+g.^4-4*f.*g-10*f-8*g;   % Objective Function
S.prob.c1 = @(a,b,c,d,e,f,g) 127-2*a^2-3*b^4-c-4*d^2-5*e;                                                 % Constraint1
S.prob.c2 = @(a,b,c,d,e,f,g) 282-7*a-3*b-10*c^2-d+e;                                                      % Constraint2
S.prob.c3 = @(a,b,c,d,e,f,g) 196-23*a-b^2-6*f^2+8*g;                                                      % Constraint3
S.prob.c4 = @(a,b,c,d,e,f,g) -4*a^2-b^2+3*a*b-2*c^2-5*f+11*g;                                             % Constraint4
S.prob.dim = 7;                                                                                           % Dimension of the problem
S.prob.maxiter = 2000;                                                                                    % Maximum iterations for Bayesopt
S.prob.max = 10;                                                                                          % Design variables' range max
S.prob.min = -10;                                                                                         % Design variables' range min
S.prob.numinitsam = 10;                                                                                   % # of initial samples

%% Hyperparameter optimization

% Variable
S.Hypopt.solver = 'fmincon';                                                                              % Maximum likelihood estimation solver ('fmincon','fminunc','ga','pso')
S.Hypopt.frequency = 5;                                                                                   % How often to calculate hyperparameters
S.Hypopt.initheta = zeros(1,S.prob.dim);                                                                  % Initial theta (guess)

% Use Default
S.Hypopt.theta = zeros(1,S.prob.dim);                                                                     % Optimize hyperparameters from MLE process
S.Hypopt.R = zeros(S.prob.numinitsam,S.prob.numinitsam,S.prob.dim);                                       % Correlation matrix between observations
S.Hypopt.r = zeros(S.prob.numinitsam,1,S.prob.dim);                                                       % Correlation matrix among observations and new data point
S.Hypopt.invR = zeros(S.prob.numinitsam,S.prob.numinitsam);                                               % Inverse matrix of R matrix
S.Hypopt.alpha = 0;                                                                                       % Alpha from kriging
S.Hypopt.sigma = 0;                                                                                       % Sigma from kriging

%% Acquisition function

% Variable
S.acqui.mode = 'EI';                                                                                      % Acquisition function for bayesian optimization ('PI', 'EI', 'LCB', 'UCB', 'MP', 'EIMP')
S.acqui.solver = 'ga+fmincon';                                                                            % Solver for the acquisition function ('fmincon','ga','ga+fmincon','pso','multi_start')
S.acqui.mindis = 0.01;                                                                                    % Minimum distance between samples

% Use Default
S.acqui.exploratio = 0;                                                                                   % Exploration ratio between exploitation and exploration
S.acqui.minobj = 0;                                                                                       % Current feasible minimum value
S.acqui.x = zeros(S.prob.dim,1);                                                                          % New sample

%% Add_point

% Variable
S.add.domscale = 100;                                                                                     % scale for the domain

% Use Default
S.add.domain = (S.prob.max-S.prob.min).*lhsdesign(S.prob.dim,S.prob.numinitsam) + (S.prob.min);           % Initial samples for Bayesopt(Latin hyper Cube sampling)
S.add.objective = transpose(S.prob.f(S.add.domain(1,:), S.add.domain(2,:),...                             % Objective function value of initial samples for Bayesopt
                                     S.add.domain(3,:), S.add.domain(4,:),...
                                     S.add.domain(5,:), S.add.domain(6,:),...
                                     S.add.domain(7,:)));
S.add.Modified_Objective = zeros(S.prob.numinitsam,1);                                                    % Modified objective function value
S.add.constraint = zeros(S.prob.numinitsam,1);                                                            % Original value of constraint function
S.add.add = zeros(S.prob.numinitsam,1);                                                                   % Modified value of constraint function
S.add.standard = zeros(S.prob.maxiter,7);                                                                 % Storage of the data for scaling
S.add.A = -min(S.add.objective)+1;                                                                        % Initial sum term of log function for objective function
S.add.B = 0;                                                                                              % Initial base of log function for objective function
S.add.C = 0;                                                                                              % Initial base of log function for constraint function
S.add.D = exp(log(S.add.domscale+1)/(S.add.domscale/2));                                                  % Initial base of log function for objective function 2
S.add.domainy = S.add.objective+S.add.A;                                                                  % y domain value for Bayesopt
S.add.minimum_Value = [0,inf];                                                                            % Record of minimum value
S.add.cnt = 0;                                                                                            % number of iteration

end