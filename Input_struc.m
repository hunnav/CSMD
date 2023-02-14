function [S] = Input_struc
% S.problem
% S.Hypopt
% S.acqui
% S.etc

%% Problem
S.prob.f = @(a,b,c,d,e,f,g) (a-10).^2+5*(b-12).^2+c.^4+3*(d-11).^2+10*e.^6+7*f.^2+g.^4-4*f.*g-10*f-8*g;   % Objective Function
S.prob.c1 = @(a,b,c,d,e,f,g) 127-2*a^2-3*b^4-c-4*d^2-5*e;                                                 % Constraint1
S.prob.c2 = @(a,b,c,d,e,f,g) 282-7*a-3*b-10*c^2-d+e;                                                      % Constraint2
S.prob.c3 = @(a,b,c,d,e,f,g) 196-23*a-b^2-6*f^2+8*g;                                                      % Constraint3
S.prob.c4 = @(a,b,c,d,e,f,g) -4*a^2-b^2+3*a*b-2*c^2-5*f+11*g;                                             % Constraint4
S.prob.dim = 7;                                                                                           % Dimension of the problem
S.prob.max = 10*ones(S.prob.dim,1);                                                                       % Design variables' range max
S.prob.min = -10*ones(S.prob.dim,1);                                                                      % Design variables' range min
S.prob.numinitsam = 10;                                                                                   % # of initial samples
S.prob.domain = (S.prob.max-S.prob.min).*lhsdesign(S.prob.dim,S.prob.numinitsam) + (S.prob.min);          % Initial samples for Bayesopt(Latin hyper Cube sampling)
S.prob.domainy = transpose(S.prob.f(S.prob.domain(1,:), S.prob.domain(2,:),...
                                    S.prob.domain(3,:), S.prob.domain(4,:),...
                                    S.prob.domain(5,:), S.prob.domain(6,:),...
                                    S.prob.domain(7,:)));

%% Hyperparameter optimization
S.Hypopt.initheta = zeros(1,S.prob.dim);                                                                  % Initial theta (guess)
S.Hypopt.solver = 'fmincon';                                                                              % Maximum likelihood estimation solver ( 'ga', 'pso', 'fmincon' )
S.Hypopt.frequency = 2;                                                                                   % How often to calculate hyperparameters (%% Orginally, "divider" %%)
S.Hypopt.theta = zeros(1,S.prob.dim);                                                                     % hyperparameters from MLE process
S.Hypopt.BLUP = 0;                                                                                        % Best linear unbiased predictor
S.Hypopt.R = zeros(S.prob.numinitsam,S.prob.numinitsam,S.prob.dim);                                       % Correlation matrix between observations
S.Hypopt.r = zeros(S.prob.numinitsam,1,S.prob.dim);                                                       % Correlation matrix among observations and new data point

%% Acquisition function
S.acqui.maxiter = 1500;                                                                                   % Maximum iterations for Bayesopt
S.acqui.solver = 'ga+fmincon';                                                                            % Solver for the acquisition function ('ga','pso','fmincon','ga+fmincon','multi_start')
S.acqui.mode = 'EI';                                                                                      % Acquisition function for bayesian optimization ( 'EI', 'PI', 'UCB', 'LCB' )
S.acqui.exploratio = 0;                                                                                   % Exploration ratio between exploitation and exploration

%% Etc.
S.etc.mindis = 0.01;                                                                                      % minimum distance between samples  %% beta %%
S.etc.cnt = 0;                                                                                            % number of iteration
S.etc.domscale = 100;                                                                                     % scale for the domain %% K %% )


Constraint = zeros(num_initial_value,1);
Add = zeros(num_initial_value,1);
Standard = zeros(max_iter,6);
A = -min(Objective)+1; 
D = exp(log(K+1)/(K/2)); 
Domain_y = Objective+A;
Minimum_Value = [0,inf];

end