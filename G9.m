clear;
clc;
rng("shuffle")
delete(gcp('nocreate'))        % returns the current pool if one exists, otherwise pool will be empty & delete it
parpool('threads')             % creates and returns a thread-based pool.

%% 1.Setting for BayesOpt

ff= @(a,b,c,d,e,f,g) (a-10).^2+5*(b-12).^2+c.^4+3*(d-11).^2+10*e.^6+7*f.^2+g.^4-4*f.*g-10*f-8*g; % function

g1 = @(a,b,c,d,e,f,g) 127-2*a^2-3*b^4-c-4*d^2-5*e;
g2 = @(a,b,c,d,e,f,g) 282-7*a-3*b-10*c^2-d+e;
g3 = @(a,b,c,d,e,f,g) 196-23*a-b^2-6*f^2+8*g;
g4 = @(a,b,c,d,e,f,g) -4*a^2-b^2+3*a*b-2*c^2-5*f+11*g;

max_iter = 1500;                   % Maximum iterations
low_Range = -10;                   % x range min
upper_Range = 10;                  % x range max 
num_initial_value = 20;            % # of initial value
Initial_theta = [0,0,0,0,0,0,0];   % Initial theta (guess)
MLE_mode = 'fmincon';              % MLE_mode
EI_mode = 'ga';                    % EI_mode ('ga','pso','fmincon','ga+fmincon','multi_start')
EI_acq_mode = 'normal';            % EI_acq_mode('normal','weighted')
divider = 2;                       % How often to calculate hyperparameters
ratio_or_weight = 0;               % Default ratio is 0 and Default weight is 0.5
beta = 0.05;                       % beta must be between 0 to 1(all)
    
miniter = 0;     % number of iteration
dim = size(Initial_theta,2);
R = zeros(num_initial_value,num_initial_value,dim);
r = zeros(num_initial_value,1,dim);
theta = zeros(1,dim);
Domain = (upper_Range - low_Range)*lhsdesign(dim,num_initial_value)+low_Range;  % initial domain
Domain_y = transpose(ff(Domain(1,:), Domain(2,:), Domain(3,:), Domain(4,:), Domain(5,:), Domain(6,:), Domain(7,:)));

%% 2.Iteration for BayesOpt

% Infilling criterion : EI process
while miniter < max_iter       % until to be maximum iterations
    tic
    min_obj = min(Domain_y);   % current minimum value

    % Hyperparameter optimization with new samples based on MLE (maximum likelyhood estimation)
    [theta,alpha_kriging,sigma,inv_R,R] = optimizeHypes(Initial_theta, theta, Domain, Domain_y, R, r, miniter, divider, MLE_mode);
    
    % EI process to extract the new point(Dom_EI)
    ratio_or_weight_modify = ratio_or_weight;
    while 1
        Dom_EI = EIval_check(Domain, Domain_y, theta, sigma, alpha_kriging, inv_R, low_Range, upper_Range, min_obj, EI_mode, EI_acq_mode, ratio_or_weight_modify);
        r = Correlation(Domain,Dom_EI);
        if min(sum(r, 3)) > beta
            break
        else
           if strcmp(EI_acq_mode,'normal')
               ratio_or_weight_modify = ratio_or_weight_modify+0.5;
           else
               ratio_or_weight_modify = ratio_or_weight_modify-0.5;
           end
        end
    end
    x = Dom_EI';
    gg1 = 0-g1(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    gg2 = 0-g2(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    gg3 = 0-g3(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    gg4 = 0-g4(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    objective = ff(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
    p0 = 10^(floor(1+log10(abs(objective))));
    
    % Add sample from the EI process
    n = size(Domain, 2);
    Domain(:,n+1) = x;   % add the new point x to domain
    Domain_y(n+1,1) = objective+p0*max([gg1,gg2,gg3,gg4,0]);    % we need to change it according to the conditon, add the new caluclated y to domain

    miniter = miniter + 1    % add the number of iteration
    Current_Minimum_Value = min_obj
    toc
end

%% 3.Result
Minimum_Value = min(Domain_y);
row = find(Domain_y==Minimum_Value);
Minimum_Value_x = Domain(:,row)'
Minimum_Value
