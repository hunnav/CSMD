clear;
clc;
delete(gcp('nocreate'))        % returns the current pool if one exists, otherwise pool will be empty & delete it
parpool('threads')             % creates and returns a thread-based pool.

%% 1.1 Setting for BayesOpt

ff= @(a,b,c,d,e,f,g) (a-10)^2+5*(b-12)^2+c^4+3*(d-11)^2+10*e^6+7*f^2+g^4-4*f*g-10*f-8*g; % function

g1 = @(a,b,c,d,e,f,g) 127-2*a^2-3*b^4-c-4*d^2-5*e;
g2 = @(a,b,c,d,e,f,g) 282-7*a-3*b-10*c^2-d+e;
g3 = @(a,b,c,d,e,f,g) 196-23*a-b^2-6*f^2+8*g;
g4 = @(a,b,c,d,e,f,g) -4*a^2-b^2+3*a*b-2*c^2-5*f+11*g;

max_iter = 1500;                   % Maximum iterations
low_Range = -10;                   % x range min
upper_Range = 10;                  % x range max 
num_initial_value = 20;            % # of initial value
initial_theta = [0,0,0,0,0,0,0];   % Initial theta (guess)
MLE_mode = 'fmincon';              % MLE_mode
EI_mode = 'pso';                   % EI_mode
EI_acq_mode = 'weighted';          % EI_acq_mode
ratio_or_weight = 0.5;             % Default ratio is 0.2 and Default weight is 0.5
divider = 1;                       % How often to calculate hyperparameters
beta = 0.01;                       % beta must be between 0 to 1(all)

Domain = (upper_Range - low_Range)*lhsdesign(size(initial_theta,2),num_initial_value)+low_Range;  % initial domain
Domain_y = zeros(size(Domain,2),1);
x = zeros(size(Domain,1),1);
for i = 1 : size(Domain,2)    % for all initial domain
    x(1:size(Domain,1)) = Domain(1:size(Domain,1),i);
    Domain_y(i,1) = ff(x(1),x(2),x(3),x(4),x(5),x(6),x(7));   % we need to change it according to the conditon, initial 'y' value when domain is given
end

R = zeros(size(Domain,2),size(Domain,2));
r = 0;
miniter = 0;     % number of iteration

%% 1.2 Iteration for BayesOpt

% Infilling criterion : EI process
while miniter < max_iter       % until to be maximum iterations
    tic
    min_obj = min(Domain_y);   % current minimum value
    
    % Hyperparameter optimization with new samples based on MLE (maximum likelyhood estimation)
    if or(rem(miniter,divider)==0,miniter<100*size(Domain,1))
        [theta,alpha_kriging,sigma,R,inv_R] = optimizeHypes(initial_theta, Domain, Domain_y, R, r, miniter, MLE_mode);
    end
       
    % EI process to extract the new point(Dom_EI)
    ratio_or_weight_modify = ratio_or_weight;
    while 1
        Dom_EI = EIval_check(Domain, Domain_y, theta, sigma, alpha_kriging, inv_R, low_Range, upper_Range, min_obj, EI_mode, EI_acq_mode, ratio_or_weight_modify);
        r=Correlation(Domain,Dom_EI,theta);
        if min(abs(log(r))) > beta
            break
        else
           if strcmp(EI_acq_mode,'normal')
               ratio_or_weight_modify = ratio_or_weight_modify+0.1;
           else
               ratio_or_weight_modify = ratio_or_weight_modify-0.1;
           end
        end
    end

    x = Dom_EI';
    add = 0;
    gg1 = g1(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    if gg1 < 0
        add = add - gg1;
    end
    gg2 = g2(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    if gg2 < 0
        add = add - gg2;
    end
    gg3 = g3(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    if gg3 < 0
        add = add - gg3;
    end
    gg4 = g4(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon
    if gg4 < 0
        add = add - gg4;
    end

    % Add sample from the EI process
    Domain = [Domain,x];                % add the new point x to domain
    Domain_y = [Domain_y;ff(x(1),x(2),x(3),x(4),x(5),x(6),x(7))+add];         % % we need to change it according to the conditon, add the new caluclated y to domain
    miniter = miniter + 1    % add the number of iteration
    Current_Minimum_Value = min_obj
    toc
end

Minimum_Value = min(Domain_y);
row = find(Domain_y==Minimum_Value);
Minimum_Value_x = Domain(:,row)'
Minimum_Value
