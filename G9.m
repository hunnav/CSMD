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

low_Range = -10;                   % x range min
upper_Range = 10;                  % x range max 
max_iter = 1500;                   % Maximum iterations
initial_theta = [0,0,0,0,0,0,0];   % Initial theta (guess)
EI_mode = 'optimoptions';          % EI_mode


Domain = (upper_Range - low_Range)*lhsdesign(7,20)+low_Range;  % initial domain
Domain_y = zeros(size(Domain,2),1);
x = zeros(7,1);
for i = 1 : size(Domain,2)    % for all initial domain
    x(1:7) = Domain(1:7,i);
    Domain_y(i,1) = ff(x(1),x(2),x(3),x(4),x(5),x(6),x(7));   % initial 'y' value when domain is given
end
miniter = 0;                  % number of iteration

%% 1.2 Iteration for BayesOpt

% Infilling criterion : EI process
while miniter < max_iter       % until to be maximum iterations
    tic
    min_obj = min(Domain_y);   % current minimum value

    % Hyperparameter optimization with new samples based on MLE (maximum likelyhood estimation)
    [theta,alpha_kriging,sigma,inv_R] = optimizeHypes(initial_theta, Domain, Domain_y , 'fmincon');
    
    % EI process to extract the new point(Dom_EI)
    Dom_EI = EIval_check(EI_mode,Domain,Domain_y,theta,sigma,alpha_kriging,inv_R,low_Range,upper_Range,min_obj);

    x = Dom_EI';
    add = 0;
    gg1 = g1(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
    if gg1 < 0
        add = add - gg1;
    end
    gg2 = g2(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
    if gg2 < 0
        add = add - gg2;
    end
    gg3 = g3(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
    if gg3 < 0
        add = add - gg3;
    end
    gg4 = g4(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
    if gg4 < 0
        add = add - gg4;
    end

    % Add sample from the EI process
    Domain = [Domain,x];                % add the new point x to domain
    Domain_y = [Domain_y;ff(x(1),x(2),x(3),x(4),x(5),x(6),x(7))+add];         % add the new caluclated y to domain
    miniter = miniter + 1    % add the number of iteration
    Current_Minimum_Value = min_obj
    toc
end

Minimum_Value = min(Domain_y);
row = find(Domain_y==Minimum_Value);
Minimum_Value_x = Domain(:,row)'
Minimum_Value
%test