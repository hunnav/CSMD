%% Q4 Isoparametric   

clear
close all
clc
tic

%% Variables

L = 20; % Length (m)
H = 5; % Height (m)
T = 0.025; % Thickness (m)
E = 210*10^9; % Young's modulus (Pa)
v = 0.3; % Poisson's ratio

% Force (N)
xyFxFy_condition = {20,5,0,-10000000}; % To add more force, you can use another row , 'a' means all coordinate

% Boundary Conditions (Fixed coordinate) (m)
xy_fix = {0,'a'};  % To add BC, you can use another row , 'a' means all coordinate

% # of mesh
i = 100; % x mesh #
j = 25; % y mesh #

%% Original Coordinate of Nodes

x_coordinate = zeros((i+1)*(j+1),1);
y_coordinate = zeros((i+1)*(j+1),1);
for Node_num = 1:(i+1)*(j+1)
    a = fix((Node_num-1)/(j+1));
    b = rem(Node_num-1,(j+1));
    x_coordinate(Node_num) = a*L/i;
    y_coordinate(Node_num) = b*H/j;
end

%% Boundary condition for fixed coordinate

Node_fixed = zeros((i+1)*(j+1),1);

for BCN = 1:numel(xy_fix)/2
    for Node_num = 1:(i+1)*(j+1)
        if or(isequal(x_coordinate(Node_num), xy_fix{BCN,1}), strcmp(xy_fix(BCN,1),'a'))
            if or(isequal(y_coordinate(Node_num), xy_fix{BCN,2}), strcmp(xy_fix(BCN,2),'a'))
                Node_fixed(Node_num) = 1;
            end
        end
    end
end

%% F

F = zeros(2*(i+1)*(j+1),1);

for FCN = 1:numel(xyFxFy_condition)/4
    for Node_num = 1:(i+1)*(j+1)
        if or(isequal(x_coordinate(Node_num), xyFxFy_condition{FCN,1}), strcmp(xyFxFy_condition(FCN,1),'a'))
            if or(isequal(y_coordinate(Node_num), xyFxFy_condition{FCN,2}), strcmp(xyFxFy_condition(FCN,2),'a'))
                F(Node_num*2-1) = xyFxFy_condition{FCN,3};
                F(Node_num*2) = xyFxFy_condition{FCN,4};
            end
        end
    end
end

%% K

% Local k & Global K

s = [-0.57735 ; -0.57735 ; 0.57735 ; 0.57735]; % 4 Points
t = [-0.57735 ; 0.57735 ; -0.57735 ; 0.57735];
W = ones(4,1);

D = E/(1-v^2)*[1 v 0; v 1 0; 0 0 (1-v)/2];

K = sparse(zeros((i+1)*(j+1)*2,(i+1)*(j+1)*2));
for Element_number = 1:i*j
    a = fix((Element_number-1)/j);
    b = rem(Element_number-1,j);
    Ele_Node_Num(1) = (j+1)*a+b+1;   %1
    Ele_Node_Num(2) = (j+1)*a+b+2+j; %7
    Ele_Node_Num(3) = (j+1)*a+b+3+j; %8
    Ele_Node_Num(4) = (j+1)*a+b+2;   %2
    
    x_1_4 = zeros(4,1);
    y_1_4 = zeros(4,1);
    x_1_4(1:4) = x_coordinate(Ele_Node_Num(1:4));  % x1, x2 .... y4
    y_1_4(1:4) = y_coordinate(Ele_Node_Num(1:4));
    
    k = zeros(8,8);
    J = zeros(4:1);
    B = zeros(3,8);
    for n = 1:4 
        for_jac = [0 1-t(n) t(n)-s(n) s(n)-1;
                   t(n)-1 0 s(n)+1 -s(n)-t(n);
                   s(n)-t(n) -s(n)-1 0 t(n)+1;
                   1-s(n) s(n)+t(n) -t(n)-1 0];
        
        J(n) = 1/8 * transpose(x_1_4)*for_jac*y_1_4;     % jacobian
        
        dNds = [  t(n)/4-1/4   % Shape Function
                  1/4-t(n)/4
                  t(n)/4+1/4
                - t(n)/4-1/4];

        dNdt = [  s(n)/4-1/4
                - s(n)/4-1/4
                  s(n)/4+1/4
                  1/4-s(n)/4];
              
        b_add = zeros(3,2,4);

        a_B = 1/4*(y_1_4(1)*(s(n)-1)+y_1_4(2)*(-1-s(n))+y_1_4(3)*(1+s(n))+y_1_4(4)*(1-s(n)));
        b_B = 1/4*(y_1_4(1)*(t(n)-1)+y_1_4(2)*(1-t(n))+y_1_4(3)*(1+t(n))+y_1_4(4)*(-1-t(n)));
        c_B = 1/4*(x_1_4(1)*(t(n)-1)+x_1_4(2)*(1-t(n))+x_1_4(3)*(1+t(n))+x_1_4(4)*(-1-t(n)));
        d_B = 1/4*(x_1_4(1)*(s(n)-1)+x_1_4(2)*(-1-s(n))+x_1_4(3)*(1+s(n))+x_1_4(4)*(1-s(n)));

        for nn = 1:4
            b_add(:,:,nn) = [a_B*(dNds(nn))-b_B*(dNdt(nn)) 0 ; 0 c_B*(dNdt(nn))-d_B*(dNds(nn)); c_B*(dNdt(nn))-d_B*(dNds(nn)) a_B*(dNds(nn))-b_B*(dNdt(nn))];
        end

        B = 1/J(n)*[b_add(:,:,1), b_add(:,:,2), b_add(:,:,3), b_add(:,:,4)];

        k = k+transpose(B)*D*B*J(n)*T*W(n)*W(n);
    end
    
    for jj = 1:4
        bejj = Ele_Node_Num(jj);
        for ii = 1:4
            beii = Ele_Node_Num(ii);
            K(bejj*2-1:bejj*2, beii*2-1:beii*2) = K(bejj*2-1:bejj*2, beii*2-1:beii*2) + k(jj*2-1:jj*2,ii*2-1:ii*2);
        end
    end

end

%% Solve

F_Known = F;
K_Part = K;
Delete_num = 0;

for Node_num = 1:(i+1)*(j+1)
    if Node_fixed(Node_num) == 1
        F_Known(2*(Node_num-Delete_num)-1:2*(Node_num-Delete_num)) = [];
        K_Part(2*(Node_num-Delete_num)-1:2*(Node_num-Delete_num),:) = [];
        K_Part(:,2*(Node_num-Delete_num)-1:2*(Node_num-Delete_num)) = [];
        Delete_num = Delete_num+1;
    end
end

% Cholesky Factorization
% R = chol(K_Part);
% d_unknown = R\(R'\F_Known);

% Inverse
d_unknown = K_Part\F_Known;

%% Result

% Displacement
d = zeros(2*(i+1)*(j+1),1);
d_unknown_num = 1;

for Node_num = 1:(i+1)*(j+1)
    if Node_fixed(Node_num) == 1
        d(2*Node_num-1:2*Node_num) = 0;
    else
        d(2*Node_num-1:2*Node_num) = d_unknown(d_unknown_num:d_unknown_num+1);
        d_unknown_num = d_unknown_num+2;
    end
end

d_x(1:(j+1)*(i+1),1) = d(2*(1:(j+1)*(i+1))-1);  % displacement
d_y(1:(j+1)*(i+1),1) = d(2*(1:(j+1)*(i+1)));

d_xx = x_coordinate + d_x;   % deformed coordinate
d_yy = y_coordinate + d_y;

% Force
F = K*d;

%% Plot
    
% Strain
for_jac = [0 1 0 -1;
          -1 0 1 0;
           0 -1 0 1;
           1 0 -1 0];

dNds = [ -1/4    % For B0 matrix
          1/4
          1/4
         -1/4 ];

dNdt = [ -1/4
         -1/4
          1/4
          1/4 ];

B_0 = zeros(3,8);
Strain = zeros(3,1,i*j);
Stress = zeros(3,1,i*j);

for Element_number = 1:i*j   
    a = fix((Element_number-1)/j);   % Element Node Num
    b = rem(Element_number-1,j);
    Ele_Node_Num(1) = (j+1)*a+b+1;   %1
    Ele_Node_Num(2) = (j+1)*a+b+2+j; %7
    Ele_Node_Num(3) = (j+1)*a+b+3+j; %8
    Ele_Node_Num(4) = (j+1)*a+b+2;   %2

    x_1_4(1:4) = x_coordinate(Ele_Node_Num(1:4));  % x1, x2 .... y4
    y_1_4(1:4) = y_coordinate(Ele_Node_Num(1:4));

    J_0 = 1/8 * transpose(x_1_4)*for_jac*y_1_4;     % Jacobian

    b_add = zeros(3,2,4);
    d_ele = zeros(8,1);

    % B0 matrix
    a_B = 1/4*(y_1_4(1)*(-1)+y_1_4(2)*(-1)+y_1_4(3)*(1)+y_1_4(4)*(1));
    b_B = 1/4*(y_1_4(1)*(-1)+y_1_4(2)*(1)+y_1_4(3)*(1)+y_1_4(4)*(-1));
    c_B = 1/4*(x_1_4(1)*(-1)+x_1_4(2)*(1)+x_1_4(3)*(1)+x_1_4(4)*(-1));
    d_B = 1/4*(x_1_4(1)*(-1)+x_1_4(2)*(-1)+x_1_4(3)*(1)+x_1_4(4)*(1));

    for nn = 1:4 
        b_add(:,:,nn) = [a_B*(dNds(nn))-b_B*(dNdt(nn)) 0 ; 0 c_B*(dNdt(nn))-d_B*(dNds(nn)); c_B*(dNdt(nn))-d_B*(dNds(nn)) a_B*(dNds(nn))-b_B*(dNdt(nn))];
        d_ele(2*nn-1:2*nn) = d((Ele_Node_Num(nn)*2-1):(Ele_Node_Num(nn)*2));   % Element displacement
    end

    B_0 = 1/J_0*[b_add(:,:,1), b_add(:,:,2), b_add(:,:,3), b_add(:,:,4)];

    Strain(:,:,Element_number) = B_0*d_ele;
    
    % Stress
    Stress(:,:,Element_number) = D*Strain(:,:,Element_number);

    % Plot
    subplot(231)  % Model
    patch(d_xx(Ele_Node_Num(1:4)),d_yy(Ele_Node_Num(1:4)),'w');
    
    subplot(232)  % x displacement
    patch(x_1_4,y_1_4,mean(d_ele(1)+d_ele(3)+d_ele(5)+d_ele(7)));  %  to make line width 0, add ",'EdgeColor', 'none'"
    
    subplot(233)  % y displacement
    patch(x_1_4,y_1_4,mean(d_ele(2)+d_ele(4)+d_ele(6)+d_ele(8)));
    
    subplot(234)  % x Stress
    patch(x_1_4,y_1_4,Stress(1,1,Element_number));

    subplot(235)  % y Stress
    patch(x_1_4,y_1_4,Stress(2,1,Element_number));

    subplot(236)  % xy Shear Stress
    patch(x_1_4,y_1_4,Stress(3,1,Element_number));
end

subplot(231)  % Model
title("Q4"); xlabel('x'); ylabel('y')
axis equal;
xlim([0 1.1*L]); ylim([-H 2*H]);

subplot(232)  % x displacement
title("x Displacement"); xlabel('x'); ylabel('y');
axis equal; 
xlim([0 1.1*L]); ylim([-H 2*H]);
colorbar;

subplot(233)  % y displacement
title("y Displacement"); xlabel('x'); ylabel('y');
axis equal;
xlim([0 1.1*L]); ylim([-H 2*H]);
colorbar;

subplot(234) % x Stress
title("x Stress"); xlabel('x'); ylabel('y');
axis equal;
xlim([0 1.1*L]); ylim([-H 2*H]);
colorbar;

subplot(235)  % y Stress
title("y Stress"); xlabel('x'); ylabel('y');
axis equal; 
xlim([0 1.1*L]); ylim([-H 2*H]);
colorbar;

subplot(236)  % xy Shear Stress
title("xy Shear Stress"); xlabel('x'); ylabel('y');
axis equal;
xlim([0 1.1*L]); ylim([-H 2*H]);
colorbar;
toc

% test