function Minimum_Value = Add_new_point(num_initial_value, Iteration, Domain, Objective, Constraint, Standard, K, f, g1, g2, g3, g4, )

    n = num_initial_value + Iteration + 1;
    Domain(:,n) = x; 
    Objective(n,1) = f(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon

    gg1 = 0-g1(x(1),x(2),x(3),x(4),x(5),x(6),x(7));            % we need to change it according to the conditon
    gg2 = 0-g2(x(1),x(2),x(3),x(4),x(5),x(6),x(7));            % we need to change it according to the conditon
    gg3 = 0-g3(x(1),x(2),x(3),x(4),x(5),x(6),x(7));            % we need to change it according to the conditon
    gg4 = 0-g4(x(1),x(2),x(3),x(4),x(5),x(6),x(7));            % we need to change it according to the conditon
    Constraint(n,1) = max([gg1,gg2,gg3,gg4,0]);
    
    [Standard,Domain_y,Modified_Objective,Add] = Scaling(Standard,Iteration,Objective,Constraint,K);

    if (Add(end) == 0 && Objective(end) < Minimum_Value(end,2))
        Minimum_Value(end+1,:) = [Iteration+1, Objective(end)];
    end

end
