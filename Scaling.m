function [Standard,Domain_y,Modified_Objective,Add] = Scaling(Standard,Iteration,Objective,Constraint,K)

    Standard(Iteration,1) = Iteration;
    if Objective(n,1)+A < 1
        A = -Objective(n,1)+1;    
    end
    Standard(Iteration,2) = median(Objective+A);
    B = exp(log(Standard(Iteration,2))/K);
    Reverse_Modified_Objective = K - log(Objective+A)/log(B);
    Reverse_Modified_Objective(Reverse_Modified_Objective < 0) = 0;
    Modified_Objective = K-log(Reverse_Modified_Objective+1)/log(D);
    non_zero_numbers = Constraint(Constraint ~= 0);
    if isempty(non_zero_numbers)
        Standard(Iteration,3) = 0;
        C = 5;
    else
        Standard(Iteration,3) = median(non_zero_numbers);
        C = exp(log(Standard(Iteration,3))/(K/2));
        Add = log(Constraint+C)/log(C)-1;
    end
    Domain_y = Modified_Objective + Add;
    Standard(Iteration,4:7) = [A,B,C,D];
end