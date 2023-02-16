function S = New_point(S)

if  S.add.minimum_Value(end,2) == inf
    S.acqui.minobj = min(S.add.domainy);
else
    S.acqui.minobj = S.add.domainy(S.prob.numinitsam + S.add.minimum_Value(end,1));
end

while 1
    while 1
        try
            S = Acq_solver(S);
            break
        catch
            disp('There is some error, repeat again.')
        end
    end
    
    S.Hypopt.r = zeros(size(S.add.domain,2),1,S.prob.dim);
    for i = 1:S.prob.dim
        S.Hypopt.r(:,:,i) = S.acqui.x(i)^2 + (S.add.domain(i,:).^ 2)' - 2*S.acqui.x(i)'*S.add.domain(i,:)';
    end

    if sqrt(min(sum(S.Hypopt.r, 3))) > S.acqui.mindis
        break
    else
        if abs(S.acqui.exploratio) >= 5
            break
        else
            S.acqui.exploratio  = S.acqui.exploratio +1;
        end
    end
end

end