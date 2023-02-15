function S = Correlation(S)

% Gaussian Process Regression
S.Hypopt.r = zeros(size(S.add.domain,2),1,S.prob.dim);
S.acqui.x = S.acqui.x';

for i = 1:S.prob.dim
    Mxx = S.acqui.x(i)^2;
    Kxx = (S.add.domain(i,:).^ 2);
    Nxy = S.acqui.x(i)'*S.add.domain(i,:);

    S.Hypopt.r(:,:,i) = Mxx + Kxx' - 2*Nxy';
end

end