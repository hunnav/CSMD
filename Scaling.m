function S = Scaling(S,n)

S.add.standard(S.add.cnt,1) = S.add.cnt;
if S.add.objective(n,1)+S.add.A < 1
    S.add.A = -S.add.objective(n,1)+1;
end
S.add.standard(S.add.cnt,2) = median(S.add.objective+S.add.A);
S.add.B = exp(log(S.add.standard(S.add.cnt,2))/S.add.domscale);
Reverse_Modified_Objective = S.add.domscale - log(S.add.objective+S.add.A)/log(S.add.B);
Reverse_Modified_Objective(Reverse_Modified_Objective < 0) = 0;
S.add.Modified_Objective  = S.add.domscale-log(Reverse_Modified_Objective+1)/log(S.add.D);
non_zero_numbers = S.add.constraint(S.add.constraint ~= 0);
if isempty(non_zero_numbers)
    S.add.standard(S.add.cnt,3) = 0;
    S.add.C = 5;
else
    S.add.standard(S.add.cnt,3) = median(non_zero_numbers);
    S.add.C = exp(log(S.add.standard(S.add.cnt,3))/(S.add.domscale /2));
    S.add.add = log(S.add.constraint+S.add.C)/log(S.add.C)-1;
end
S.add.domainy = S.add.Modified_Objective + S.add.add;
S.add.standard(S.add.cnt,4:7) = [S.add.A,S.add.B,S.add.A,S.add.D];

end