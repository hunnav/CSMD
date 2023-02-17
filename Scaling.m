function S = Scaling(S,n)

S.add.standard(S.add.cnt,1) = S.add.cnt;

if strcmp(S.acqui.mode,'MP')
    value = S.add.objective;
else
    value = S.add.obPLUScon;
end

if value(n,1)+S.add.A < 1
    S.add.A = -value(n,1)+1;
end

S.add.standard(S.add.cnt,2) = median(value+S.add.A);
S.add.B = exp(log(S.add.standard(S.add.cnt,2))/S.add.domscale);
Reverse_Modified_Objective = S.add.domscale - log(value+S.add.A)/log(S.add.B);
Reverse_Modified_Objective(Reverse_Modified_Objective < 0) = 0;
S.add.domainy = S.add.domscale-log(Reverse_Modified_Objective+1)/log(S.add.C);
S.add.standard(S.add.cnt,3:5) = [S.add.A,S.add.B,S.add.C];

end