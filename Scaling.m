function S = Scaling(S,n)

S.add.standard(S.add.cnt,1) = S.add.cnt;

if S.add.obPLUScon(n,1)+S.add.A < 1
    S.add.A = -S.add.obPLUScon(n,1)+1;
end

S.add.standard(S.add.cnt,2) = median(S.add.obPLUScon+S.add.A);
S.add.B = exp(log(S.add.standard(S.add.cnt,2))/S.add.domscale);
if S.add.reversescale == true
    Reverse_Modified_Objective = S.add.domscale - log(S.add.obPLUScon+S.add.A)/log(S.add.B);
    Reverse_Modified_Objective(Reverse_Modified_Objective < 0) = 0;
    S.add.domainy = S.add.domscale-log(Reverse_Modified_Objective+1)/log(S.add.C);
else
    S.add.domainy = log(S.add.obPLUScon+S.add.A-1+S.add.B^30)/log(S.add.B)-30;
end

S.add.standard(S.add.cnt,3:5) = [S.add.A,S.add.B,S.add.C];

end