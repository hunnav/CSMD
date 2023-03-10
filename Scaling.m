function S = Scaling(S,n)

S.add.standard(S.add.cnt,1) = S.add.cnt;

switch S.add.scaling
    case 'boxcox'
        [S.add.domainy,lambda] = boxcox(S.add.obPLUScon);
        S.add.standard(S.add.cnt,2:5) = [lambda,0,0,0];
    case 'yeojohnson'
        [S.add.domainy,lambda] = palm_yeojohnson(S.add.obPLUScon);
        S.add.standard(S.add.cnt,2:5) = [lambda,0,0,0];
    otherwise
        if S.add.obPLUScon(n,1) + S.add.A < 1
            S.add.A = -S.add.obPLUScon(n,1)+1;
        end
        S.add.standard(S.add.cnt,2) = median(S.add.obPLUScon+S.add.A);
        S.add.B = exp(log(S.add.standard(S.add.cnt,2))/S.add.log_domscale);
        if strcmp(S.add.scaling,'log2')
            Reverse_Modified_Objective = S.add.log_domscale - log(S.add.obPLUScon+S.add.A)/log(S.add.B);
            Reverse_Modified_Objective(Reverse_Modified_Objective < 0) = 0;
            S.add.domainy = S.add.log_domscale-log(Reverse_Modified_Objective+1)/log(S.add.C);
        else
            S.add.domainy = log(S.add.obPLUScon+S.add.A-1+S.add.B^30)/log(S.add.B)-30;
        end
        S.add.standard(S.add.cnt,3:5) = [S.add.A,S.add.B,S.add.C];
end

if S.prob.surconst == 1
    S.add.constdomainy = S.add.original_constdomainy;
    for i = 1 : S.prob.numconstraint
        switch S.add.scaling
            case 'boxcox'
                S.add.constdomainy(:,i) = boxcox(S.add.original_constdomainy(:,i));
            case 'yeojohnson'
                S.add.constdomainy(:,i) = palm_yeojohnson(S.add.original_constdomainy(:,i));
            otherwise
                A = -min(S.add.original_constdomainy(:,i))+1;
                B = exp(log(median(S.add.original_constdomainy(:,i)+A))/S.add.log_domscale);
                if strcmp(S.add.scaling,'log2')
                    Reverse_Modified_Objective = S.add.log_domscale - log(S.add.original_constdomainy(:,i)+A)/log(B);
                    Reverse_Modified_Objective(Reverse_Modified_Objective < 0) = 0;
                    S.add.constdomainy(:,i) = S.add.log_domscale-log(Reverse_Modified_Objective+1)/log(S.add.C);
                else
                    S.add.constdomainy(:,i) = log(S.add.original_constdomainy(:,i)+A-1+B^30)/log(B)-30;
                end
        end
    end
end

end