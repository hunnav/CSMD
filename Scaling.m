function S = Scaling(S,n)

S.add.standard(S.add.cnt,1) = S.add.cnt;

switch S.add.scaling
    case 'minmax'
        minvalue = min(S.add.obPLUScon);
        maxmin = max(S.add.obPLUScon)-minvalue;
        S.add.domainy = (S.add.obPLUScon-minvalue)/maxmin;
    case 'boxcox'
        [S.add.domainy,lambda] = boxcox(S.add.obPLUScon);
        S.add.standard(S.add.cnt,2:5) = [lambda,0,0,0];
    case 'yeojohnson'
        [S.add.domainy,lambda] = palm_yeojohnson(S.add.obPLUScon,-3);
        S.add.standard(S.add.cnt,2:5) = [lambda,0,0,0];
    case 'none'
        S.add.domainy = S.add.obPLUScon;
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
        switch S.add.constscaling
            case 'minmax'
                value = [S.add.original_constdomainy(:,i);0];
                minvalue = min(value);
                maxmin = max(value)-minvalue;
                result= (value-minvalue)/maxmin;
                S.add.constdomainy(:,i) = result(1:end-1);
                S.add.surconst_standard(i) = result(end);
            case 'boxcox'
                value = [S.add.original_constdomainy(:,i);0]-1;
                sign = -1.*(value < 0) + 1.*(value >= 0);
                result = sign.*boxcox(abs(value));
                S.add.constdomainy(:,i) = result(1:end-1);
                S.add.surconst_standard(i) = result(end);
%                 S.add.constdomainy(:,i) = boxcox(S.add.original_constdomainy(:,i));
            case 'yeojohnson'
                value = [S.add.original_constdomainy(:,i);0]-1;
                sign = -1.*(value < 0) + 1.*(value >= 0);
                result = sign.*palm_yeojohnson(abs(value));
                S.add.constdomainy(:,i) = result(1:end-1);
                S.add.surconst_standard(i) = result(end);
%                 S.add.constdomainy(:,i) = palm_yeojohnson(S.add.original_constdomainy(:,i));
            case 'none'
                S.add.constdomainy(:,i) = S.add.original_constdomainy(:,i);
            otherwise
                A = -min(S.add.original_constdomainy(:,i))+1;
                B = exp(log(median(S.add.original_constdomainy(:,i)+A))/S.add.log_domscale);
                if strcmp(S.add.scaling,'log2')
                    Reverse_Modified_Objective = S.add.log_domscale - log([S.add.original_constdomainy(:,i);0]+A)/log(B);
                    Reverse_Modified_Objective(Reverse_Modified_Objective < 0) = 0;
                    value = S.add.log_domscale-log(Reverse_Modified_Objective+1)/log(S.add.C);
                    S.add.constdomainy(:,i) = value(1:end-1);
                    S.add.surconst_standard(i) = value(end);
                else
                    value = log([S.add.original_constdomainy(:,i);0]+A-1+B^30)/log(B)-30;
                    S.add.constdomainy(:,i) = value(1:end-1);
                    S.add.surconst_standard(i) = value(end);
                end
        end
    end
end

end