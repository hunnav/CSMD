function S = Add_point(S)

S.add.cnt = S.add.cnt + 1;
n = S.prob.numinitsam + S.add.cnt;
S.add.domain(:,n) = S.acqui.x;
S.add.objective(n,1) = S.prob.f(S.acqui.x(1),S.acqui.x(2),S.acqui.x(3),S.acqui.x(4),S.acqui.x(5),S.acqui.x(6),S.acqui.x(7));    % we need to change it according to the conditon

p0 = 10^(floor(-1+log10(abs(S.add.objective(n,1)))));

gg1 = 0-S.prob.c1(S.acqui.x(1),S.acqui.x(2),S.acqui.x(3),S.acqui.x(4),S.acqui.x(5),S.acqui.x(6),S.acqui.x(7));              % we need to change it according to the conditon
gg2 = 0-S.prob.c2(S.acqui.x(1),S.acqui.x(2),S.acqui.x(3),S.acqui.x(4),S.acqui.x(5),S.acqui.x(6),S.acqui.x(7));              % we need to change it according to the conditon
gg3 = 0-S.prob.c3(S.acqui.x(1),S.acqui.x(2),S.acqui.x(3),S.acqui.x(4),S.acqui.x(5),S.acqui.x(6),S.acqui.x(7));              % we need to change it according to the conditon
gg4 = 0-S.prob.c4(S.acqui.x(1),S.acqui.x(2),S.acqui.x(3),S.acqui.x(4),S.acqui.x(5),S.acqui.x(6),S.acqui.x(7));              % we need to change it according to the conditon
if S.prob.surconst == 1
    S.add.original_constdomainy(n,:) = [gg1,gg2,gg3,gg4];                                                                   % we need to change it according to the conditon
end
S.add.constraint(n,1) = p0*max([gg1,gg2,gg3,gg4,0]);                                                                        % we need to change it according to the conditon
S.add.obPLUScon(n,1) = S.add.objective(n,1) + S.add.constraint(n,1);

S = Scaling(S,n);

if (S.add.constraint(end) == 0 && S.add.objective(end) < S.add.minimum_Value(end,2))
    S.add.minimum_Value(end+1,1:2) = [S.add.cnt, S.add.objective(end)];
end

end
