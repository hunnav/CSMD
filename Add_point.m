function S = Add_point(S)

S.add.cnt = S.add.cnt + 1;
n = S.prob.numinitsam + S.add.cnt;
S.add.domain(:,n) = x;
S.add.objective(n,1) = f(x(1),x(2),x(3),x(4),x(5),x(6),x(7));    % we need to change it according to the conditon

gg1 = 0-S.prob.c1(x(1),x(2),x(3),x(4),x(5),x(6),x(7));            % we need to change it according to the conditon
gg2 = 0-S.prob.c2(x(1),x(2),x(3),x(4),x(5),x(6),x(7));            % we need to change it according to the conditon
gg3 = 0-S.prob.c3(x(1),x(2),x(3),x(4),x(5),x(6),x(7));            % we need to change it according to the conditon
gg4 = 0-S.prob.c4(x(1),x(2),x(3),x(4),x(5),x(6),x(7));            % we need to change it according to the conditon
S.add.constraint(n,1) = max([gg1,gg2,gg3,gg4,0]);

S = Scaling(S);

if (S.add.add(end) == 0 && S.add.objective(end) < S.add.minimum_Value(end,2))
    S.add.minimum_Value(end+1,:) = [S.add.cnt, S.add.objective(end)];
end

end
