function S = Print(S)

fprintf('\n Iteration : %g',S.add.cnt);
fprintf('\n Current minimum value : %g (%g)\n\n',S.add.minimum_Value(end,2),S.add.minimum_Value(end,1));

S.print.x(S.add.cnt) = S.add.cnt;
S.print.y_min(S.add.cnt) = norm([2.3331998;1.942837;-0.479546;4.387811;-0.632617;1.039775;1.600996]-[S.add.domain(:,S.add.objective==S.add.minimum_Value(end,2))],'inf');
S.print.y_current(S.add.cnt) = norm([2.3331998;1.942837;-0.479546;4.387811;-0.632617;1.039775;1.600996]-[S.add.domain(:,end)]);
subplot(1,2,1); title('Min'); addpoints(animatedline,S.print.x,S.print.y_min); drawnow;
subplot(1,2,2); title('Current'); addpoints(animatedline,S.print.x,S.print.y_current); drawnow;

if S.print.y_min(S.add.cnt) < S.print.near(end,end) && S.add.constraint(end) == 0
    S.print.near(end+1,:) = [S.add.domain(:,end)', S.add.objective(end),S.print.y_min(S.add.cnt)];
end

end