except_num = [1,2];
R = ones(4,4);
idx = true(size(R));
idx(except_num,:) = false;
idx(:,except_num) = false;
idx

