clc
for i=1:8
    [crit_cal(i),ind] = max(sub(i).crit_cal);
    chan{i} = sub(i).chan{ind};
    t_start(i) = sub(i).t_start(ind);
    t_end(i) = sub(i).t_end(ind);
    disp(['individu : ' num2str(i)])
    disp(['chan : ' num2str(chan{i})])
    disp(['interval de temp : ' num2str([t_start(i) t_end(i)])])
    disp(['pwc :' num2str(crit_cal(i))])
    disp(' ')
end

disp(['moyenne global : ' num2str(mean(crit_cal))])