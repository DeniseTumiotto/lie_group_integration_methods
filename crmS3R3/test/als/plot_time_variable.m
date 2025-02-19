function plot_time_variable(t_h)

n = max(size(t_h));

figure()
hold on
plot(t_h(:,1),t_h(:,2),'k-', 'LineWidth',0.5)
for i = 1:n
    if t_h(i,3)>1
        plot(t_h(i,1),t_h(i,2),'rx','MarkerSize',5)
    elseif i>1 && t_h(i-1,3)>1
        plot(t_h(i,1),t_h(i,2),'go','MarkerSize',5)
    end
end
end