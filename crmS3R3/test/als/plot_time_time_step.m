function plot_time_time_step(sols)

[m,n] = size(sols);

if ~isfield(sols{1}, 'h')
    sols = create_h(sols);
end

the_size = max(size(sols));
colors = linspecer(the_size+1);

figure()
for i = 1:max(m,n)
    N_h = max(size(sols{i}.h));
    N_t = max(size(sols{i}.rslt.t));
    if N_h+1 ~= N_t
        sols{i}.h = sols{i}.h * ones(N_t-1,1);
    end
    semilogy(sols{i}.rslt.t,[0; sols{i}.h],'LineWidth',1.5,'DisplayName',['N=' num2str(sols{i}.n)],'Color',colors(i,:))
    hold on
end

legend('Location','best','Interpreter','latex','FontSize',14)
title('\bf Time step size vs Time','Interpreter','latex','FontSize',18)
xlabel('Time ($s$)','Interpreter','latex','FontSize',18)
ylabel('Time step size ($h_{new}$)','Interpreter','latex','FontSize',18)
grid on
axis square
box on

end