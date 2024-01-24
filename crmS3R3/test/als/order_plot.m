function order_plot(sol2, sol1)

h = sol2(:, 1);

% order 2
max_RTR_2 = sol2(:,2);
max_res_2 = sol2(:,3);
err_x_2 = sol2(:,4);
err_R_2 = sol2(:,5);
err_u_2 = sol2(:,6);
err_O_2 = sol2(:,7);
err_l_2 = sol2(:,8);


% order 4
max_RTR_4 = sol2(:,9);
max_res_4 = sol2(:,10);
err_x_4 = sol2(:,11);
err_R_4 = sol2(:,12);
err_u_4 = sol2(:,13);
err_O_4 = sol2(:,14);
err_l_4 = sol2(:,15);


% order 5
max_RTR_5 = sol2(:,16);
max_res_5 = sol2(:,17);
err_x_5 = sol2(:,18);
err_R_5 = sol2(:,19);
err_u_5 = sol2(:,20);
err_O_5 = sol2(:,21);
err_l_5 = sol2(:,22);

colors = linspecer(4);


h_se = sol1(:, 1);

% order 2
max_RTR_2_se = sol1(:,2);
max_res_2_se = sol1(:,3);
err_x_2_se = sol1(:,4);
err_R_2_se = sol1(:,5);
err_u_2_se = sol1(:,6);
err_O_2_se = sol1(:,7);
err_l_2_se = sol1(:,8);


% order 4
max_RTR_4_se = sol1(:,9);
max_res_4_se = sol1(:,10);
err_x_4_se = sol1(:,11);
err_R_4_se = sol1(:,12);
err_u_4_se = sol1(:,13);
err_O_4_se = sol1(:,14);
err_l_4_se = sol1(:,15);


% order 5
max_RTR_5_se = sol1(:,16);
max_res_5_se = sol1(:,17);
err_x_5_se = sol1(:,18);
err_R_5_se = sol1(:,19);
err_u_5_se = sol1(:,20);
err_O_5_se = sol1(:,21);
err_l_5_se = sol1(:,22);



figure()
subplot(1,2,1)
loglog(h_se, err_x_2_se, 'LineWidth', 1.5, 'Color', colors(1,:), 'DisplayName', 'p=2')
hold on
loglog(h_se, err_x_4_se, 'LineWidth', 1.5, 'Color', colors(2,:), 'DisplayName', 'p=4')
loglog(h_se, err_x_5_se, 'LineWidth', 1.5, 'Color', colors(3,:), 'DisplayName', 'p=5')
title('Convergence order at position level', 'FontSize', 16)
xlabel('Time step size (h)', 'FontSize', 14)
ylabel('Global error', 'FontSize', 14)
grid on
legend('Location','northwest')

subplot(1,2,2)
loglog(h, max_res_5, 'LineWidth', 1.5, 'Color', colors(1,:), 'DisplayName', 'SO(3)xR^3')
hold on
loglog(h_se, max_res_5_se, 'LineWidth', 1.5, 'Color', colors(2,:), 'DisplayName', 'SE(3)')
title('Drift-off effect', 'FontSize', 16)
xlabel('Time step size (h)', 'FontSize', 14)
ylabel('||\Phi(R,x)||_2', 'FontSize', 14)
grid on
legend('Location','northwest')

end