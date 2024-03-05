function order_plot_ht(sol2, sol1)

h = sol2(:, 1);

% order 2
max_RTR_2 = sol2(:,2);
max_res_2 = sol2(:,3);
err_x_2 = sol2(:,4);
err_R_2 = sol2(:,5);
err_u_2 = sol2(:,6);
err_O_2 = sol2(:,7);
err_l_2 = sol2(:,8);


% order 3
max_RTR_3 = sol2(:,9);
max_res_3 = sol2(:,10);
err_x_3 = sol2(:,11);
err_R_3 = sol2(:,12);
err_u_3 = sol2(:,13);
err_O_3 = sol2(:,14);
err_l_3 = sol2(:,15);


% order 4
max_RTR_4 = sol2(:,16);
max_res_4 = sol2(:,17);
err_x_4 = sol2(:,18);
err_R_4 = sol2(:,19);
err_u_4 = sol2(:,20);
err_O_4 = sol2(:,21);
err_l_4 = sol2(:,22);


% order 5
max_RTR_5 = sol2(:,23);
max_res_5 = sol2(:,24);
err_x_5 = sol2(:,25);
err_R_5 = sol2(:,26);
err_u_5 = sol2(:,27);
err_O_5 = sol2(:,28);
err_l_5 = sol2(:,29);

colors = linspecer(5);


h_se = sol1(:, 1);

% order 2
max_RTR_2_se = sol1(:,2);
max_res_2_se = sol1(:,3);
err_x_2_se = sol1(:,4);
err_R_2_se = sol1(:,5);
err_u_2_se = sol1(:,6);
err_O_2_se = sol1(:,7);
err_l_2_se = sol1(:,8);


% order 3
max_RTR_3_se = sol1(:,9);
max_res_3_se = sol1(:,10);
err_x_3_se = sol1(:,11);
err_R_3_se = sol1(:,12);
err_u_3_se = sol1(:,13);
err_O_3_se = sol1(:,14);
err_l_3_se = sol1(:,15);


% order 4
max_RTR_4_se = sol1(:,16);
max_res_4_se = sol1(:,17);
err_x_4_se = sol1(:,18);
err_R_4_se = sol1(:,19);
err_u_4_se = sol1(:,20);
err_O_4_se = sol1(:,21);
err_l_4_se = sol1(:,22);

% order 5
max_RTR_5_se = sol1(:,23);
max_res_5_se = sol1(:,24);
err_x_5_se = sol1(:,25);
err_R_5_se = sol1(:,26);
err_u_5_se = sol1(:,27);
err_O_5_se = sol1(:,28);
err_l_5_se = sol1(:,29);


figure()
subplot(1,2,1)
loglog(h_se, err_x_2_se, 'LineWidth', 1.5, 'Color', colors(1,:), 'DisplayName', 'p=2')
hold on
loglog(h_se, err_x_3_se, 'LineWidth', 1.5, 'Color', colors(2,:), 'DisplayName', 'p=3')
loglog(h_se, err_x_4_se, 'LineWidth', 1.5, 'Color', colors(3,:), 'DisplayName', 'p=4')
loglog(h_se, err_x_5_se, 'LineWidth', 1.5, 'Color', colors(4,:), 'DisplayName', 'p=5')
title('Convergence order: q', 'FontSize', 16)
xlabel('Time step size (h)', 'FontSize', 14)
ylabel('Global error', 'FontSize', 14)
grid on
legend('Location','northwest')
axis([10^-4 10^-2 10^-10 10^0])
axis square;

subplot(1,2,2)
loglog(h_se, err_l_2_se, 'LineWidth', 1.5, 'Color', colors(1,:), 'DisplayName', 'p=2')
hold on
loglog(h_se, err_l_3_se, 'LineWidth', 1.5, 'Color', colors(2,:), 'DisplayName', 'p=3')
loglog(h_se, err_l_4_se, 'LineWidth', 1.5, 'Color', colors(3,:), 'DisplayName', 'p=4')
loglog(h_se, err_l_5_se, 'LineWidth', 1.5, 'Color', colors(4,:), 'DisplayName', 'p=5')
title('Convergence order: \lambda', 'FontSize', 16)
xlabel('Time step size (h)', 'FontSize', 14)
ylabel('Global error', 'FontSize', 14)
grid on
legend('Location','northwest')
axis([10^-4 10^-2 10^-8 10^2])
axis square;

figure()
subplot(1,2,1)
loglog(h_se, err_x_2, 'LineWidth', 1.5, 'Color', colors(1,:), 'DisplayName', 'p=2')
hold on
loglog(h_se, err_x_3, 'LineWidth', 1.5, 'Color', colors(2,:), 'DisplayName', 'p=3')
loglog(h_se, err_x_4, 'LineWidth', 1.5, 'Color', colors(3,:), 'DisplayName', 'p=4')
loglog(h_se, err_x_5, 'LineWidth', 1.5, 'Color', colors(4,:), 'DisplayName', 'p=5')
title('Convergence order: q', 'FontSize', 16)
xlabel('Time step size (h)', 'FontSize', 14)
ylabel('Global error', 'FontSize', 14)
grid on
legend('Location','northwest')
axis([10^-4 10^-2 10^-10 10^0])
axis square;

subplot(1,2,2)
loglog(h_se, err_l_2, 'LineWidth', 1.5, 'Color', colors(1,:), 'DisplayName', 'p=2')
hold on
loglog(h_se, err_l_3, 'LineWidth', 1.5, 'Color', colors(2,:), 'DisplayName', 'p=3')
loglog(h_se, err_l_4, 'LineWidth', 1.5, 'Color', colors(3,:), 'DisplayName', 'p=4')
loglog(h_se, err_l_5, 'LineWidth', 1.5, 'Color', colors(4,:), 'DisplayName', 'p=5')
title('Convergence order: \lambda', 'FontSize', 16)
xlabel('Time step size (h)', 'FontSize', 14)
ylabel('Global error', 'FontSize', 14)
grid on
legend('Location','northwest')
axis([10^-4 10^-2 10^-8 10^2])
axis square;

figure()
subplot(1,2,1)
loglog(h, max_res_5, 'LineWidth', 1.5, 'Color', colors(1,:), 'DisplayName', 'SO(3)xR^3')
title('Drift-off effect in SO(3)xR^3', 'FontSize', 16)
xlabel('Time step size (h)', 'FontSize', 14)
ylabel('||\Phi(R,x)||_2', 'FontSize', 14)
grid on
axis([10^-4 10^-2 10^-15 10^0])
axis square;

subplot(1,2,2)
loglog(h_se, max_res_5_se, 'LineWidth', 1.5, 'Color', colors(2,:), 'DisplayName', 'SE(3)')
title('Drift-off effect in SE(3)', 'FontSize', 16)
xlabel('Time step size (h)', 'FontSize', 14)
ylabel('||\Phi(R,x)||_2', 'FontSize', 14)
grid on
axis([10^-4 10^-2 10^-15 10^0])
axis square;

end