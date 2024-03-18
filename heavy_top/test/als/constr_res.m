function constr_res(sols)

size_s = max(size(sols));

figure()
subplot(1,2,2)
hold on
subplot(1,2,1)
hold on

for i = 1:size_s
    if floor(sols{i}.problemset/10) == 6
        subplot(1,2,2)
    else
        subplot(1,2,1)
    end
    plot(sols{i}.rslt.t,vecnorm(sols{i}.rslt.Phi),'DisplayName',['a = ' num2str(sols{i}.a_baumgarte)],'LineWidth',1.5)
end

subplot(1,2,2)
title('SE3')
legend()

subplot(1,2,1)
title('SO3xR3')
legend()

end