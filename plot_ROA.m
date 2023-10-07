function a = plot_ROA(f_x, rad)
%Plotting
    figure
    min = -1*rad;
    max = rad;
    step = 0.5;
    tspan = [0, 10];
    x0 = [0; 0];
    [t_plot, x_plot] = ode45(f_x, tspan, x0);
    x_1_plt = x_plot(:, 1);
    x_2_plt = x_plot(:, 2);
    plot(x_1_plt, x_2_plt, '-c')
    hold on
    grid
    for i = min:step:max
        for j = min:step:max
            x0 = [i; j];
            [t_plot, x_plot] = ode45(f_x, tspan, x0);
            x_1_plt = x_plot(:, 1);
            x_2_plt = x_plot(:, 2);
    %         plot(x_1_plt, x_2_plt)
            if x_1_plt(end) < 1 && x_1_plt(end) > -1 && x_2_plt(end) < 1 && x_2_plt(end) > -1
                plot(x0(1), x0(2), 'r--o')
                plot(x_1_plt, x_2_plt, '-c')
            end
        end
    end
%     tspan = [0, 10];
%     x0 = [-2; 0];
%     [t_plot, x_plot] = ode45(f_x, tspan, x0);
%     x_1_plt = x_plot(:, 1);
%     x_2_plt = x_plot(:, 2);
%     plot(x_1_plt, x_2_plt)
%     hold on
%     grid
%     x0 = [2; 0];
%     [t_plot, x_plot] = ode45(f_x, tspan, x0);
%     x_1_plt = x_plot(:, 1);
%     x_2_plt = x_plot(:, 2);
%     plot(x_1_plt, x_2_plt)
end