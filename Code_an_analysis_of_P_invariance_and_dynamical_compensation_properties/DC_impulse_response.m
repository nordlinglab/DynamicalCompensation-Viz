%% Noisy input and disturbance
% System:
%     dy/dt = by(t)+d(t)-sz(t)(lr(t)-y(t))
%     dz/dt = -cz(t)(r(t)-y(t))
clear all, close all, clc
color = cbrewer2('qual', 'Set1', 5);
rng(1);
%% Create input and disturbance
% Create noisy input
r = 11; sys_d = 0.01;

input_shift = 50;
r_step = [ones(1, input_shift)*r, ones(1, 50)*r, ones(1, 250)*(r+0.5*r), ones(1, 50)*(r+0.25*r)];
r_noise = [zeros(1, input_shift), randn(1, 50), randn(1, 50), zeros(1, 150), randn(1, 50), zeros(1, 50)];
r_input = r_step + r_noise;
d_step = [ones(1, input_shift)*sys_d, ones(1, 200)*sys_d, ones(1, 100)*10, ones(1, 50)*5];
d_noise = [zeros(1, input_shift), zeros(1, 150), randn(1, 150), zeros(1, 50)];
d_input = d_step + d_noise;
d_input(d_input < 0) = 0;
tf = length(d_input);
t = 0:0.1:tf;
%% Plot input
% Plot input
figure(200)
plot_input(r_input, d_input, t)
xlabel('Time', 'Fontsize', 12, 'Interpreter', 'Latex')

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [9.5, 9.5])
saveas(gcf,'Akram_2023_DC_noisy_step_input.pdf')
%% Case 1 changing s
sys_a = 0;  sys_l = 0.7; sys_s = 0.25; sys_b = 0.3; sys_c = 2;
parm1 = [sys_a, sys_l, sys_s, sys_b, sys_c];
x0 = [r; -(sys_b+sys_d/r)*(1/(sys_s*(sys_l-1)))];
x0_shift = [10, 6];

% Plot input
figure(1)
plot_input(r_input, d_input, t)
% Plot response comparison: changing s
sys_a = 0; sys_l = 0.7; sys_s = 0.25*6; sys_b = 0.3; sys_c = 2;
parm2 = [sys_a, sys_l, sys_s, sys_b, sys_c];
plot_response(parm1, parm2, x0_shift, tf, r_input, d_input, color)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 12], 'PaperUnits', 'Inches', 'PaperSize', [10-1.5, 12-1.5])
saveas(gcf,'Akram_2023_DC_noisy_step_input_change_s.pdf')
%% Case 1 changing b
sys_a = 0;  sys_l = 0.7; sys_s = 0.25; sys_b = 0.3; sys_c = 2;
parm1 = [sys_a, sys_l, sys_s, sys_b, sys_c];
x0 = [r; -(sys_b+sys_d/r)*(1/(sys_s*(sys_l-1)))];

% Plot input
figure(2)
plot_input(r_input, d_input, t)
% Plot response comparison: changing s
sys_a = 0; sys_d = 0.01; sys_l = 0.7; sys_s = 0.25; sys_b = 0.3*2; sys_c = 2;
parm2 = [sys_a, sys_l, sys_s, sys_b, sys_c];
plot_response(parm1, parm2, x0_shift, tf, r_input, d_input, color)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 12], 'PaperUnits', 'Inches', 'PaperSize', [10-1.5, 12-1.5])
saveas(gcf,'Akram_2023_DC_noisy_step_input_change_b.pdf')
%% Case 3 changing c
sys_a = 0;  sys_l = 0.7; sys_s = 0.25; sys_b = 0.3; sys_c = 2;
parm1 = [sys_a, sys_l, sys_s, sys_b, sys_c];
x0 = [r; -(sys_b+sys_d/r)*(1/(sys_s*(sys_l-1)))];

% Plot input
figure(3)
plot_input(r_input, d_input, t)
% Plot response comparison: changing s
sys_a = 0; sys_d = 0.01; sys_l = 0.7; sys_s = 0.25; sys_b = 0.3; sys_c = 2*2;
parm2 = [sys_a, sys_l, sys_s, sys_b, sys_c];
plot_response(parm1, parm2, x0_shift, tf, r_input, d_input, color)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 12], 'PaperUnits', 'Inches', 'PaperSize', [10-1.5, 12-1.5])
saveas(gcf,'Akram_2023_DC_noisy_step_input_change_c.pdf')
%% Plot input with noise
function plot_input(r_input, d_input, t)
    % Plot input
    subplot(4,1,1)
    u = interp1([0:length(r_input)-1], r_input, t);
    plot(t, u, '-k', 'LineWidth', 2), hold on
    ylabel('$r(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
    
    subplot(4,1,2)
    d = interp1([0:length(d_input)-1], d_input, t);
    plot(t, d, '-k', 'LineWidth', 2), hold on
    % xlabel('Time', 'Fontsize', 12, 'Interpreter', 'Latex')
    ylabel('$d(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
end
%% Plot state reponses of y and z
function plot_response(parm1, parm2, x0_shift, tf, r_input, d_input, color)
    sys_a = parm1(1); sys_l = parm1(2); sys_s = parm1(3); sys_b = parm1(4); sys_c = parm1(5);
    sys_b_case1 = sys_b; sys_s_case1 = sys_s; sys_c_case1 = sys_c;
    % ODE solver
    my_opts = odeset('RelTol',1e-8,'AbsTol',1e-10);  
    f = @(t,x)DC_2D_generalized_model(t,x,r_input,...
        sys_a,sys_b,d_input,sys_s,sys_l,sys_c); %[y; z]
    [t,y] = ode45(f, [0:0.8:tf], x0_shift, my_opts);
    y_old = y;
    
    subplot(4,1,3)
    p1 = plot(t, y(:,1), '-', 'LineWidth', 2, 'Color', color(2,:)), hold on
    ylabel('$y(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
    
    subplot(4,1,4)
    plot(t, y(:,2), '-', 'LineWidth', 2, 'Color', color(2,:)), hold on
    xlabel('Time', 'Fontsize', 12, 'Interpreter', 'Latex')
    ylabel('$z(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
    hold on

    sys_a = parm2(1); sys_l = parm2(2); sys_s = parm2(3); sys_b = parm2(4); sys_c = parm2(5);
    sys_b_case2 = sys_b; sys_s_case2 = sys_s; sys_c_case2 = sys_c;
    % ODE solver
    my_opts = odeset('RelTol',1e-8,'AbsTol',1e-10);  
    f = @(t,x)DC_2D_generalized_model(t,x,r_input,...
        sys_a,sys_b,d_input,sys_s,sys_l,sys_c); %[y; z]
    [t,y] = ode45(f, [0:0.8:tf], x0_shift, my_opts);
    difference = y-y_old;
    
    subplot(4,1,3)
    p2 = plot(t, y(:,1), '--', 'LineWidth', 2, 'Color', color(1,:)), hold on
    p3 = plot(t, difference(:,1), '-.', 'LineWidth', 2, 'Color', color(3,:))
    ylabel('$y(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
    ylim([-7, 24])
    legend([p3], 'Difference')

    subplot(4,1,4)
    plot(t, y(:,2), '--', 'LineWidth', 2, 'Color', color(1,:)), hold on
    xlabel('Time', 'Fontsize', 12, 'Interpreter', 'Latex')
    ylabel('$z(t)$', 'Fontsize', 12, 'Interpreter', 'Latex')
    set(gca,'FontSize',12)
    grid on
    if sys_b_case1 ~= sys_b_case2
        legend(['b=', num2str(sys_b_case1)], ['b=', num2str(sys_b_case2)])
    elseif sys_s_case1 ~= sys_s_case2
        legend(['s=', num2str(sys_s_case1)], ['s=', num2str(sys_s_case2)])
    elseif sys_c_case1 ~= sys_c_case2
        legend(['c=', num2str(sys_c_case1)], ['c=', num2str(sys_c_case2)])
    else
        disp("Parameter setting error")
    end
    
    ylim([0, 20])
    
end