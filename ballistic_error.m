% SLIPモデル

clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultAxesFontSize',16);
set(0,'defaultAxesFontName', 'Times new roman');
set(0,'defaultTextFontSize',16);
set(0,'defaultTextFontName', 'Times new roman');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% parameter setting
param.m = 1.0;
param.g = 9.81;

function du = eom(u, p)
    % flight phase dynamics
    % disp("flight phase")
    du(1) = u(3);
    du(2) = u(4);
    du(3) = 0.0;
    du(4) = -p.g;
end

function [tout, uout] =  RK4(u, p, dt, tspan)
    tout = [];
    uout = [];
    t = tspan(1);
    for i_t = 1:(tspan(2)-tspan(1))/dt
        k1 = eom(u, p);
        k2 = eom(u + 0.5 * dt * k1, p);
        k3 = eom(u + 0.5 * dt * k2, p);
        k4 = eom(u + dt * k3, p);
        u = u + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;

        tout = [tout; t];
        uout = [uout; u];
    end
end

function xout = calc_exact(t, u0, p)
    % exact solution during flight phase
    xout = zeros(length(t), 4);
    for i = 1:length(t)
        ti = t(i);
        xout(i,1) = u0(1) + u0(3)*ti;
        xout(i,2) = u0(2) + u0(4)*ti - 0.5*p.g*ti^2;
        xout(i,3) = u0(3);
        xout(i,4) = u0(4) - p.g*ti;
    end
end

y0 = 1.0;
v0 = 2.0;
u0 = [0.0, y0, v0, 0.0];
dtset = -1:-0.1:-4;
tspan = [0.0, 1.0];

% figure
for i_dt = 1:length(dtset)
    dt = 10^dtset(i_dt);
    [tout, uout] = RK4(u0, param, dt, tspan);
    x_exact = calc_exact(tout, u0, param);
    error_begin(i_dt) = abs(uout(2,2) - x_exact(2,2));
    error_end(i_dt) = abs(uout(end,2) - x_exact(end,2));

    % subplot(2,4,i_dt);
    % plot(tout, uout(:,2), '-', 'DisplayName', 'RK4 solution');
    % hold on;
    % plot(tout, x_exact(:,2), '--', 'DisplayName', 'Exact solution');
    % xlabel('Time (s)');
    % ylabel('Vertical position (m)');
    % title(sprintf('dt = %.1e', dt));
    % grid on;
    % legend('Location', 'best');
    % subplot(2,4,i_dt+4);
    % plot(tout, (uout(:,2) - x_exact(:,2)), '-', 'DisplayName', 'Position error');
    % xlabel('Time (s)');
    % ylabel('Position error (m)');
    % title(sprintf('dt = %.1e', dt));
    % grid on;
    % % legend('Location', 'best');
end
for i = 4:5
    for i_dt = 1:length(dtset)
        dt = 10^dtset(i_dt);
        refset(i, i_dt) = (tspan(end)/dt) * dt^i;
        % refset2(i, i_dt) = dt^i;
    end
end

% plot
figure;
dtset_=10.^dtset;
plot(dtset_, error_end, '+','DisplayName', 'RK4 error');
hold on;
% plot(dtset, error_begin, 'o','DisplayName', 'RK4 error (begin)');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time step size (s)');
ylabel('Position error (m)');
% grid on;
for i = 4:5
    plot(dtset_, refset(i,:), '--','DisplayName', sprintf('$$O(dt^{%d})$$', i));
end

legend('Location', 'best', 'Interpreter','latex');