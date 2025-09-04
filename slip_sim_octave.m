% SLIPモデル

clear
close all

% parameter setting
param.m = 80.0;
param.k = 20*1e3;
param.l0 = 1.0;
param.g = 9.81;
param.x0 = 0.0;
param.gamma = 0.0;
param.phaseflag = 0;


function du = eom(u, p)
    % SLIP model dynamics
    if p.phaseflag == 0
        % flight phase dynamics
        % disp("flight phase")
        du(1) = u(3);
        du(2) = u(4);
        du(3) = 0.0;
        du(4) = -p.g;
    else
        % stance phase dynamics
        % disp("stance phase")
        x = u(1);
        y = u(2);
        l = sqrt((x - p.x0)^2 + y^2);
        du(1) = u(3);
        du(2) = u(4);
        du(3) = -p.k/p.m * (l - p.l0) * (x - p.x0) / l;
        du(4) = -p.k/p.m * (l - p.l0) * y / l - p.g;
    end
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

        x = u(1);
        y = u(2);
        gamma = p.gamma;
        l0 = p.l0;
        % event function
        if p.phaseflag == 0
            % flight phase -> stance phase
            out(1) = y - l0*cos(gamma);
            out(2) = y;
            if out(1) < 0
                % flight phase -> stance phase
                p.phaseflag = 1;
                p.x0 = u(1) + l0*sin(gamma);
                disp("touch down")
            end
        else
            % stance phase -> flight phase
            x0 = p.x0;
            l = sqrt((x-x0)^2 + y^2);
            out(1) = l - l0;
            out(2) = y;
            if out(1) > 0
                % stance phase -> flight phase
                p.phaseflag = 0;
                p.x0 = u(1);
                disp("lift off")
            end
        end
        if out(2)<0
            disp("fall down")
            break
        end
    end
end

y0 = 1.0;
v0 = 7.0;
u0 = [0.0, y0, v0, 0.0];
param.gamma = pi/6;
dt = 1e-2;
tspan = [0.0, 1.0];

[tout, uout] = RK4(u0, param, dt, tspan);

% plot
figure;
plot(uout(:,1), uout(:,2), 'b-');
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title('SLIP model simulation');

