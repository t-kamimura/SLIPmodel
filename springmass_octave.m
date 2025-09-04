clear
close all

global m = 1.0
global k = 50.0
global c = 10.0

function dy = eom(y)
  global m k c
  x = y(1);
  v = y(2);
  dy(1) = v;
  dy(2) = (-k*x -c*v)/m;
end

function [tout, yout] = Eular(y0, dt, tspan)
  n = (tspan(end)-tspan(1))/dt;

  t = tspan(1);
  y = y0;
  tout(1) = t;
  yout(1,:) = y;

  for i = 2:n
    t = t + dt
    y = y + eom(y)*dt;
    tout = [tout;t];
    yout = [yout;y];
  endfor
end

function [tout, yout] = RK4(y0, dt, tspan)
  n = (tspan(end)-tspan(1))/dt;

  t = tspan(1);
  y = y0;
  tout(1) = t;
  yout(1,:) = y;

  for i = 2:n
    t = t + dt
    k1 = eom(y);
    k2 = eom(y + 0.5*k1*dt);
    k3 = eom(y + 0.5*k2*dt);
    k4 = eom(y + k3*dt);
    y = y + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
    tout = [tout;t];
    yout = [yout;y];
  endfor
end


y0 = [1.0, 0.0];
dt = 5*1e-2;
tspan = [0.0, 0.5];

[tout1, yout1] = Eular(y0, dt, tspan);
[tout2, yout2] = RK4(y0, dt, tspan);

figure
plot(tout1,yout1,"o-","markersize",3)
hold on
plot(tout2,yout2,"^-","markersize",3)

xlabel("t","fontsize",12)
ylabel("x, v","fontsize",12)
legend({"x(Eular)", "v(Eular)","x(RK4)","v(RK4)"},"fontsize",12)
