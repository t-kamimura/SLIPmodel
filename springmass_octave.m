clear
close all

global m k c
m = 1.0;
k = 50.0;
c = 10.0;

function dy = eom(y)
  % 運動方程式
  global m k c
  x = y(1);
  v = y(2);
  dy(1) = v;
  dy(2) = (-k*x -c*v)/m;
end

function [tout, yout] = Eular(y0, dt, tend)
  % オイラー法
  t = 0.0;
  y = y0;
  tout(1) = t;
  yout(1,:) = y;

  while t < tend
    t = t + dt;
    y = y + eom(y)*dt;
    tout = [tout;t]; % 結果を積み上げる
    yout = [yout;y]; % 結果を積み上げる
  end
end

function [tout, yout] = RK4(y0, dt, tend)
  % 4次のルンゲクッタ法
  t = 0.0;
  y = y0;
  tout(1) = t;
  yout(1,:) = y;

  while t < tend
    t = t + dt;
    k1 = eom(y);
    k2 = eom(y + 0.5*k1*dt);
    k3 = eom(y + 0.5*k2*dt);
    k4 = eom(y + k3*dt);
    y = y + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
    tout = [tout;t]; % 結果を積み上げる
    yout = [yout;y]; % 結果を積み上げる
  end
end

function xout = analytical(tout, y0)
  % 解析解を計算する．ただしv0=0とする
  global m k c
  zeta = c/(2*m);
  w = sqrt(k/m - zeta^2);
  xout = y0(1)*exp(-zeta*tout).*(cos(w*tout) + (zeta/w)*sin(w*tout));
end

% 初期値、刻み幅、シミュレーション時間の定義
y0 = [1.0, 0.0];
dt = 5*1e-2;
tend = 0.5;

% 数値積分の実行
[tout, yout] = Eular(y0, dt, tend);
% [tout, yout] = RK4(y0, dt, tend);

% 結果のプロット
figure
plot(tout,yout,"o-","markersize",3)
xlabel("t","fontsize",12)
ylabel("x, v","fontsize",12)
legend({"x", "v"},"fontsize",12)

% 解析解の計算
x_analy = analytical(tout, y0);
hold on
plot(tout, x_analy, "r-","linewidth",1)
legend({"x (num)", "v (num)", "x (ana)"},"fontsize",12)

% 誤差の計算
dtset = [1e-4, 1e-3, 1e-2, 1e-1, 1];
figure
for i_dt = 1:length(dtset)
  tend = dtset(i_dt);
  [tout1, yout1] = Eular(y0, dtset(i_dt), tend);
  [tout2, yout2] = RK4(y0, dtset(i_dt), tend);
  x_analy = analytical(tout1, y0);
  error1(i_dt) = abs(yout1(end,1) - x_analy(end));
  error2(i_dt) = abs(yout2(end,1) - x_analy(end));
end
loglog(dtset, error1, "+r", "Displayname", "Eular")
hold on
loglog(dtset, error2, "+b", "Displayname", "RK4")
hold on
for i = 1:5
  orderset = dtset.^i;
  loglog(dtset, orderset, "Displayname", num2str(i))
end
legend()
