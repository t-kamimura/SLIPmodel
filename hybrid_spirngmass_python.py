# -*- coding: utf-8 -*-

import numpy as np # 数値計算ライブラリの読み込み
import matplotlib.pyplot as plt # 描画ライブラリの読み込み

# パラメータの定義
m = 1.0
k = 10
c = 0.4
kappa = 0.45
epsilon = 0.15
omega0 = np.sqrt(k/m)
omega = omega0/kappa

def eom(t,y):
    # 運動方程式
    x = y[0]
    v = y[1]
    if (t*omega//np.pi) % 2 == 1:
        omega_ = omega*(1 + epsilon)
    else:
        omega_ = omega*(1 - epsilon)
    k_ = m*omega_**2
    dy = np.array([v, (-k_*x - c*v)/m])
    return dy

def RK4(y0, dt, tend):
    t = 0.0
    y = np.copy(y0) # y=y0だと、同一アドレスが共有されてしまう
    tout = np.copy(t)
    yout = np.copy(y)
    while t < tend:
        t += dt
        k1 = eom(t,y)
        k2 = eom(t,y + 0.5*k1*dt)
        k3 = eom(t,y + 0.5*k2*dt)
        k4 = eom(t,y + k3*dt)
        y += (k1 + 2*k2 + 2*k3 + k4)*dt/6
        tout = np.vstack((tout, t)) # 結果を積み上げる
        yout = np.vstack((yout, y)) # 結果を積み上げる
    return tout, yout


# 以下，メイン処理
# 初期値、刻み幅、シミュレーション時間の定義
y0 = np.array([0.1, 0.0])
dt = 1e-2
tend = 10.0

# 数値積分の実行
tout, yout = RK4(y0, dt, tend)
# x_exact = exact(y0, tout)

# 結果の描画
plt.figure()
plt.plot(tout,yout)
# plt.plot(tout,x_exact,"k--")
plt.xlabel("t")
plt.ylabel("x,v")
plt.legend(("x","v"), loc="best")
plt.show()