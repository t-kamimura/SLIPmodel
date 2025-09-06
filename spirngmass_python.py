# -*- coding: utf-8 -*-

import numpy as np # 数値計算ライブラリの読み込み
import matplotlib.pyplot as plt # 描画ライブラリの読み込み


def eom(y):
    # 運動方程式
    x = y[0]
    v = y[1]
    dy = np.array([v, (-k*x - c*v)/m])
    return dy

def Eular(y0, dt, tend):
    # オイラー法
    t = 0.0
    y = np.copy(y0) # y=y0だと、同一アドレスが共有されてしまう
    tout = np.copy(t)
    yout = np.copy(y)
    while t < tend:
        t += dt
        y += eom(y)*dt
        tout = np.vstack((tout, t)) # 結果を積み上げる
        yout = np.vstack((yout, y)) # 結果を積み上げる
    return tout, yout

def RK4(y0, dt, tend):
    t = 0.0
    y = np.copy(y0) # y=y0だと、同一アドレスが共有されてしまう
    tout = np.copy(t)
    yout = np.copy(y)
    while t < tend:
        t += dt
        k1 = eom(y)
        k2 = eom(y + 0.5*k1*dt)
        k3 = eom(y + 0.5*k2*dt)
        k4 = eom(y + k3*dt)
        y += (k1 + 2*k2 + 2*k3 + k4)*dt/6
        tout = np.vstack((tout, t)) # 結果を積み上げる
        yout = np.vstack((yout, y)) # 結果を積み上げる
    return tout, yout
    
# パラメータの定義
m = 1.0
k = 100.0
c = 10.0

# 初期値、刻み幅、シミュレーション時間の定義
y0 = np.array([1.0, 0.0])
dt = 1e-2
tend = 3.0

# 数値積分の実行
# tout, yout = Eular(y0,dt,tend)
tout, yout = RK4(y0, dt, tend)

# 結果の描画
plt.figure()
plt.plot(tout,yout)
plt.xlabel("t")
plt.ylabel("x,v")
plt.legend(("x","v"), loc="best")
plt.show()