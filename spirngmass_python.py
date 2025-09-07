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
    
def exact(y0,tout):
    # 厳密解
    omega = np.sqrt(k/m - (c/(2*m))**2)
    A = np.copy(y0[0])
    B = (y0[1] + c/(2*m)*y0[0]) / omega
    xout = np.zeros(len(tout))
    for i in range(len(tout)):
        t = tout[i].item() # tout[i]はnumpyのスカラー型なので、.item()でPythonのスカラー型に変換
        xout[i] = np.exp(-c/(2*m)*t) * (A * np.cos(omega*t) + B * np.sin(omega*t))
    return xout

# パラメータの定義
m = 1.0
k = 100.0
c = 10.0

# 初期値、刻み幅、シミュレーション時間の定義
y0 = np.array([1.0, 0.0])
dt = 5*1e-2
tend = 3.0

# 数値積分の実行
# tout, yout = Eular(y0,dt,tend)
tout, yout = RK4(y0, dt, tend)
x_exact = exact(y0, tout)

# 結果の描画
plt.figure()
plt.plot(tout,yout)
plt.plot(tout,x_exact,"k--")
plt.xlabel("t")
plt.ylabel("x,v")
plt.legend(("x","v"), loc="best")
plt.show()