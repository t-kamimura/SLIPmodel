# -*- coding: utf-8 -*-

import numpy as np # 数値計算ライブラリの読み込み
import matplotlib.pyplot as plt # 描画ライブラリの読み込み

def eom(z, phase, x0):
    # 運動方程式
    x = z[0]
    y = z[1]
    u = z[2]
    v = z[3]
    dz = np.array([u, v, 0.0, 0.0])
    if phase == 0:
        # flight phase
        dz[2] = 0.0
        dz[3] = -g
    elif phase == 1:
        # stance phase
        l = np.sqrt((x - x0)**2 + y**2)
        dz[2] = -k/m * (l - l0) * (x - x0) / l
        dz[3] = -k/m * (l - l0) * y / l - g

    return dz

def eventHandler(z, phase, x0):
    if phase == 0:
        # flight phase
        toe = z[1] - l0 * np.cos(gamma) # 脚先高さ
        if toe < 0.0: # 接地判定
            print("touch down")
            x0 = z[0] + l0 * np.sin(gamma) # 接地点の更新
            phase = 1
    elif phase == 1:
        # stance phase
        l = np.sqrt((z[0] - x0)**2 + z[1]**2)
        extension = l - l0 # 脚ばねの伸び．負の値なら圧縮
        if extension > 0.0: # 離地判定
            print("lift off")
            phase = 0

    if z[1] < 0.0: # 転倒判定
        print("fall down")
        phase = -1

    return phase, x0

def RK4(z0, dt, tend):
    t = 0.0
    z = np.copy(z0) # y=y0だと、同一アドレスが共有されてしまう
    x0 = 0.0 # flightからはじめる

    tout = np.copy(t)
    zout = np.copy(z0)
    phase  = 0 # 0: flight, 1: stance
    
    while t < tend:
        t += dt
        k1 = eom(z, phase, x0)
        k2 = eom(z + 0.5*k1*dt, phase, x0)
        k3 = eom(z + 0.5*k2*dt, phase, x0)
        k4 = eom(z + k3*dt, phase, x0)
        z += (k1 + 2*k2 + 2*k3 + k4)*dt/6

        # イベント判定
        phase, x0 = eventHandler(z, phase, x0)
        if phase == -1:
            break # 転倒したら終了

        tout = np.vstack((tout, t)) # 結果を積み上げる
        zout = np.vstack((zout, z)) # 結果を積み上げる

    return tout, zout

# パラメータ設定
m = 80.0
k = 20*1e3
l0 = 1.0
g = 9.81

# 初期条件設定
y0 = 1.0
v0 = 7.0
q0 = np.array([0.0, y0, v0, 0.0]) # [x, y, u, v]
gamma = np.deg2rad(27.0) # 脚の角度

dt = 1e-2
tend = 5.0

# シミュレーション実行
tout, zout = RK4(q0, dt, tend)

# 描画
plt.figure()
plt.plot(tout, zout)
plt.legend(['x', 'y', 'dx', 'dy'])
plt.xlabel('time [s]')
plt.show()

plt.figure()
plt.plot(zout[:,0], zout[:,1])
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()