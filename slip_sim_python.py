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
            x0 = z[0] + z[1] * np.tan(gamma) # 接地点の更新(ここが数値例2のポイントかも)
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
    z = np.copy(z0) # z=z0だと、同一アドレスが共有されてしまう
    x0 = 0.0 # flightからはじめる

    tout = np.copy(t)
    zout = np.copy(z0)
    phase  = 0 # 0: flight, 1: stance
    pout = np.copy(phase)
    toeout = np.array([z[0]+l0*np.sin(gamma), z[1]-l0*np.cos(gamma)]) # 脚先位置
    
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
        pout = np.vstack((pout, phase)) # 結果を積み上げる
        if phase == 0:
            toe = np.array([z[0]+l0*np.sin(gamma), z[1]-l0*np.cos(gamma)]) # 脚先位置
        else:
            toe = np.array([x0, 0.0])
        toeout = np.vstack((toeout, toe)) # 脚先位置を積み上げる

    return tout, zout, pout, toeout

# パラメータ設定
m = 1.0
k = 200
l0 = 1.0
g = 9.81

# 初期条件設定
# # 数値例1
# y0 = 1.05
# v0 = 6.455
# q0 = np.array([0.0, y0, v0, 0.0]) # [x, y, u, v]
# gamma = 0.54 # 脚の角度

# 数値例2
y0 = 1.05
v0 = 4.4023
q0 = np.array([0.0, y0, v0, 0.0]) # [x, y, u, v]
gamma = 0.43017 # 脚の角度

dt = 1e-2
tend = 10.0

# シミュレーション実行
tout, zout, pout, toeout = RK4(q0, dt, tend)

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

# # アニメーションの作成
# import matplotlib.animation as animation
# def make_animation(tout, zout, toeout):
#     fig, ax = plt.subplots()
#     ax.set_xlim(0, 10.0)
#     ax.set_ylim(-0.1, 1.5)
#     ax.set_aspect('equal')
#     ax.set_title("SLIP model")
#     ax.axis('on')

#     # バネと質点の描画要素
#     spring_line, = ax.plot([], [], 'k-', lw=2)
#     mass_patch = plt.Circle((0, 0), 0.1, fc='k')
#     toe_patch = plt.Circle((0, 0), 0.01, fc='k')
#     ax.add_patch(mass_patch)
#     ax.add_patch(toe_patch)

#     def init():
#         mass_patch.center = (zout[0, 0], zout[0, 1])
#         toe_patch.center = (toeout[0, 0], toeout[0, 1])
#         spring_line.set_data([zout[0, 0], toeout[0, 0]], [zout[0, 1], toeout[0, 1]])
#         return mass_patch, toe_patch, spring_line

#     def update(frame):
#         x = zout[frame, 0]
#         y = zout[frame, 1]
#         mass_patch.center = (x, y)
#         toe_patch.center = (toeout[frame, 0], toeout[frame, 1])
#         spring_line.set_data([x, toeout[frame, 0]], [y, toeout[frame, 1]])
#         ax.set_xlim(x - 1, x + 1)
#         return mass_patch, toe_patch, spring_line

#     ani = animation.FuncAnimation(fig, update, frames=len(tout),
#                                   init_func=init, blit=True, interval=dt*10)
    
#     plt.close(fig)
    
#     return ani

# ani = make_animation(tout, zout, toeout)

# # ローカル環境でアニメーションを保存する場合
# ani.save("SLIP.gif", writer='Pillow', fps=30)

# # # google colabでアニメーションを表示する場合
# # from IPython.display import HTML
# # HTML(ani.to_jshtml())