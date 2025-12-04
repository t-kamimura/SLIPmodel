# -*- coding: utf-8 -*-

import numpy as np # 数値計算ライブラリの読み込み
import matplotlib.pyplot as plt # 描画ライブラリの読み込み


# パラメータの定義
m = 1.0
k = 100.0
c = 9.8

def eom(z, phase, x0):
    # 運動方程式
    u = z[2]
    v = z[3]
    dz = np.array([u, v, 0.0, 0.0])
    dz[2] = 0.0
    dz[3] = -g

    return dz

def RK4(z0, dt, tend):
    t = 0.0
    z = np.copy(z0) # z=z0だと、同一アドレスが共有されてしまう
    x0 = 0.0 # flightからはじめる

    tout = np.copy(t)
    zout = np.copy(z0)
    phase  = 0 # 0: flight, 1: stance
    pout = np.copy(phase)
    
    while t < tend:
        t += dt
        k1 = eom(z, phase, x0)
        k2 = eom(z + 0.5*k1*dt, phase, x0)
        k3 = eom(z + 0.5*k2*dt, phase, x0)
        k4 = eom(z + k3*dt, phase, x0)
        z += (k1 + 2*k2 + 2*k3 + k4)*dt/6

        tout = np.vstack((tout, t)) # 結果を積み上げる
        zout = np.vstack((zout, z)) # 結果を積み上げる

    return tout, zout

    
def exact(y0,tout):
    # 厳密解
    x0 = np.copy(y0[0])
    
    xout = np.zeros(len(tout))
    for i in range(len(tout)):
        t = tout[i].item() # tout[i]はnumpyのスカラー型なので、.item()でPythonのスカラー型に変換
        xout[i] = x0
    return xout

# 以下，メイン処理
# 初期値、刻み幅、シミュレーション時間の定義
y0 = np.array([1.0, 0.1])

eout = np.zeros([5,2])
ref = np.zeros([5,5])
dtset = np.array([0.0005, 0.001, 0.002, 0.005, 0.01])
for i in range(5):
    dt = dtset[i]
    tout1, yout1 = Eular(y0, dt, 2*dt)
    error_eular = np.abs(yout1[1,0] - exact(y0, np.array([dt])))
    # print(error_eular)
    tout4, yout4 = RK4(y0, dt, 2*dt)
    exact4 = exact(y0, np.array([dt]))
    error_rk4 = yout4[1,0] - exact4.item()
    # print(error_rk4)
    eout[i,0] = error_eular
    eout[i,1] = error_rk4
    ref[i,0] = dt**1
    ref[i,1] = dt**2
    ref[i,2] = dt**3
    ref[i,3] = dt**4
    ref[i,4] = dt**5

# 結果の描画
plt.figure()                  # 新しい図を作成
plt.plot(dtset,eout[:,0], '+', label = "Eular")           # 数値解の描画
plt.plot(dtset,eout[:,1], '+', label = "RK4")           # 数値解の描画
plt.plot(dtset,ref[:,0],'--', label = 'dt^1')           # 誤差の描画
plt.plot(dtset,ref[:,1],'--', label = 'dt^2')           # 誤差の描画
plt.plot(dtset,ref[:,2],'--', label = 'dt^3')           # 誤差の描画
plt.plot(dtset,ref[:,3],'--', label = 'dt^4')           # 誤差の描画
plt.plot(dtset,ref[:,4],'--', label = 'dt^5')           # 誤差の描画
plt.xscale('log')
plt.yscale('log')
plt.legend(loc="best") # 凡例の表示
plt.show()                    # 描画の実行

# # アニメーションの作成
# import matplotlib.animation as animation
# def make_animation(tout, yout):
#     fig, ax = plt.subplots()
#     ax.set_xlim(-1.5, 1.5)
#     ax.set_ylim(-0.5, 0.5)
#     ax.set_aspect('equal')
#     ax.set_title("Spring-Mass-Damper System")
#     ax.axis('off')

#     # バネと質点の描画要素
#     spring_line, = ax.plot([], [], 'k-', lw=2)
#     mass_patch = plt.Circle((0, 0), 0.1, fc='k')
#     ax.add_patch(mass_patch)

#     def init():
#         spring_line.set_data([], [])
#         mass_patch.center = (yout[0, 0], 0)
#         return spring_line, mass_patch

#     def update(frame):
#         x = yout[frame, 0]
#         spring_line.set_data([0, x], [0, 0])
#         mass_patch.center = (x, 0)
#         return spring_line, mass_patch

#     ani = animation.FuncAnimation(fig, update, frames=len(tout),
#                                   init_func=init, blit=True, interval=dt*1000)
    
#     plt.close(fig)
    
#     return ani

# ani = make_animation(tout, yout)

# # ローカル環境でアニメーションを保存する場合
# ani.save("spring_mass_damper.gif", writer='Pillow', fps=30)

# # google colabでアニメーションを表示する場合
# from IPython.display import HTML
# HTML(ani.to_jshtml())