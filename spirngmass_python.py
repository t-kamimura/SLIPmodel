# -*- coding: utf-8 -*-

import numpy as np # 数値計算ライブラリの読み込み
import matplotlib.pyplot as plt # 描画ライブラリの読み込み


# パラメータの定義
m = 1.0
k = 100.0
c = 1.0

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

# 以下，メイン処理
# 初期値、刻み幅、シミュレーション時間の定義
y0 = np.array([1.0, 0.0])
dt = 5*1e-2
tend = 3.0

# 数値積分の実行
# tout, yout = Eular(y0,dt,tend)
tout, yout = RK4(y0, dt, tend)
x_exact = exact(y0, tout)

# 結果の描画
plt.figure()                  # 新しい図を作成
plt.plot(tout,yout)           # 数値解の描画
plt.plot(tout,x_exact,"k--")  # 厳密解の描画
plt.xlabel("t")               # x軸ラベルの設定
plt.ylabel("x,v")             # y軸ラベルの設定
plt.legend(("x","v"), loc="best") # 凡例の表示
plt.show()                    # 描画の実行

# アニメーションの作成
import matplotlib.animation as animation
def make_animation(tout, yout):
    fig, ax = plt.subplots()
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_aspect('equal')
    ax.set_title("Spring-Mass-Damper System")
    ax.axis('off')

    # バネと質点の描画要素
    spring_line, = ax.plot([], [], 'k-', lw=2)
    mass_patch = plt.Circle((0, 0), 0.1, fc='k')
    ax.add_patch(mass_patch)

    def init():
        spring_line.set_data([], [])
        mass_patch.center = (yout[0, 0], 0)
        return spring_line, mass_patch

    def update(frame):
        x = yout[frame, 0]
        spring_line.set_data([0, x], [0, 0])
        mass_patch.center = (x, 0)
        return spring_line, mass_patch

    ani = animation.FuncAnimation(fig, update, frames=len(tout),
                                  init_func=init, blit=True, interval=dt*1000)
    
    plt.close(fig)
    
    return ani

ani = make_animation(tout, yout)

# ローカル環境でアニメーションを保存する場合
ani.save("spring_mass_damper.gif", writer='Pillow', fps=30)

# # google colabでアニメーションを表示する場合
# from IPython.display import HTML
# HTML(ani.to_jshtml())