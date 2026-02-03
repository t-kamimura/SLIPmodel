# -*- coding: utf-8 -*-

using LinearAlgebra
using Plots
using Printf

# グローバルパラメータ
const m = 1.0
const k = 200
const l0 = 1.0
const g = 9.81

"""
運動方程式
"""
function eom(z, phase, x0)
    x = z[1]
    y = z[2]
    u = z[3]
    v = z[4]
    dz = [u, v, 0.0, 0.0]
    
    if phase == 0
        # flight phase
        dz[3] = 0.0
        dz[4] = -g
    elseif phase == 1
        # stance phase
        l = sqrt((x - x0)^2 + y^2)
        dz[3] = -k/m * (l - l0) * (x - x0) / l
        dz[4] = -k/m * (l - l0) * y / l - g
    end
    
    return dz
end

"""
イベント判定（接地・離地・転倒）
"""
function eventHandler(z, phase, x0)
    if phase == 0
        # flight phase
        toe = z[2] - l0 * cos(gamma)  # 脚先高さ
        if toe < 0 && z[4] < 0  # 接地判定
            println("touch down")
            x0 = z[1] + z[2] * tan(gamma)  # 接地点の更新
            phase = 1
        end
    elseif phase == 1
        # stance phase
        l = sqrt((z[1] - x0)^2 + z[2]^2)
        extension = l - l0  # 脚ばねの伸び．負の値なら圧縮
        if extension > 0 && z[4] > 0  # 離地判定
            println("lift off")
            phase = 0
        end
    end
    
    if z[2] < 0.0  # 転倒判定
        println("fall down")
        phase = -1
    end
    
    return phase, x0
end

"""
エネルギー計算
"""
function calc_energy(z, phase, x0)
    x = z[1]
    y = z[2]
    u = z[3]
    v = z[4]
    
    KE = 0.5 * m * (u^2 + v^2)
    PE = m * g * y
    E = KE + PE
    
    if phase == 1
        l = sqrt((x - x0)^2 + y^2)
        SE = 0.5 * k * (l - l0)^2
        E += SE
    end
    
    return E
end

"""
4次Runge-Kuttaの実装によるシミュレーション
"""
function RK4(z0, dt, tend)
    t = 0.0
    z = z0
    x0 = 0.0  # flightからはじめる
    
    tout = [copy(t)]
    zout = [copy(z)]
    phase = 0  # 0: flight, 1: stance
    pout = [phase]
    Eout = [calc_energy(z, phase, x0)]
    toeout = [[z[1] + l0*sin(gamma), z[2] - l0*cos(gamma)]]  # 脚先位置
    
    while t < tend
        t += dt
        k1 = eom(z, phase, x0)
        k2 = eom(z .+ 0.5*k1*dt, phase, x0)
        k3 = eom(z .+ 0.5*k2*dt, phase, x0)
        k4 = eom(z .+ k3*dt, phase, x0)
        z .+= (k1 .+ 2*k2 .+ 2*k3 .+ k4) .* dt / 6
        
        # イベント判定
        phase, x0 = eventHandler(z, phase, x0)
        if phase == -1
            break  # 転倒したら終了
        end
        
        # 結果を積み上げる
        push!(tout, copy(t))
        push!(zout, copy(z))
        push!(pout, phase)
        
        if phase == 0
            toe = [z[1] + l0*sin(gamma), z[2] - l0*cos(gamma)]  # 脚先位置
        else  # Juliaではデフォルトで値型
            toe = [x0, 0.0]
        end
        push!(toeout, toe)
        push!(Eout, calc_energy(z, phase, x0))  # エネルギーを積み上げる
    end
    
    return tout, zout, pout, toeout, Eout
end

"""
アニメーション作成
"""
function make_animation(tout, zout, toeout)
    fig = @layout grid(1, 1)
    
    anim = @animate for i in 1:length(tout)
        plot(fig, xlim=(0, 10), ylim=(-0.1, 1.5), aspect_ratio=:equal,
             title="SLIP model", legend=false)
        
        # 質点の描画
        scatter!([zout[i][1]], [zout[i][2]], color=:black, markersize=10, label="mass")
        
        # 脚先の描画
        scatter!([toeout[i][1]], [toeout[i][2]], color=:black, markersize=3, label="toe")
        
        # ばねの描画
        plot!([zout[i][1], toeout[i][1]], [zout[i][2], toeout[i][2]], color=:black, lw=2, label="spring")
    end
    
    return anim
end

# 初期条件設定
# 数値例1
# y0 = 1.05
# v0 = 6.455
# q0 = [0.0, y0, v0, 0.0]  # [x, y, u, v]
# const gamma = 0.54  # 脚の角度

# # 数値例2
y0 = 1.05
v0 = 4.4023
q0 = [0.0, y0, v0, 0.0]  # [x, y, u, v]
const gamma = 0.43017  # 脚の角度

dt = 1e-3
tend = 10.0

# シミュレーション実行
println("Simulation running...")
tout, zout, pout, toeout, Eout = RK4(q0, dt, tend)

@printf("Energy error: %.6e\n", Eout[1] - Eout[end])

# 描画
# zout をマトリックスに変換
z_matrix = hcat(zout...)'

# 状態量のプロット
p1 = plot(tout, z_matrix, label=["x" "y" "dx" "dy"], 
          xlabel="time [s]", title="State variables")
display(p1)

# 軌跡のプロット
p2 = plot(z_matrix[:, 1], z_matrix[:, 2], 
          xlabel="x [m]", ylabel="y [m]", title="Trajectory", legend=false)
display(p2)

# エネルギーのプロット
p3 = plot(tout, Eout, xlabel="time [s]", ylabel="Energy [J]", 
          title="Energy over time", legend=false)
display(p3)

# フェーズのプロット
p4 = plot(tout, pout, xlabel="time [s]", ylabel="phase", 
          yticks=([0, 1], ["flight", "stance"]), title="Phase", legend=false)
display(p4)

# # アニメーション作成（オプション）
# println("Creating animation...")
# anim = make_animation(tout, zout, toeout)
# gif(anim, "slip.gif", fps=30)
