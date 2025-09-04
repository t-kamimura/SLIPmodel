import DifferentialEquations as DE
using Plots

mutable struct ModelParam
    const m::Float64    # モデルの質量
    const g::Float64    # 重力加速度
    const k::Float64    # バネ定数
    const l0::Float64   # バネの自然長
    x0::Float64         # 脚の接地位置
    γ::Float64          # 接地角度
    phaseflag::Bool     # 接地状態を表すフラグ（0:flight, 1:stance）
    function ModelParam(γ)
        m = 80.0
        k = 20.0*1e3
        l0 = 1.0
        g = 9.81
        x0 = 0.0
        phaseflag = 0
        new(m, g, k, l0, x0, γ, phaseflag)
    end
end

function eom!(du, u, p, t)
    # システムの運動方程式
    # u = [x, y, ẋ, ẏ̇]
    if p.phaseflag == 0
        # flight phase
        du[1] = u[3]
        du[2] = u[4]
        du[3] = 0.0
        du[4] = -p.g
    else
        # stance phase
        x = u[1]
        y = u[2]
        l = sqrt((x-p.x0)^2 + y^2)
        du[1] = u[3]
        du[2] = u[4]
        du[3] =      - (p.k/p.m)*(l-p.l0)*(x - p.x0)/l
        du[4] = -p.g - (p.k/p.m)*(l-p.l0)*y/l
    end
end

function condition!(out,u, i, integrator)
    # イベント条件
    if integrator.p.phaseflag == 0
        # flight phase -> stance phase
        y = u[2]
        l0 = integrator.p.l0
        γ = integrator.p.γ 
        out[1] = y - l0*cos(γ)
        out[2] = y
    else
        # stance phase -> flight phase
        x = u[1]
        y = u[2]
        x0 = integrator.p.x0
        l = sqrt((x-x0)^2 + y^2)
        l0 = integrator.p.l0
        out[1] = l - l0
        out[2] = y
    end
end

function affect!(integrator,idx)
    if idx == 1
        if integrator.p.phaseflag == 0
            # flight phase -> stance phase
            integrator.p.phaseflag = 1
            x = integrator.u[1]
            l0 = integrator.p.l0
            γ = integrator.p.γ
            integrator.p.x0 = x + l0 * sin(γ)
            println("touchdown at x = ", integrator.p.x0)
        else
            # stance phase -> flight phase
            integrator.p.phaseflag = 0
            println("liftoff")
        end
    else
        # fall down
        DE.terminate!(integrator) #この時点で数値積分終了
        println("fall down")
    end
end

function myplot(sol)
    # 状態ベクトルの次元を取得
    n = length(sol.u[1])  # 例：位置と速度なら n = 2

    # 各成分の時間変化を抽出
    ts = sol.t
    us = [u[i] for u in sol.u, i in 1:n]  # us[i] が i番目の成分の時系列

    # subplot風に並べる
    plot_layout = @layout [a{0.5h}; b{0.5h}]  # 2行レイアウト（高さ比指定）

    p = plot(layout=plot_layout, size=(600, 400))  # 全体サイズ指定

    for i in 1:n
        plot!(p, ts, us[i], label="u[$i]", xlabel="t", ylabel="u[$i]", subplot=i)
    end

    display(p)

end

function main()
    # 初期条件の設定
    y0 = 1.0
    v0 = 7.0
    u0 = [0.0, y0, v0, 0.0]
    γ = pi/6
    p = ModelParam(γ)
    tspan = (0.0, 1.0)

    # 常微分方程式の定義
    prob = DE.ODEProblem(eom!, u0, tspan, p)
    cb = DE.VectorContinuousCallback(condition!, affect!,2, save_positions=(true,true)) # イベントコールバックの定義
    # 解の計算
    sol = DE.solve(prob, DE.Tsit5(), callback=cb)
    println("最終時刻: ", sol.t[end])
    println("最終位置: ", sol.u[end][2])

    # 結果のプロット
    # myplot(sol)
    # plot([u[1] for u in sol.u], [u[2] for u in sol.u])
    ts = range(sol.t[1], sol.t[end], length=100)
    xs = [sol(t)[1] for t in ts]  # 補間された位置
    ys = [sol(t)[2] for t in ts]  # 補間された位置

    plot(xs, ys, label="補間された位置", xlabel="x", ylabel="y")


end
main()