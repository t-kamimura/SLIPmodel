// SLIP (Spring Loaded Inverted Pendulum) モデルシミュレーション
// Python版からの移植

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// パラメータ設定
const double m = 1.0;      // 質量
const double k = 200.0;    // バネ定数
const double l0 = 1.0;     // 自然長
const double g = 9.81;     // 重力加速度

// 状態ベクトル構造体
struct State {
    double x;  // 水平位置
    double y;  // 垂直位置
    double u;  // 水平速度
    double v;  // 垂直速度
    
    State() : x(0), y(0), u(0), v(0) {}
    State(double x_, double y_, double u_, double v_) 
        : x(x_), y(y_), u(u_), v(v_) {}
    
    State operator+(const State& s) const {
        return State(x + s.x, y + s.y, u + s.u, v + s.v);
    }
    
    State operator*(double scalar) const {
        return State(x * scalar, y * scalar, u * scalar, v * scalar);
    }
};

// 運動方程式
State eom(const State& z, int phase, double x0) {
    State dz;
    dz.x = z.u;
    dz.y = z.v;
    
    if (phase == 0) {
        // flight phase
        dz.u = 0.0;
        dz.v = -g;
    } else if (phase == 1) {
        // stance phase
        double l = std::sqrt((z.x - x0) * (z.x - x0) + z.y * z.y);
        dz.u = -k/m * (l - l0) * (z.x - x0) / l;
        dz.v = -k/m * (l - l0) * z.y / l - g;
    }
    
    return dz;
}

// イベントハンドラー
void eventHandler(const State& z, int& phase, double& x0, double gamma) {
    if (phase == 0) {
        // flight phase
        double toe = z.y - l0 * std::cos(gamma);  // 脚先高さ
        if (toe < 0 && z.v < 0) {  // 接地判定
            std::cout << "touch down" << std::endl;
            x0 = z.x + z.y * std::tan(gamma);  // 接地点の更新
            phase = 1;
        }
    } else if (phase == 1) {
        // stance phase
        double l = std::sqrt((z.x - x0) * (z.x - x0) + z.y * z.y);
        double extension = l - l0;  // 脚ばねの伸び
        if (extension > 0 && z.v > 0) {  // 離地判定
            std::cout << "lift off" << std::endl;
            phase = 0;
        }
    }
    
    if (z.y < 0.0) {  // 転倒判定
        std::cout << "fall down" << std::endl;
        phase = -1;
    }
}

// エネルギー計算
double calc_energy(const State& z, int phase, double x0) {
    double KE = 0.5 * m * (z.u * z.u + z.v * z.v);
    double PE = m * g * z.y;
    double E = KE + PE;
    
    if (phase == 1) {
        double l = std::sqrt((z.x - x0) * (z.x - x0) + z.y * z.y);
        double SE = 0.5 * k * (l - l0) * (l - l0);
        E += SE;
    }
    
    return E;
}

// RK4積分法によるシミュレーション
void RK4(const State& z0, double dt, double tend, double gamma,
         std::vector<double>& tout, std::vector<State>& zout, 
         std::vector<int>& pout, std::vector<std::pair<double, double>>& toeout,
         std::vector<double>& Eout) {
    
    double t = 0.0;
    State z = z0;
    double x0 = 0.0;  // flightから始める
    int phase = 0;    // 0: flight, 1: stance
    
    // 初期値を保存
    tout.push_back(t);
    zout.push_back(z);
    pout.push_back(phase);
    toeout.push_back(std::make_pair(z.x + l0 * std::sin(gamma), 
                                     z.y - l0 * std::cos(gamma)));
    Eout.push_back(calc_energy(z, phase, x0));
    
    while (t < tend) {
        t += dt;
        
        // RK4積分
        State k1 = eom(z, phase, x0);
        State k2 = eom(z + k1 * (0.5 * dt), phase, x0);
        State k3 = eom(z + k2 * (0.5 * dt), phase, x0);
        State k4 = eom(z + k3 * dt, phase, x0);
        z = z + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
        
        // イベント判定
        eventHandler(z, phase, x0, gamma);
        if (phase == -1) {
            break;  // 転倒したら終了
        }
        
        // 結果を保存
        tout.push_back(t);
        zout.push_back(z);
        pout.push_back(phase);
        
        if (phase == 0) {
            toeout.push_back(std::make_pair(z.x + l0 * std::sin(gamma), 
                                             z.y - l0 * std::cos(gamma)));
        } else {
            toeout.push_back(std::make_pair(x0, 0.0));
        }
        
        Eout.push_back(calc_energy(z, phase, x0));
    }
}

// 結果をCSVファイルに保存
void save_to_csv(const std::string& filename, 
                 const std::vector<double>& tout,
                 const std::vector<State>& zout,
                 const std::vector<int>& pout,
                 const std::vector<std::pair<double, double>>& toeout,
                 const std::vector<double>& Eout) {
    
    std::ofstream file(filename);
    file << std::fixed << std::setprecision(6);
    file << "time,x,y,u,v,phase,toe_x,toe_y,energy\n";
    
    for (size_t i = 0; i < tout.size(); ++i) {
        file << tout[i] << ","
             << zout[i].x << ","
             << zout[i].y << ","
             << zout[i].u << ","
             << zout[i].v << ","
             << pout[i] << ","
             << toeout[i].first << ","
             << toeout[i].second << ","
             << Eout[i] << "\n";
    }
    
    file.close();
    std::cout << "Results saved to " << filename << std::endl;
}

int main() {
    // 初期条件設定
    // 数値例1
    // double y0 = 1.05;
    // double v0 = 6.455;
    // State q0(0.0, y0, v0, 0.0);
    // double gamma = 0.54;
    
    // 数値例2
    double y0 = 1.05;
    double v0 = 4.4023;
    State q0(0.0, y0, v0, 0.0);
    double gamma = 0.43017;
    
    double dt = 1e-3;
    double tend = 10.0;
    
    // シミュレーション結果を格納する配列
    std::vector<double> tout;
    std::vector<State> zout;
    std::vector<int> pout;
    std::vector<std::pair<double, double>> toeout;
    std::vector<double> Eout;
    
    // シミュレーション実行
    std::cout << "Starting simulation..." << std::endl;
    RK4(q0, dt, tend, gamma, tout, zout, pout, toeout, Eout);
    
    // エネルギー誤差の表示
    std::cout << "Energy error: " << (Eout[0] - Eout[Eout.size() - 1]) << std::endl;
    
    // 結果をCSVファイルに保存
    save_to_csv("slip_results.csv", tout, zout, pout, toeout, Eout);
    
    std::cout << "Simulation completed. Total steps: " << tout.size() << std::endl;
    
    return 0;
}
