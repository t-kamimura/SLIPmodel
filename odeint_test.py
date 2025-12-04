# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import RK45
import matplotlib.pyplot as plt


#1階常微分方程式（1次反応、放射性崩壊、ロジスティック方程式）
def func_dydt(y, t):
    # リッカチ方程式
    dydt = -2*y - np.tanh(t)*y**2+np.tanh(t)

    return dydt


#2d可視化
def plot2d(t_list, y_list, t_label, y_label):
    plt.xlabel(t_label)  #x軸の名前
    plt.ylabel(y_label)  #y軸の名前
    plt.grid()  #点線の目盛りを表示
    plt.plot(t_list, y_list)

    plt.show()


#メイン実行部
if (__name__ == '__main__'):
    #常微分方程式(リッカチ方程式)
    y_init = 0.0  #初期値
    t_init = -10.0  #初期時刻
    t_bound = (-10.0, 10.0)
    y_list = RK45(func_dydt,t_init,y_init,t_bound)
    print(y_list)

    #可視化
    plot2d(t_list, y_list[:, 0], "$t$", "$y(t)$")
