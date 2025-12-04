# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt # 描画ライブラリの読み込み

t = np.linspace(0, 2*np.pi, 100) # 0から2πまでの100点
y1 = np.sin(t) # 正弦波
y2 = np.cos(t) # 余弦波

plt.figure()            # 新しい図を作成
plt.plot(t, y1)         # tを横軸、y1を縦軸にプロット
plt.plot(t, y2)         # tを横軸、y2を縦軸にプロット
plt.legend(['sin', 'cos']) # 凡例の表示
plt.xlabel('Time [s]')  # x軸のラベル
plt.ylabel('Amp')       # y軸のラベル
plt.title('Sine and Cosine Waves')  # 図のタイトル
plt.show()              # 図を表示