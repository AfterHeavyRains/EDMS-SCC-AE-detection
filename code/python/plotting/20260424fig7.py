import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# 数据处理
methods = ['Hilbert Envelope', 'Improved AIC', 'Spectral Kurtosis',
           'Kurtogram', 'Shannon Entropy', 'EDMS (Proposed)']

# 计算 4 种厚度的平均值
AGR = [np.mean([0.948, 0.921, 0.877, 0.924]), np.mean([0.673, 0.760, 0.467, 0.656]),
       np.mean([0.511, 0.645, 0.349, 0.493]), np.mean([0.545, 0.637, 0.299, 0.504]),
       np.mean([0.705, 0.732, 0.385, 0.639]), np.mean([0.956, 0.966, 0.910, 0.940])]

Cov = [np.mean([0.752, 0.721, 0.510, 0.686]), np.mean([0.649, 0.389, 0.582, 0.541]),
       np.mean([0.991, 0.884, 0.953, 0.970]), np.mean([0.994, 0.890, 0.988, 0.983]),
       np.mean([0.677, 0.487, 0.538, 0.608]), np.mean([0.929, 0.867, 0.930, 0.927])]

EOR = [np.mean([0.052, 0.079, 0.123, 0.076]), np.mean([0.327, 0.240, 0.533, 0.344]),
       np.mean([0.489, 0.355, 0.651, 0.507]), np.mean([0.455, 0.363, 0.701, 0.496]),
       np.mean([0.295, 0.268, 0.615, 0.361]), np.mean([0.071, 0.133, 0.070, 0.073])]

bubble_sizes = [(1 - e) * 800 for e in EOR]

fig, ax = plt.subplots(figsize=(6, 5.5))

# 隐藏默认框线
for side in ['top', 'right', 'left', 'bottom']:
    ax.spines[side].set_visible(False)

# 绘制数据：所有方法都用圆形
for i in range(5):
    ax.scatter(Cov[i], AGR[i], s=bubble_sizes[i], color='#4A79A5', alpha=0.7, zorder=3)

# EDMS 改回圆形，加大尺寸并使用红色，边缘加粗
ax.scatter(Cov[5], AGR[5], s=bubble_sizes[5], color='#EE3224',
           edgecolors='black', linewidths=1.2, zorder=5)

# --- 关键修改：手动绘制带箭头的坐标轴，止于 1.0 ---
ax.annotate('', xy=(1.02, 0.42), xytext=(0.45, 0.42),
            arrowprops=dict(arrowstyle='->', lw=1.5, color='black', shrinkA=0, shrinkB=0))
ax.annotate('', xy=(0.45, 1.02), xytext=(0.45, 0.42),
            arrowprops=dict(arrowstyle='->', lw=1.5, color='black', shrinkA=0, shrinkB=0))

# 设置刻度：只显示到 1.0
ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

# 标注文字位置
label_params = [
    {'text': 'Hilbert Envelope', 'pos': (0.61, 0.86), 'ha': 'left'},
    {'text': 'Improved AIC', 'pos': (0.50, 0.67), 'ha': 'left'},
    {'text': 'Spectral Kurtosis', 'pos': (0.94, 0.53), 'ha': 'left'},
    {'text': 'Kurtogram', 'pos': (0.85, 0.49), 'ha': 'left'},
    {'text': 'Shannon Entropy', 'pos': (0.54, 0.57), 'ha': 'left'},
    {'text': 'EDMS (Proposed)', 'pos': (0.97, 0.89), 'ha': 'right', 'weight': 'bold'},
]

for lp in label_params:
    ax.text(lp['pos'][0], lp['pos'][1], lp['text'], fontsize=10,
            ha=lp['ha'], fontweight=lp.get('weight', 'normal'))

# 设置物理限制，留出箭头和圆形边缘的空间
ax.set_xlim(0.45, 1.08)
ax.set_ylim(0.42, 1.08)

ax.set_xlabel('Detection Coverage (Cov)', fontsize=12, loc='right')
ax.set_ylabel('Detection Agreement (AGR)', fontsize=12, loc='center')
ax.tick_params(direction='out', labelsize=10)

plt.tight_layout()
out_dir = Path(__file__).resolve().parents[3] / "results" / "figures" / "check_run"
out_dir.mkdir(parents=True, exist_ok=True)
plt.savefig(out_dir / 'Fig9_Arrow_Circle_Style.png', dpi=600, bbox_inches='tight')
plt.show()
