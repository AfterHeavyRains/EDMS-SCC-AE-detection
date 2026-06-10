import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import matplotlib
from matplotlib.ticker import AutoMinorLocator

matplotlib.use('Agg')
import matplotlib as mpl
# === ✅ 新版：四周边框 + 向内主次刻度 ===
from matplotlib.ticker import AutoMinorLocator

def add_ieee_border(ax):
    """四周边框可见，主刻度和次刻度均向内，仅左下显示刻度"""
    # 四周边框
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.9)

    # 主次刻度线向内（仅左下）
    ax.tick_params(axis='both', which='major', direction='in', length=4, width=0.8,
                   top=False, right=False, bottom=True, left=True)
    ax.tick_params(axis='both', which='minor', direction='in', length=2, width=0.6,
                   top=False, right=False, bottom=True, left=True)

    # 添加次刻度
    # ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    # 保证上下边框都显示，但刻度只在左下
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)


# === 路径 ===
repo_root = Path(__file__).resolve().parents[3]
root = repo_root / "results" / "eval_exports"
save_dir = root / "fig_ieee_full"
save_dir.mkdir(exist_ok=True)

# files = {
#     '2mm': root / 'per_sample_2mm.csv',
#     '4mm': root / 'per_sample_4mm.csv',
#     '6mm': root / 'per_sample_6mm.csv',
#     '8mm': root / 'per_sample_8mm.csv'
# }

files = {
    '2mm': root / 'per_sample_2mm_full_eval_withCovOL.csv',
    '4mm': root / 'per_sample_4mm_full_eval_withCovOL.csv',
    '6mm': root / 'per_sample_6mm_full_eval_withCovOL.csv',
    '8mm': root / 'per_sample_8mm_full_eval_withCovOL.csv'
}

# === 数据 ===
dfs = []
for thick, file in files.items():
    df = pd.read_csv(file)
    df['Thickness'] = thick
    dfs.append(df)
data = pd.concat(dfs, ignore_index=True)

cols = ['AGRo','AGRl','EORo','EORl','CovO','CovL']
summary = data.groupby('Thickness')[cols].mean().round(3)
print(summary)

# === 全局绘图风格 ===
palette_ieee = ["#1F77B4", "#00A6D6", "#7FC8F8", "#A7B6C2"]
sns.set_theme(context="paper", style="white", font="Times New Roman", font_scale=1.2)
# === IEEE / Elsevier 图表全局格式 ===
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.titlesize": 13,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "axes.linewidth": 0.9,
    "lines.linewidth": 0.9,
    "lines.markersize": 5,
    "grid.linewidth": 0.6,
    "grid.alpha": 0.5,
    "axes.edgecolor": "black",
    "axes.labelpad": 4,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "figure.dpi": 650,
    "savefig.dpi": 650,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

sns.set_palette(palette_ieee)

# # === 函数：统一IEEE边框+主次刻度 ===
# def add_ieee_border(ax):
#     for spine in ax.spines.values():
#         spine.set_visible(True)
#         spine.set_color("black")
#         spine.set_linewidth(0.9)
#     ax.tick_params(axis='both', which='major', direction='out', length=4, width=0.8)
#     ax.tick_params(axis='both', which='minor', direction='out', length=2, width=0.6)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.spines['top'].set_visible(True)
#     ax.spines['right'].set_visible(True)

# === 标签映射 ===
label_dict = {
    "AGRl": "Detection Accuracy (AGR)",
    "EORl": "False Alarm Rate (EOR)",
    "CovL": "Detection Coverage (Cov)"
}

# =============================================================
# 1️⃣ AGR 条形图
# =============================================================
fig, ax = plt.subplots(figsize=(7, 5))
sns.barplot(
    data=data.melt(id_vars='Thickness', value_vars=['AGRo','AGRl']),
    x='Thickness', y='value', hue='variable',
    palette=["#8FB9A8", "#F4B183"], width=0.6, saturation=0.9, ax=ax
)
add_ieee_border(ax)
ax.set_xlabel("Specimen Thickness", fontsize=10)
ax.set_ylabel(label_dict["AGRl"], fontsize=10)
ax.legend(title="", loc='upper right', frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig(save_dir / "Fig1_AGR_bar_ieee_full2.png", dpi=650, bbox_inches='tight')
plt.close()

# =============================================================
# 2️⃣ EOR 条形图
# =============================================================
fig, ax = plt.subplots(figsize=(7, 5))
sns.barplot(
    data=data.melt(id_vars='Thickness', value_vars=['EORo','EORl']),
    x='Thickness', y='value', hue='variable',
    palette=["#1f77b4", "#00a0b0"], width=0.6, saturation=0.9, ax=ax
)
add_ieee_border(ax)
ax.set_xlabel("Specimen Thickness", fontsize=10)
ax.set_ylabel(label_dict["EORl"], fontsize=10)
ax.legend(title="", loc='upper right', frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig(save_dir / "Fig2_EOR_bar_ieee_full2.png", dpi=650, bbox_inches='tight')
plt.close()

# =============================================================
# # 3️⃣ Pareto Cloud
# # =============================================================
# fig, ax = plt.subplots(figsize=(6.5, 5.5))
# sns.scatterplot(
#     data=data, x='EORl', y='AGRl', hue='Thickness',
#     palette=palette_ieee, s=55, alpha=0.9,
#     edgecolor='black', linewidth=0.3, ax=ax
# )
# add_ieee_border(ax)
# ax.set_xlabel(label_dict["EORl"], fontsize=12, style='italic')
# ax.set_ylabel(label_dict["AGRl"], fontsize=12, style='italic')
# ax.legend(title="Thickness", loc='lower left', frameon=False)
# plt.tight_layout()
# plt.savefig(save_dir / "Fig3_ParetoCloud_ieee_full.png", dpi=600, bbox_inches='tight')
# plt.close()
#
# # =============================================================
# # 4️⃣ Coverage vs AGR
# # =============================================================
# fig, ax = plt.subplots(figsize=(6.5, 5.5))
# sns.scatterplot(
#     data=data, x='CovL', y='AGRl', hue='Thickness', style='Thickness',
#     palette=palette_ieee, s=70, alpha=0.9,
#     edgecolor='black', linewidth=0.4, ax=ax
# )
# add_ieee_border(ax)
# ax.set_xlabel(label_dict["CovL"], fontsize=12, style='italic')
# ax.set_ylabel(label_dict["AGRl"], fontsize=12, style='italic')
# ax.legend(title="Thickness", loc='lower right', frameon=False)
# plt.tight_layout()
# plt.savefig(save_dir / "Fig4_Cov_vs_AGR_ieee_full.png", dpi=600, bbox_inches='tight')
# plt.close()
fig, ax = plt.subplots(figsize=(6.2, 5.2))
palette_ieee2 = {
    '2mm': '#103C73',   # 深蓝 — 代表薄板
    '4mm': '#1D6FA3',   # 蓝青 — 平衡点
    '6mm': '#3FA7D6',   # 湖蓝 — 视觉重心
    '8mm': '#9DC3C2'    # 灰青 — 厚板、低能级
}

# 背景密度可视化（增强分布形态）
sns.kdeplot(
    data=data, x='EORl', y='AGRl', hue='Thickness',
    fill=True, alpha=0.18, levels=8, thresh=0.06,
    palette=palette_ieee2, linewidths=0.4, ax=ax
)


# 前景散点
sns.scatterplot(
    data=data, x='EORl', y='AGRl', hue='Thickness',
    palette=palette_ieee2, s=60, alpha=0.85,
    edgecolor='black', linewidth=0.3, ax=ax, legend=False
)


# 坐标范围 & 样式
ax.set_xlim(-0.02, 0.75)
ax.set_ylim(0.74, 1.01)
ax.set_xlabel(r'False Alarm Rate (EOR)', fontsize=15, fontname='Times New Roman')
ax.set_ylabel(r'Detection Accuracy (AGR)', fontsize=15, fontname='Times New Roman')
ax.grid(True, linestyle=':', linewidth=0.6, alpha=0.6)
add_ieee_border(ax)

# 图例单独置底
sns.move_legend(ax, "lower right", title="Thickness", frameon=False)
plt.tight_layout()
plt.savefig(save_dir / "Fig3_ParetoCloud_opt2.png", dpi=650, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(6.2, 5.2))

sns.kdeplot(
    data=data, x='CovL', y='AGRl', hue='Thickness',
    fill=True, alpha=0.18, levels=8, thresh=0.06,
    palette=palette_ieee2, linewidths=0.4, ax=ax
)

sns.scatterplot(
    data=data, x='CovL', y='AGRl',hue='Thickness',
    palette=palette_ieee2, s=70, alpha=0.85,
    edgecolor='black', linewidth=0.3, ax=ax, legend=False
)

ax.set_xlim(0.45, 1.01)
ax.set_ylim(0.74, 1.01)
ax.set_xlabel(r'Detection Coverage (Cov)', fontsize=15, fontname='Times New Roman')
ax.set_ylabel(r'Detection Accuracy (AGR)', fontsize=15, fontname='Times New Roman')

ax.grid(True, linestyle=':', linewidth=0.6, alpha=0.6)
add_ieee_border(ax)
sns.move_legend(ax, "lower left", title="Thickness", frameon=False)
plt.tight_layout()
plt.savefig(save_dir / "Fig4_Cov_vs_AGR_opt2.png", dpi=650, bbox_inches='tight')
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13

# =============================================================
# 5️⃣ KDE 分布图
# =============================================================
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
metrics = ['AGRl', 'EORl', 'CovL']
titles = ['(a) Accuracy', '(b) False Alarm', '(c) Coverage']
for i, ax in enumerate(axes):
    sns.kdeplot(data=data, x=metrics[i], hue='Thickness',
                fill=True, alpha=0.35, linewidth=1.0, ax=ax)
    add_ieee_border(ax)
    ax.set_title(titles[i], fontsize=13)
    ax.set_xlabel(label_dict[metrics[i]], fontsize=10, style='italic')
    ax.set_ylabel("Density", fontsize=12)
plt.tight_layout()
plt.savefig(save_dir / "Fig5_KDE_ieee_full2.png", dpi=650, bbox_inches='tight')
plt.close()

# =============================================================
# 6️⃣ 小提琴图
# =============================================================
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
metrics = ['AGRl', 'EORl', 'CovL']
for i, ax in enumerate(axes):
    sns.violinplot(data=data, x='Thickness', y=metrics[i],
                   palette=palette_ieee, inner="quartile", linewidth=0.9, cut=0, ax=ax)
    add_ieee_border(ax)
    ax.set_xlabel("Specimen Thickness", fontsize=15)
    ax.set_ylabel(label_dict[metrics[i]], fontsize=15)
plt.tight_layout()
plt.savefig(save_dir / "Fig6_Violin_ieee_full2.png", dpi=650, bbox_inches='tight')
plt.close()

print(f"\nAll 6 IEEE-style figures saved to: {save_dir}")
