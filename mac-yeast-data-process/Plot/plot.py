"""
最简易懂的散点图脚本（面向初学者）。

假设：
- CSV 文件与本脚本在同一目录，文件名为 `test_predictions.csv`。
- CSV 的第一列是预测值（pred），第二列是真实值（true）。

脚本功能：
- 读取 CSV 的前两列（pred, true），将 true 放在横轴，pred 放在纵轴。
- 绘制散点图，坐标轴范围固定为 0 到 6。
- 将图片保存为 `scatter.png`（与脚本同目录），并在终端打印保存路径。

如何运行：
在终端中执行：
    python plot.py

如果系统缺少依赖，先执行：
    pip install pandas matplotlib

代码注释已写在每一步，方便不会编程的人也能理解。
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === 可调参数（改这里即可） ===
# 画布尺寸（单位：英寸）。导出像素 = 英寸 × DPI，例如 (10, 10) 配合 300 DPI -> 3000×3000 像素。
FIG_SIZE = (10, 10)
# 导出分辨率（每英寸像素点数）。数值越大，图片越清晰（文件也越大）。
SAVE_DPI = 500
# 标题/坐标轴标题字体大小（满足“如何调整 xy 轴标题和图表标题 fontsize”的需求）
TITLE_FONTSIZE = 16
LABEL_FONTSIZE = 14
# 坐标轴刻度字体大小（可选）
TICK_FONTSIZE = 12
# 图例和统计信息字体大小（可选）
LEGEND_FONTSIZE = 12
STATS_FONTSIZE = 12


def main():
    """主函数：读取 CSV（假设第一列 pred，第二列 true），绘制绿色散点；
    在左上角显示统计信息（Pearson, Spearman, R^2, MAE, N），在右下角显示拟合直线方程。

    代码尽量写明每一步，便于初学者阅读和理解。
    """

    # 1) 找到 CSV 文件（与脚本在同一目录）
    csv_path = Path(__file__).with_name("test_predictions.csv")

    # 友好提示：如果没有文件则退出
    if not csv_path.exists():
        print(f"错误：未找到文件 {csv_path}\n请确认 test_predictions.csv 与 plot.py 在同一目录。")
        return

    # 2) 读取 CSV（将表格装入 pandas.DataFrame）
    df = pd.read_csv(csv_path)

    # 3) 取出预测值和真实值：按位置索引，第一列 pred，第二列 true
    #    iloc[:, 0] 表示第 1 列， iloc[:, 1] 表示第 2 列
    pred = df.iloc[:, 0].astype(float)
    true = df.iloc[:, 1].astype(float)

    # 4) 计算一些常用统计量，便于展示在图上
    #    - Pearson 相关系数：测量线性相关程度
    #    - Spearman 相关系数：基于秩的相关，耐异常值
    #    - 拟合一条直线（最小二乘），并计算 R^2 和 MAE
    n = len(true)
    pearson_r = np.corrcoef(true, pred)[0, 1]

    # Spearman 用秩（rank）后计算皮尔逊相关作为近似实现（无需额外依赖）
    spearman_rho = np.corrcoef(true.rank(), pred.rank())[0, 1]

    # 用 numpy.polyfit 做一次线性拟合：得到 slope（斜率）和 intercept（截距）
    slope, intercept = np.polyfit(true, pred, deg=1)

    # 用拟合参数计算预测值并得出 R^2（决定系数）和 MAE（平均绝对误差）
    fitted = slope * true + intercept
    ss_res = np.sum((pred - fitted) ** 2)
    ss_tot = np.sum((pred - pred.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else float("nan")
    mae = np.mean(np.abs(pred - true))

    # 5) 开始绘图：把 true 放横轴，pred 放纵轴，点颜色改为绿色（易辨认）
    # 创建一个正方形画布，便于观察 0-6 的范围（尺寸由 FIG_SIZE 控制）
    plt.figure(figsize=FIG_SIZE)

    # scatter 参数说明（帮助初学者）
    # - s: 点的大小（标量或数组），值越大点越大，常见范围 1-50（点很多时建议 1-5）
    # - alpha: 透明度，0=完全透明，1=不透明；大量点时建议设置较小的值（如 0.1-0.4）防止遮盖
    # - c: 颜色，可以是常见颜色名 ('green')、十六进制颜色 ('#297270') 或数值数组
    # - edgecolors: 点边缘颜色，设置为 'none' 可以避免边缘绘制
    # - zorder: 图层顺序，值越大绘制在越上方
    plt.scatter(
        true,
        pred,
        s=1.5,                # 点的大小。数值越大点越大；这里适中即可
        alpha=0.15,          # 透明度。点多时用 0.2~0.4 能减少遮挡感
        c="#1b7c3d", 
        edgecolors="none", # 关闭边缘描边，避免点过小出现黑边
        zorder=1,
    )

    # 6) 画参考线 y = x（红色虚线），方便看出预测是否偏离真实值
    plt.plot([0, 6], [0, 6], "r--", linewidth=1, label="y=x")

    # 7) 画拟合直线（实线，深绿色），显示拟合的趋势
    x_line = np.linspace(0, 6, 100)
    y_line = slope * x_line + intercept
    # 为拟合线添加图例标签，直接显示成 y=ax+b 的形式，便于与 y=x 一起放在右下角。
    fit_label = f"fit: y={slope:.3f}x+{intercept:.3f}"
    plt.plot(x_line, y_line, color="green", linewidth=1.2, zorder=3, label=fit_label)

    # 添加一个图标题（用户要求），并可通过 TITLE_FONTSIZE 控制字体大小
    plt.title("Test Set Scatter: Predictions vs True", fontsize=TITLE_FONTSIZE)

    # 8) 添加轴标签与范围
    # 通过 LABEL_FONTSIZE 控制坐标轴标题字体大小
    plt.xlabel("True expression", fontsize=LABEL_FONTSIZE)
    plt.ylabel("Predicted expression", fontsize=LABEL_FONTSIZE)
    # 可选：统一调大/调小刻度文字
    plt.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    plt.xlim(0, 6)
    plt.ylim(0, 6)

    # 9) 左上角文本框：放统计信息，使用轴坐标（0.02, 0.98）并设置背景框样式
    stats_text = (
        f"N = {n}\n"
        f"Pearson r = {pearson_r:.4f}\n"
        f"Spearman ρ = {spearman_rho:.4f}\n"
        f"R² = {r2:.4f}\n"
        f"MAE = {mae:.4f}"
    )
    plt.gca().text(
        0.02,
        0.98,
        stats_text,
        transform=plt.gca().transAxes,
        fontsize=STATS_FONTSIZE,
        va="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#f0e6c8", alpha=0.9, edgecolor="#a98f5a"),
    )

    # 10) 图例：将 y=x 与拟合线的说明放到右下角（lower right），放在一起更好对比。
    #     如果你希望图例半透明，可调 framealpha；想固定白底可指定 facecolor。
    plt.legend(loc="lower right", fontsize=LEGEND_FONTSIZE, framealpha=0.95, edgecolor="#444444")

    # 12) 调整布局并保存图片
    plt.tight_layout()
    out_path = Path(__file__).with_name("scatter.png")
    # 保存时使用 SAVE_DPI；像素大小 = FIG_SIZE(英寸) × SAVE_DPI
    plt.savefig(out_path, dpi=SAVE_DPI)

    # 13) 给用户提示保存位置
    print(f"图片已保存到: {out_path}")


if __name__ == "__main__":
    main()