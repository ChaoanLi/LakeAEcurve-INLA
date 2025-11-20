项目：基于栅格数据的湖泊 Area–Elevation (A–E) 曲线估计

## ⚠️ 注意：模型简化说明（2025更新）

**当前实现采用简化策略**：
- **主模型**：仅使用水出现频率数据（Binomial likelihood）拟合空间场
- **空间场解释**：表示相对深度（非绝对高程）
- **校准步骤**：模型拟合后，使用岸边 DEM 观测校准空间场，转换为绝对高程
- **优点**：避免多 likelihood 模型的复杂性，更稳定收敛
- **输出**：最终得到的 A-E 曲线仍基于绝对高程

详见下文第3.3节（水出现频率模型）和新增的第3.4节（校准步骤）。

---

1. 研究目标与总体思路

目标：对德州的两个湖泊，利用 3 类栅格数据，重建湖底数字高程模型（bathymetry / underwater DEM），进而计算湖泊的 Area–Elevation (A–E) 曲线，即水位高程与被水覆盖面积之间的函数关系。

使用方法（简化版）：采用课程中讲过的 潜在高斯场 + 单源观测 的层次空间模型框架，结合

高斯随机场 / SPDE–GMRF（Lec3 / Lec5 / INLA materials）

非高斯观测（Binomial，Lec4 / Lec5）

空间数据失配与多源数据融合（Lec10 Spatial Misalignment）

使用 R + INLA 完成实际拟合

实现工具：R 中使用以下包（可按需增减）

terra 或 raster：读写栅格数据

sf：处理矢量边界、多边形

INLA：潜在高斯模型 + SPDE

sp 或 sf + inlabru（可选）：帮助构造 mesh 和投影

ggplot2 / tmap：可视化

2. 数据结构与预处理假设
2.1 三类栅格数据

假设你已经有以下 3 类栅格（全部在同一坐标系 CRS 下，分辨率和网格对齐）：

岸边高程栅格：DEM

变量名示例：dem_rast

类型：浮点栅格

单位：高程（m），绝对高程（例如相对于海平面）

特征：

陆地区域（包括湖岸及周边）高程有效

水域内部可能为 NA 或错误值（可视为不可靠）

水出现频率栅格：Water Occurrence / Frequency

变量名示例：water_freq_rast

类型：浮点栅格，值域在 [0, 1] 或 [0, 100]

每个像元的含义：在观测期内被水覆盖的频率（概率）

例如：0.0 = 从未有水，1.0 = 每次观测都有水

若原始数据是百分比 [0, 100]，需转换为 [0, 1]

永久水域栅格：Permanent Water Mask

变量名示例：perm_water_rast

类型：逻辑栅格或 0/1 栅格

含义：

1：在时间序列中始终有水（permanent water）

0：非永久水域

该区域可视为湖中最深区域的一部分

2.2 衍生变量与观测计数

为了将水出现“频率”转成 Binomial 观测：

对每个像元 s，构造

N_s：该像元在观测期内被观测的次数（例如 100，统一设定即可）

y_water(s)：被水覆盖的次数，令
y_water(s) = round(p(s) * N_s)，其中 p(s) 为水出现频率（0–1）

对于永久水域像元：

设 p(s) = 1

则 y_water(s) = N_s

如果需要，也可以对从未有水的区域（p = 0）设置

p(s) = 0, y_water(s) = 0，用来约束最高高程附近区域

2.3 湖区边界与分析 mask

假设有湖区或研究区域的边界多边形：

变量名示例：lake_boundary_sf（sf 对象）

用于：

裁剪栅格，定义分析区域（mask）

构造 SPDE mesh 边界

3. 潜在空间过程与观测模型（统计建模规范）

我们在湖区内定义一个潜在的连续高程场：

对每个位置 s（在 2D 平面坐标中）：
记真实高程为 
𝑍
(
𝑠
)
Z(s)

3.1 潜在场先验：高斯随机场 / SPDE–GMRF

假设：

𝑍
(
𝑠
)
Z(s) 是一个高斯随机场（Gaussian Random Field, GRF）

协方差使用 Matérn 结构（INLA 标准 SPDE 模型），具有：

方差参数 
𝜎
2
σ
2

range 参数 
𝜌
ρ（控制空间相关尺度）

实际实现：

使用 INLA 的 SPDE 模型：inla.spde2.pcmatern 或类似接口

在湖区构造三角网格 mesh，将连续场离散为一个 GMRF

形式上：

𝑍
=
(
𝑍
(
𝑠
1
)
,
…
,
𝑍
(
𝑠
𝑚
)
)
⊤
∼
𝑁
(
0
,
𝑄
−
1
(
𝜃
)
)
Z=(Z(s
1
	​

),…,Z(s
m
	​

))
⊤
∼N(0,Q
−1
(θ))，
其中 
𝑄
(
𝜃
)
Q(θ) 为由 SPDE 定义的精度矩阵，
𝜃
θ 为超参数。

3.2 岸边 DEM 观测模型（Gaussian likelihood）

在岸边及陆地区域，我们有直接高程观测：

对于 DEM 中有效像元位置 
𝑠
s，观测为

𝑌
elev
(
𝑠
)
=
𝑍
(
𝑠
)
+
𝜀
(
𝑠
)
,
𝜀
(
𝑠
)
∼
𝑁
(
0
,
𝜏
elev
−
1
)
Y
elev
	​

(s)=Z(s)+ε(s),ε(s)∼N(0,τ
elev
−1
	​

)

解释：

𝑌
elev
(
𝑠
)
Y
elev
	​

(s)：DEM 提供的高程观测值

𝑍
(
𝑠
)
Z(s)：潜在真实高程

𝜏
elev
τ
elev
	​

：精度（方差的倒数），可设较大（即误差较小）

INLA 设置：

Likelihood："gaussian"

线性预测器：

𝜂
elev
(
𝑠
)
=
𝑍
(
𝑠
)
η
elev
	​

(s)=Z(s)

不必额外固定效应，除非需要趋势项

3.3 水出现频率观测模型（Binomial + logit link）—— 简化实现

**⚠️ 当前实现注意**：我们采用简化策略，仅使用水频率数据拟合模型，岸边 DEM 用于后处理校准。

水出现频率与相对深度关系的假设：

地势越低（相对越深），越容易积水 → 水出现频率 
𝑝
(
𝑠
)
p(s) 越高

地势越高（接近岸边），
𝑝
(
𝑠
)
p(s) 越低

简化模型（当前实现）：

对每个湖内像元 s，定义观测

𝑌
water
(
𝑠
)
∼
Binomial
(
𝑁
𝑠
,
𝑝
(
𝑠
)
)
Y
water
	​

(s)∼Binomial(N
s
	​

,p(s))

logit 链接函数：

logit
(
𝑝
(
𝑠
)
)
=
𝛼
+
𝑓
(
𝑠
)
logit(p(s))=α+f(s)

其中

𝛼
α：截距

𝑓
(
𝑠
)
f(s)：空间随机场（SPDE），表示相对深度

**注意**：f(s) 本身不是绝对高程，而是相对深度的指标（正值表示更深）

INLA 实现：

Likelihood："binomial"

线性预测器：

𝜂
water
(
𝑠
)
=
𝛼
+
𝑓
(
𝑠
)
η
water
	​

(s)=α+f(s)

Formula: y ~ -1 + intercept + f(field, model = spde)

3.3.1 校准步骤：从相对深度到绝对高程

模型拟合完成后，执行以下校准步骤：

1. 提取岸边 DEM 观测位置 
𝑠
shore
s
shore
	​

 的真实高程 
𝑍
obs
(
𝑠
shore
)
Z
obs
	​

(s
shore
	​

)

2. 提取这些位置的空间场后验预测 
𝑓
^
(
𝑠
shore
)
f
^
​

(s
shore
	​

)

3. 计算校准偏移量：

offset
=
median
(
𝑍
obs
(
𝑠
shore
)
−
𝑓
^
(
𝑠
shore
)
)
offset=median(Z
obs
	​

(s
shore
	​

)−
f
^
​

(s
shore
	​

))

（使用中位数进行鲁棒估计）

4. 对所有位置，计算绝对高程：

𝑍
calibrated
(
𝑠
)
=
𝑓
^
(
𝑠
)
+
offset
Z
calibrated
	​

(s)=
f
^
​

(s)+offset

**优点**：
- 避免多 likelihood 模型的复杂性
- 模型更稳定，易于收敛
- 最终仍能得到绝对高程和准确的 A-E 曲线

3.4 永久水域的处理

对于永久水域栅格中标记为 1 的像元：

设定 p(s) ≈ 1，在 Binomial 模型中：

𝑌
water
(
𝑠
)
=
𝑁
𝑠
Y
water
	​

(s)=N
s
	​


此时模型将强烈约束这些位置的 
𝑍
(
𝑠
)
Z(s) 必须足够小，使得预测的 
𝑝
(
𝑠
)
p(s) 接近 1。

数值实现注意：

为避免 logit(1) → 无穷大，可在代码中根据需要将 y = N_s 保持原状即可，由 Binomial + logit 模型处理；

或者在后验预测时，对极端概率进行截断。

3.5 先验设置

对空间场超参数（range, σ）采用 PC priors 或弱信息先验；

对 
𝛼
α 和 
𝛽
β 使用适度宽松的正态先验，例如

𝛼
∼
𝑁
(
0
,
10
2
)
α∼N(0,10
2
)

𝛽
∼
𝑁
(
0
,
10
2
)
β∼N(0,10
2
)，并在实际实现中通过适当 reparameterization 或正值约束体现“
𝛽
>
0
β>0”的先验知识（如在代码中让 
𝛽
=
exp
⁡
(
𝛾
)
β=exp(γ)，对 
𝛾
γ 设正态先验）。

4. INLA 实现结构（需要 Cursor 写出具体 R 代码）
4.1 栅格读取与对齐

R 端要求：

使用 terra 或 raster 读取三个栅格：

dem_rast

water_freq_rast

perm_water_rast

检查：

CRS 是否一致

分辨率是否一致

extent 是否一致

若不一致，统一重投影和重采样（如使用 project 或 resample）。

使用湖区边界 lake_boundary_sf 将栅格裁剪并 mask 到研究区域。

4.2 构建观测数据表

构建三个数据源对应的 data frame，并给出坐标 (x, y)，以便投影到 SPDE mesh：

高程观测表 elev_df：

选择陆地区域和岸边（即 DEM 有效且在湖区边界内）

提取所有像元的：

像元中心点坐标 (x, y)

高程值 elev

形成数据框：

elev_df:
  x: numeric
  y: numeric
  elev: numeric


水频率观测表 water_df：

对湖内（包括部分岸边附近）的 water_freq 像元：

提取中心坐标 (x, y)

频率值 p ∈ [0, 1]

设定统一 N_s，例如 100：

N_s = 100

y_water = round(p * N_s)

形成数据框：

water_df:
  x: numeric
  y: numeric
  y_water: integer
  N: integer


对永久水域像元：p = 1 → y_water = N。

对 p = 0 的岸边/高地像元：y_water = 0，可选择性包含。

4.3 构建 SPDE mesh

使用湖区边界多边形构建 mesh：

使用 inla.mesh.2d() 或其他合适函数

mesh 要求：

覆盖整个湖区范围

在湖岸与复杂结构附近适当加密

输出：

mesh 对象，例如 mesh

然后基于 mesh 定义 SPDE 模型：

使用 inla.spde2.pcmatern() 构建 SPDE 对象，例如 spde

4.4 将观测投影到 mesh 节点（A 矩阵）

使用 inla.spde.make.A() 将 (x, y) 坐标投影到 mesh 空间：

对高程观测

A_elev = inla.spde.make.A(mesh, loc = coords_elev)

对水频率观测

A_water = inla.spde.make.A(mesh, loc = coords_water)

4.5 构建 INLA stack（多 likelihood）

构建两个 stack：

elevation stack（Gaussian）：

响应：y = elev

观测矩阵：A = A_elev

线性预测器：仅包含空间场 Z（记为某个 index，如 field）

water stack（Binomial）：

响应：y = cbind(y_water, N)（INLA 对 binomial 的写法）

观测矩阵：A = A_water

线性预测器：包含截距 + 与 Z 相关的项。理想情况是设置

eta = alpha + beta * (-Z(s))

最后将两个 stack 合并，例如

stack_full = inla.stack(elev_stack, water_stack)

4.6 模型公式（R 中的形式）

需要 Cursor 生成一个可行的 formula，核心要求：

共有一个空间随机场 field（对应 Z(s)）

在 elevation 部分：

线性预测器约为

eta_elev = field(s)

在 water 部分：

线性预测器约为

eta_water = alpha + beta * (- field(s))

即在 formula 中引入截距与一个系数对应的随机场项（例如通过多索引、不同 copy of field 或适当构造 design matrix 实现）

由于 INLA 的公式语法较复杂，具体实现细节交给 Cursor 自动生成代码即可。要求：

family 向量包含两个元素：c("gaussian", "binomial")

在 control.family 中设置 link for binomial 为 "logit"

在 control.predictor 中允许预测所有 mesh 节点

4.7 拟合与诊断

使用 inla() 进行拟合：

输入：

formula

data 来自 inla.stack.data(stack_full)

family

control.xxx 等

输出：

后验参数估计（包括 alpha, beta, space hyperparameters）

后验场 Z(s) 在 mesh 节点上的分布

要求 Cursor 在代码中：

打印模型摘要（summary(res)）

提取并可视化：

alpha, beta 的后验均值与标准差

空间场超参数（range, sigma）

5. 从潜在场构建湖底 DEM 与 A–E 曲线
5.1 从 mesh 后验场插值回栅格

步骤：

获取 mesh 节点处的后验均值 Z_hat（以及标准差）

对湖区栅格（与原始栅格同分辨率）中每个像元中心坐标 (x, y)，构造 A 矩阵 A_pred = inla.spde.make.A(mesh, loc = coords_grid)

预测栅格上像元的高程估计：

Z_grid = A_pred %*% Z_hat

将 Z_grid 转成 terra / raster 栅格对象，得到湖底 DEM 估计（bathymetry）

可以同时构建标准差栅格，用于表示不确定性。

5.2 计算 Area–Elevation 曲线

设

Z_grid：湖区所有像元的高程估计

像元面积为 cell_area（可通过 CRS 单位计算）

定义一组水位高度水平：

例如，从 min(Z_grid) 到 max(Z_grid)，步长 delta_h（例如 0.1 m）

对每个水位高度 h：

找出所有满足 Z_grid <= h 的像元（被淹没区域）

计算面积

Area
(
ℎ
)
=
(
#
{
𝑍
(
𝑠
)
≤
ℎ
}
)
×
cell_area
Area(h)=(#{Z(s)≤h})×cell_area

收集所有 (h, Area(h))，形成 A–E 曲线数据框：

ae_df:
  elevation: numeric  # h
  area: numeric       # Area(h)


绘图：

使用 ggplot2 画 A–E 曲线：

x 轴：elevation

y 轴：area

可叠加（可选）

若通过后验样本重复计算 AE 曲线，可画出均值曲线 + 例如 95% 区间

6. 结果输出与作业可解释性要求

为了达到 A+ 作业标准，最终代码与结果应包含：

数据读取和质量检查：

打印 CRS、分辨率、extent

验证三类栅格对齐

网格和 SPDE 设置：

mesh 可视化（mesh + 湖区边界）

模型拟合结果：

alpha, beta 的后验估计解释：验证“高程越高 → 水频率越低”的关系

空间场 range, sigma 的量级解释

湖底 DEM 结果：

地图：彩色高程图（重建的湖底）+ 湖岸线 + 永久水域范围

标注最低和最高高程

A–E 曲线：

曲线图，解释：

高程较低时面积增加较快 → 盆地平缓

高程较高时面积变化缓慢 → 岸壁较陡

不确定性：

标出高不确定性区域（标准差栅格）

说明不确定性的来源（远离岸边且水频率信息弱的区域）