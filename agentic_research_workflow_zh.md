# 智能体辅助研究工作流分析与可复现搭建指南

**cqed_sim + VS Code + 多智能体模型**

> 本报告为外部研究人员提供全面的技术文档，使其能够理解、复现并改编
> 为电路量子电动力学（cQED）仿真研究开发的自主智能体辅助研究工作流。

---

## 目录

1. [执行摘要](#1-执行摘要)
2. [第一部分 — cqed_sim：核心仿真器分析](#第一部分--cqed_sim核心仿真器分析)
   - [1.1 架构](#11-架构)
   - [1.2 物理与约定](#12-物理与约定)
   - [1.3 仿真引擎](#13-仿真引擎)
   - [1.4 正确性与验证](#14-正确性与验证)
   - [1.5 智能体集成](#15-智能体集成)
3. [第二部分 — cqed_based_study：智能体研究层](#第二部分--cqed_based_study智能体研究层)
   - [2.1 课题结构](#21-课题结构)
   - [2.2 执行流水线](#22-执行流水线)
   - [2.3 智能体角色](#23-智能体角色)
   - [2.4 输出产物](#24-输出产物)
   - [2.5 故障处理](#25-故障处理)
4. [第三部分 — 智能体工作流设计](#第三部分--智能体工作流设计)
   - [3.1 多智能体系统设计](#31-多智能体系统设计)
   - [3.2 迭代循环](#32-迭代循环)
   - [3.3 提示策略](#33-提示策略)
   - [3.4 AGENTS.md 系统](#34-agentsmd-系统)
5. [第四部分 — VS Code 搭建指南（可复现）](#第四部分--vs-code-搭建指南可复现)
   - [4.1 所需工具](#41-所需工具)
   - [4.2 模型使用策略](#42-模型使用策略)
   - [4.3 VS Code 推荐工作流](#43-vs-code-推荐工作流)
   - [4.4 自动化循环（高级）](#44-自动化循环高级)
   - [4.5 终端/执行集成](#45-终端执行集成)
6. [第五部分 — 向其他领域的推广](#第五部分--向其他领域的推广)
7. [第六部分 — 建议与改进](#第六部分--建议与改进)
8. [附录 A — 文件树参考](#附录-a--文件树参考)
9. [附录 B — 完整智能体定义文件](#附录-b--完整智能体定义文件)
10. [附录 C — 术语表](#附录-c--术语表)

---

## 1. 执行摘要

本仓库实现了一套**全自主、智能体辅助的研究工作流**，用于电路量子电动力学（cQED）计算物理课题。该系统整合了：

- **`cqed_sim`** — 基于 QuTiP 的硬件保真、脉冲级 cQED 仿真器
- **`cqed_based_study`** — 包含智能体定义、技能、工具和课题模板的结构化研究工作区
- **多模型 AI 智能体** — 通过 VS Code 的 GitHub Copilot 基础设施协调，使用 Claude Opus 4.6、OpenAI Codex 和 Copilot Chat

该工作流以**持续研究循环**方式运行：用户提出研究问题，系统自主规划实验、编写仿真代码、运行仿真、验证结果，并生成出版质量的 LaTeX 报告——所有过程均具备结构化状态管理，支持中断恢复和可续执行。

**六项已完成的课题**展示了该工作流的能力：
1. 色散读出脉冲优化（4 个递进子课题）
2. 灰箱自适应控制 cQED 系统
3. 热噪声腔体传感
4. SQR 门设计
5. 混合量子比特-腔体控制
6. 文献驱动的选择性脉冲基元

所有课题均已重新运行并端到端验证，在测试套件中达到 69/69 项验证检查全部通过。

---

## 第一部分 — cqed_sim：核心仿真器分析

### 1.1 架构

#### 设计理念

`cqed_sim` 是一款**硬件保真、时域、脉冲级电路量子电动力学仿真器**。它以完整的色散哈密顿量建模，包含显式驱动波形、硬件失真效应和 Lindblad 开放系统动力学。设计优先考虑：

1. **物理保真度** — 模型包含高阶色散项（下降阶乘展开）、多模耦合和硬件信号链效应
2. **可组合性** — 冻结数据类模型、缓存算符和构建器模式，允许物理组件的灵活组合
3. **可复现性** — 确定性求解器后端、单位一致的约定、JSON 可序列化的门序列
4. **可扩展性** — 插件式后端系统（QuTiP、NumPy 密集、JAX 密集）、模块化子包

#### 模块结构

框架按清晰分离的层次组织：

```
cqed_sim/
├── core/           ← 物理模型、希尔伯特空间、理想态/门
│   ├── models.py          # DispersiveTransmonCavityModel, UniversalCQEDModel
│   ├── frames.py          # FrameSpec, 旋转坐标系定义
│   ├── ideal_gates.py     # 理想酉门库
│   ├── state_prep.py      # 初态构造
│   └── energy_spectrum.py # 修饰能级求解器
│
├── pulses/         ← 脉冲构造与校准
│   ├── pulse.py           # Pulse 数据类（冻结）
│   ├── envelope.py        # GaussianEnvelope, FlatTopEnvelope 等
│   ├── builders.py        # build_displacement_pulse, build_rotation_pulse, ...
│   ├── calibration.py     # 解析脉冲校准公式
│   └── hardware.py        # HardwareConfig（IQ 失真、ZOH、DAC）
│
├── sequence/       ← 波形编译流水线
│   └── compiler.py        # SequenceCompiler：脉冲列表 → 采样时间线
│
├── sim/            ← 仿真引擎
│   ├── engine.py          # simulate_sequence(), SimulationSession
│   ├── hamiltonian.py     # 由模型 + 驱动组装哈密顿量
│   ├── noise.py           # NoiseSpec, Lindblad 衰减算符
│   ├── result.py          # SimulationResult 含态提取器
│   └── diagnostics.py     # 运行时检查与收敛监控
│
├── backends/       ← 可插拽求解器后端
│   ├── base.py            # BaseBackend 抽象基类
│   ├── numpy_backend.py   # 密集矩阵 expm 传播
│   └── jax_backend.py     # JAX 加速密集后端
│
├── measurement/    ← 读出建模
│   ├── qubit.py           # 理想、混淆矩阵和 IQ 测量
│   ├── readout_chain.py   # 谐振腔 + Purcell 滤波器 + 放大器
│   └── continuous.py      # 基于 SME 的连续监测
│
├── floquet/        ← 周期驱动 Floquet 分析
├── analysis/       ← 参数转换（裸态 → 色散态）
├── calibration/    ← SQR 门校准、多色调验证
├── calibration_targets/  ← 代理实验（光谱学、Rabi、Ramsey 等）
├── gates/          ← 理想酉门库（100+ 门）
├── operators/      ← 缓存 Pauli、升降算符、Fock 投影算符
├── observables/    ← Bloch、Fock 分辨、相位、Wigner 诊断
├── plotting/       ← 出版质量可视化
├── tomo/           ← 态/过程层析、泄漏矩阵校准
├── unitary_synthesis/  ← 门序列最优控制
├── optimal_control/    ← 含硬件感知前向模型的直接 GRAPE
├── rl_control/     ← 强化学习环境
├── system_id/      ← 校准信息先验
├── quantum_algorithms/ ← 全息算法工具
└── io/             ← 门序列 JSON 输入/输出
```

#### 关注点分离

| 层 | 职责 | 关键类 |
|---|------|--------|
| **物理模型** | 哈密顿量参数、耦合规格、旋转坐标系 | `UniversalCQEDModel`, `DispersiveTransmonCavityModel`, `FrameSpec` |
| **脉冲构造** | 波形生成、包络成形、校准 | `Pulse`, `PulseBuilder`, `HardwareConfig` |
| **数值求解器** | 时间演化、态传播、噪声通道 | `simulate_sequence()`, `SimulationSession`, `NoiseSpec` |
| **分析 / API** | 结果提取、绘图、层析、优化 | `SimulationResult`, `UnitarySynthesizer`, `GrapeSolver` |

### 1.2 物理与约定

#### 哈密顿量构造

旋转坐标系中双模（量子比特 + 腔体）系统的色散哈密顿量：

$$H_0/\hbar = \delta_c \hat{n}_c + \delta_q \hat{n}_q + \frac{\alpha}{2} \hat{b}^{\dagger 2}\hat{b}^2 + \frac{K}{2}\hat{n}_c(\hat{n}_c - 1) + \chi \hat{n}_c \hat{n}_q + \chi_2 \hat{n}_c(\hat{n}_c - 1)\hat{n}_q + \cdots$$

其中：
- $\delta_c, \delta_q$ = 相对旋转坐标系的失谐量
- $\alpha$ = transmon 非谐性（负值，通常 −255 MHz）
- $K$ = 腔体自 Kerr 效应（通常 −28 kHz）
- $\chi$ = 色散频移（参考参数集中为 −2.84 MHz）
- $\chi_2$ = 二阶色散频移

**高阶项**使用**下降阶乘**形式：$\chi_{\text{higher}}[i]$ 乘以 $\hat{n}(\hat{n}-1)\cdots(\hat{n}-i)$，而非 $\hat{n}^{i+1}$。

#### 符号约定

| 物理量 | 约定 | 物理含义 |
|--------|------|----------|
| $\chi$ | 正 $\chi$ → 量子比特频率随光子数**增加** | $\omega_q(n) = \omega_q + \chi \cdot n$ |
| $\alpha$ | transmon 为负值 | $\omega_{ef} = \omega_{ge} + \alpha$ |
| $K$ | 负值（腔体自 Kerr） | 光子数依赖的频率牵引 |
| 驱动载波 | $\omega_{\text{carrier}} = -\omega_{\text{transition}}$ | 使用 $e^{+i\omega t}$ 约定 |

#### 驱动哈密顿量

$$H_{\text{drive}} = \epsilon(t) \hat{O}^+ + \epsilon^*(t) \hat{O}^-$$

其中 $\hat{O}^{\pm}$ 是从字符串目标（`"qubit"`, `"storage"`, `"sideband"`）或结构化规格（`TransmonTransitionDriveSpec`, `SidebandDriveSpec`）解析得到的升降算符。

#### 旋转坐标系

`FrameSpec` 为每个模式定义旋转坐标系。默认坐标系以模式频率旋转，使色散频移成为主要能量尺度。所有脉冲载波均相对于该坐标系定义。

#### 与物理约定文档的一致性

软件包中包含 `physics_and_conventions/conventions.py`，提供枚举类型（`UnitType`, `DetuningSign`, `TensorOrdering`）及强制装饰器，确保所有模块间符号和顺序的一致性。物理约定文档（`physics_conventions_report.tex`）提供了从 Jaynes-Cummings 模型到色散极限的完整推导链。

### 1.3 仿真引擎

#### 求解器后端

**主路径（QuTiP）：**
- 纯态：`qutip.sesolve`（薛定谔方程）
- 混合态 / 开放系统：`qutip.mesolve`（Lindblad 主方程）
- 稳态：`qutip.steadystate`（平衡态计算）

**密集后端路径（NumPy / JAX）：**
- 分段常值矩阵指数传播（`scipy.linalg.expm` 或 `jax.scipy.linalg.expm`）
- 适用于小系统和校验对比
- JAX 后端支持 GPU 加速的 GRAPE 最优控制

#### 时间演化流水线

```
┌──────────────┐     ┌──────────────┐     ┌──────────────────┐     ┌───────────────┐
│  物理模型     │     │  脉冲对象     │     │  序列编译器       │     │  仿真引擎     │
│              │────▶│              │────▶│                  │────▶│               │
│ • 参数       │     │ • 包络       │     │ • 采样脉冲       │     │ • 构建 H(t)   │
│ • 耦合       │     │ • 载波       │     │ • 施加硬件效应   │     │ • 求解 ODE    │
│ • 坐标系     │     │ • 相位       │     │ • ZOH，低通，   │     │ • 提取可       │
│              │     │ • DRAG       │     │   量化           │     │   观测量       │
└──────────────┘     └──────────────┘     └──────────────────┘     └───────────────┘
                                                                          │
                                                                          ▼
                                                                   ┌───────────────┐
                                                                   │  仿真结果     │
                                                                   │               │
                                                                   │ • 态 States(t)│
                                                                   │ • 布居数      │
                                                                   │ • Bloch 坐标  │
                                                                   │ • Wigner 函数  │
                                                                   │ • 保真度      │
                                                                   └───────────────┘
```

**各阶段详解：**

1. **模型构造**：冻结数据类，包含所有哈密顿量参数。三种模型类封装 `UniversalCQEDModel`：

   | 模型 | 模式 | 应用场景 |
   |------|------|----------|
   | `DispersiveTransmonCavityModel` | 量子比特 + 腔体 | 单模读出或控制 |
   | `DispersiveReadoutTransmonStorageModel` | 量子比特 + 存储 + 读出 | 完整三模实验 |
   | `UniversalCQEDModel` | 任意 N 模 | 通用多模仿真 |

2. **脉冲构造**：`Pulse` 数据类，含包络（解析或预采样）、载波频率、相位、幅度、DRAG 系数和采样率。波形：$\epsilon(t) = \text{amp} \cdot \text{env}(t_{\text{rel}}) \cdot e^{i(\text{carrier} \cdot t + \text{phase})}$

3. **序列编译**：`SequenceCompiler.compile()` 接收 `Pulse` 对象列表，生成采样波形时间线，施加硬件效应：零阶保持（ZOH）、低通滤波、DAC 量化、时序量化、IQ 失真和串扰混合。

4. **仿真**：`simulate_sequence()` 组装含时哈密顿量 $H_0 + H_{\text{drive}}(t)$，添加来自 `NoiseSpec` 的 Lindblad 衰减算符，调用求解器。

5. **结果提取**：`SimulationResult` 提供偏迹、Bloch 坐标、Fock 布居数、Wigner 函数、保真度计算等。

#### 高通量执行

- `SimulationSession` / `prepare_simulation()` 预计算哈密顿量一次
- `run_many()` / `simulate_batch()` 通过 `ProcessPoolExecutor` 提供并行执行
- 会话复用模式在参数扫描中避免重复算符构造

#### 噪声建模

`NoiseSpec` 支持：

| 通道 | 参数 | 衰减算符 |
|------|------|---------|
| Transmon $T_1$ | 聚合或逐跃迁 | $\sqrt{\gamma_1} \hat{b}$ |
| Transmon 纯退相 | $T_\phi$（速率 $\gamma_\phi = 1/2T_\phi$） | $\sqrt{\gamma_\phi} \hat{\sigma}_z$ |
| 腔体损耗 | $\kappa$, $n_{\text{th}}$ | $\sqrt{\kappa(n_{\text{th}}+1)} \hat{a}$, $\sqrt{\kappa n_{\text{th}}} \hat{a}^\dagger$ |
| 腔体退相 | $T_\phi$（速率 $\gamma_\phi = 1/T_\phi$） | $\sqrt{\gamma_\phi} \hat{n}$ |

支持 `split_collapse_operators()` 在随机主方程（SME）仿真中分离受监控与未监控通道。

### 1.4 正确性与验证

#### 测试基础设施

框架包含 **58+ 测试文件**，按物理主题组织：

| 测试类别 | 文件 | 测试内容 |
|---------|------|----------|
| 自由演化与完备性检查 | `test_01` | 无驱动恒等演化、能量守恒 |
| 腔体驱动与 Kerr | `test_02` | 位移保真度、Kerr 旋转 |
| 色散与 Ramsey | `test_03` | χ 依赖相位积累 |
| χ 约定 | `test_10` | 符号约定一致性 |
| 模型不变量 | `test_11` | 参数往返、冻结数据类 |
| 脉冲语义 | `test_12` | 载波/相位/包络组合 |
| 耗散 | `test_14` | T1 衰减率、稳态光子数 |
| 收敛回归 | `test_07` | 步长收敛、截断稳定性 |
| 硬件效应 | `test_08` | 时间线编译、ZOH、量化 |
| 门库 | `test_17`, `test_42` | 100+ 门的酉正确性 |
| 三模模型 | `test_27`, `test_28` | 多模哈密顿量组装 |
| 通用模型 | `test_34` | 任意耦合规格 |
| GRAPE 最优控制 | `test_40`, `test_41`, `test_51`, `test_52` | 收敛性、梯度精度 |
| Floquet 分析 | `test_58` | 准能量、多光子共振 |
| API 完整性 | `test_37` | 所有公共符号正确导出 |

额外的专门测试子目录：`tests/analysis/`、`tests/calibration_targets/`、`tests/conventions/`、`tests/experiment/`、`tests/golden/`、`tests/quantum_algorithms/`、`tests/rl_control/`、`tests/sim/`、`tests/unitary_synthesis/`

#### 验证方法

1. **解析极限情况**：零驱动 → 恒等演化；弱耦合 → Jaynes-Cummings；谐振极限 → 精确相干态
2. **跨后端校验**：QuTiP 求解器结果与密集 NumPy/JAX 后端对比
3. **黄金文件回归**：关键仿真输出存储为参考数据；CI 检查漂移
4. **收敛扫描**：希尔伯特空间截断、时间步细化和优化器迭代次数
5. **约定对账**：专门的 `tests/conventions/` 确保所有符号定义自洽

### 1.5 智能体集成

智能体通过结构化技能系统与 `cqed_sim` 交互：

#### 智能体如何编写仿真代码

1. **先查 API**：`cqed-sim-lookup` 技能要求智能体在编写任何仿真代码前先阅读 API 参考文档。这确保智能体使用现有功能而非重复造轮子。

2. **覆盖度评估**：智能体将所需功能分类为：
   - **完全支持** → 直接使用 `cqed_sim`
   - **部分支持** → 以 `cqed_sim` 为基础扩展，记录缺失功能
   - **不支持** → 编写独立代码，记录差距，建议上游合并

3. **代码生成**：智能体编写导入 `cqed_sim` 类并遵循既定模式的 Python 脚本：
   ```python
   from cqed_sim import DispersiveTransmonCavityModel, simulate_sequence, NoiseSpec
   
   model = DispersiveTransmonCavityModel(
       omega_c=5.241e9 * 2 * np.pi,
       omega_q=6.150e9 * 2 * np.pi,
       chi=-2.84e6 * 2 * np.pi,
       ...
   )
   result = simulate_sequence(model, pulses, config)
   ```

#### 迭代优化循环

```
智能体编写脚本
       │
       ▼
智能体运行脚本  ──── 报错？ ──── 智能体读取回溯信息
       │                                    │
       │                               分类错误
       │                               (ENVIRONMENT / DEPENDENCY /
       │                                SYNTAX / RUNTIME / PHYSICS)
       │                                    │
       │                               应用修复（最多 3 次尝试）
       │                                    │
       ▼                                    ▼
产出结果        ◄──────────────── 重新运行脚本
       │
       ▼
智能体验证结果
  (完备性检查、收敛性、文献对比)
       │
       ▼
结果通过？── 否 ──▶ 记录到 IMPROVEMENTS.md，继续迭代
       │
      是
       │
       ▼
生成图表 + 报告
```

---

## 第二部分 — cqed_based_study：智能体研究层

### 2.1 课题结构

每项课题遵循标准化的目录布局：

```
studies/<study_name>/
├── README.md           ← 课题定义（目标、方法、状态）
├── IMPROVEMENTS.md     ← 动态日志：局限性、改进思路、失败尝试
├── study_state.json    ← 机器可读状态（用于智能体协调）
├── scripts/            ← Python 仿真与分析代码
├── data/               ← 原始与处理后的数值输出（.npz, .json）
├── figures/            ← 图表，同时保存 .png（300 dpi）和 .pdf（矢量）
└── report/
    ├── report.tex      ← 遵循 AGENTS.md 模板的 LaTeX 报告
    ├── references.bib  ← BibTeX 参考文献
    └── report.pdf      ← 编译后的 PDF
```

#### README.md — 课题合约

每项课题的 README 包含强制性章节，作为人类研究人员与智能体之间的合约：

| 章节 | 目的 |
|------|------|
| **问题分类** | OPT / REP / DES / ANA — 决定适用的验证检查和附录内容 |
| **研究动机** | 为何此课题重要；REP 类需链接论文 |
| **目标** | 编号的、具体的、可证伪的目标 |
| **方法** | 将使用哪些 `cqed_sim` 模块；已记录的功能差距 |
| **预期结果** | 量化的成功标准 |
| **已知局限** | 全程更新；导入报告 |
| **状态** | ACTIVE / COMPLETE / BLOCKED |

#### IMPROVEMENTS.md — 动态改进日志

此文件是当前工作与未来工作之间的桥梁。智能体在实施过程中**实时更新**，而非仅在结束时。结构如下：

- **关键缺陷（P1）** — 可能导致结果在定性上出错的问题
- **推荐改进（P2）** — 有意义的精度或范围提升
- **锦上添花（P3）** — 较低优先级的增强
- **开放问题** — 未解决的物理观测
- **已尝试但未成功的方案** — 记录死胡同以避免重复
- **计算与资源备注** — 挂钟时间、内存使用、瓶颈

### 2.2 执行流水线

#### 标准工作流

```
步骤 1：初始化  →  创建课题文件夹、README、IMPROVEMENTS.md
步骤 2：规划    →  查阅 cqed_sim API、识别差距、陈述假设
步骤 3：实施    →  编写脚本、运行仿真、保存数据、生成图表
步骤 4：验证    →  完备性检查 ✓  收敛性 ✓  文献对比 ✓
步骤 5：报告    →  编写含强制附录的 report.tex → 编译 PDF
```

#### 参数扫描与优化

课题通过多种模式调用 `cqed_sim`：

| 模式 | 应用场景 | 示例 |
|------|----------|------|
| **单次仿真** | 验证特定参数点 | `simulate_sequence(model, pulses, config)` |
| **参数扫描** | 绘制保真度景观 | 循环遍历 $\chi$、$\kappa$、幅度网格 |
| **批量执行** | 高通量并行扫描 | `SimulationSession.run_many()` |
| **优化** | 寻找最优控制脉冲 | `GrapeSolver`, `UnitarySynthesizer`, `scipy.optimize` |
| **校准** | 定向参数提取 | `cqed_sim.calibration_targets` 代理实验 |

#### 示例：灰箱自适应控制课题

该课题展示了完整的流水线：

1. **Phase 4**：在 chi 失配水平（0–40%）下比较标称、灰箱、全知和黑箱控制
2. **Phase 5**：在噪声、读出混淆、探测预算削减、漂移和哈密顿量遗漏条件下的鲁棒性压力测试
3. **验证**：12/12 项检查通过；3 个种子的多起点 GRAPE；扩展希尔伯特空间下收敛性已验证
4. **输出**：9 张图表、7 个数据文件、经验证的 LaTeX 报告

### 2.3 智能体角色

系统采用**双模型架构**，职责清晰分离：

#### 科学主管（Codex / GPT）

| 能力 | 描述 |
|------|------|
| **角色** | 科学大脑 — 推理物理，不关心实现 |
| **调用时机** | 规划阶段（设计实验）和审查阶段（评估结果） |
| **输入** | 紧凑的结构化状态（study_state.json + 图表摘要 + 结果摘要） |
| **输出** | 含有序行动项的结构化 SCIENCE_DIRECTIVE.md |
| **优势** | 物理推理、假设生成、实验设计、质量判断 |
| **不做** | 编写代码、运行仿真、调试错误 |

#### 执行工程师（Claude Opus 4.6）

| 能力 | 描述 |
|------|------|
| **角色** | 研究工程师 + 技术写作 — 处理所有实现工作 |
| **调用时机** | 引导、实施、验证和报告阶段 |
| **输入** | SCIENCE_DIRECTIVE.md + 完全文件访问权限 |
| **输出** | 代码、数据、图表、文档、EXECUTION_SUMMARY.md |
| **优势** | 代码生成、调试、文档编写、LaTeX 报告、结构化推理 |
| **不做** | 做物理判断或决定研究方向 |

#### Copilot Chat（Codex-medium）

| 能力 | 描述 |
|------|------|
| **角色** | 快速迭代助手，用于快速编辑和实验 |
| **调用时机** | 交互式编辑、快速调试、临时查询 |
| **优势** | 速度快、低延迟、适合常规代码更改 |
| **局限** | 对深度物理推理或复杂多步任务可靠性较低 |

### 2.4 输出产物

每项完成的课题产出：

| 产物 | 格式 | 用途 |
|------|------|------|
| **数据文件** | `.npz`, `.json` | 可复现的数值结果 |
| **图表** | `.png`（300 dpi）+ `.pdf`（矢量） | 出版质量图表，采用色盲友好配色 |
| **LaTeX 报告** | `report.tex` + `report.pdf` | 同行评审风格文档：摘要、方法、结果、验证、讨论、局限性、附录 |
| **IMPROVEMENTS.md** | Markdown | 面向未来智能体/研究人员的可执行交接文档 |
| **study_state.json** | JSON | 用于智能体协调的机器可读状态 |

#### 报告结构（强制）

```
摘要 → 引言 → 系统与方法 → 结果 → 验证
→ 讨论 → 结论 → 局限性与未来工作 → 参考文献 → 附录
```

附录是**必需的**，包含支持正文发现和解读的详细数据（脉冲形状、完整参数表、扫描数据、代价景观）。

### 2.5 故障处理

#### 自调试协议（4 级升级）

```
第 1 级：检查（< 30 秒）
  ├─ 读取错误回溯信息
  ├─ 分类：ENVIRONMENT | DEPENDENCY | SYNTAX | RUNTIME | PHYSICS | ASSUMPTION
  └─ 应用定向修复

第 2 级：修复（每个故障最多 3 次尝试）
  ├─ ENVIRONMENT → 检查路径、权限、Python 版本
  ├─ DEPENDENCY  → pip install --user，检查版本
  ├─ SYNTAX      → 修复代码，重新运行
  ├─ RUNTIME     → 检查数组形状、参数范围、NaN/Inf
  ├─ PHYSICS     → 检查哈密顿量、单位、量级
  └─ ASSUMPTION  → 记录为科学问题，提交主管审查

第 3 级：记录并升级（3 次尝试失败后）
  ├─ 在 BLOCKERS.md 中记录完整回溯信息
  ├─ 添加到 study_state.json 的 failed_tasks
  ├─ 继续下一个未阻塞的任务
  └─ 在 EXECUTION_SUMMARY.md 中标记以供审查

第 4 级：停止（仅在所有剩余任务均被阻塞时）
  ├─ 编写全面的阻塞报告
  ├─ 保存所有部分结果
  └─ 设置 status = "BLOCKED"
```

#### 中断恢复

基于文件的状态系统实现无缝恢复：
- `study_state.json` 跟踪已完成、失败、待处理和被阻塞的任务
- `TASK_CHECKLIST.md` 作为执行的唯一事实来源
- `PROGRESS_LOG.md` 提供仅追加的检查点
- `autonomous-resume` 智能体从文件重建状态，从中断处精确继续

---

## 第三部分 — 智能体工作流设计

### 3.1 多智能体系统设计

#### 角色分化

| 角色 | 模型 | 职责 | 调用方式 |
|------|------|------|----------|
| **科学主管** | OpenAI Codex / GPT | 物理推理、实验设计、结果审查、质量判断 | `@science-director study=... phase=plan\|review` |
| **执行工程师** | Claude Opus 4.6 | 代码编写、仿真执行、调试、文档、报告 | `@execution-engineer study=... phase=bootstrap\|implement\|validate\|report` |
| **研究循环协调器** | 组合（角色切换） | 端到端自主协调 | `@research-loop study=... goal='...'` |
| **自主规划器** | 通用 | 将任务文档转换为执行计划 | `/Autonomous Plan task=... run=...` |
| **自主实施器** | 通用 | 带状态跟踪地执行清单项目 | `/Autonomous Implement task=... run=...` |
| **自主恢复器** | 通用 | 从文件状态恢复中断的任务 | `/Autonomous Resume task=... run=...` |

#### 为什么使用双模型？

双模型分离源于互补优势：

| 维度 | 科学主管（Codex/GPT） | 执行工程师（Opus） |
|------|----------------------|-------------------|
| **物理推理** | 深度领域知识、假设生成 | 执行指令，不质疑物理判断 |
| **代码生成** | 仅伪代码 / API 引用 | 完整实现与调试 |
| **工具访问** | 只读（搜索、读取文件） | 完全访问（读、写、执行、编辑） |
| **上下文效率** | 查看摘要，不看原始数据 | 查看完整文件，编写更新 |
| **决策权限** | 做研究决策（继续/修正/验证/停止） | 执行决策，上报问题 |

这种分离防止单一模型既设计糟糕的实验又以确认偏差"验证"它们。

### 3.2 迭代循环

#### 完整研究循环协议

```
┌──────────────────────────────────────────┐
│           用户触发课题                    │
│    （目标 + 可选约束）                    │
└──────────────┬───────────────────────────┘
               ▼
┌──────────────────────────────────────────┐
│     阶段 0：引导（Opus）                  │
│  • 创建课题文件夹结构                     │
│  • 初始化 README、IMPROVEMENTS.md         │
│  • 创建 study_state.json                 │
│  • 查阅 cqed_sim API                    │
└──────────────┬───────────────────────────┘
               ▼
┌──────────────────────────────────────────┐
│     阶段 1：科学规划（Codex）             │
│  • 分类问题（OPT/REP/DES/ANA）           │
│  • 提出假设                              │
│  • 设计实验                              │
│  • 定义成功标准                           │
│  • 产出 SCIENCE_DIRECTIVE.md             │
└──────────────┬───────────────────────────┘
               ▼
┌──────────────────────────────────────────┐
│   阶段 2：实施与执行（Opus）              │
│  • 读取 SCIENCE_DIRECTIVE.md             │
│  • 编写仿真脚本                           │
│  • 运行仿真，保存数据                     │
│  • 生成图表                              │
│  • 自行调试故障                           │
│  • 编写 EXECUTION_SUMMARY.md             │
└──────────────┬───────────────────────────┘
               ▼
┌──────────────────────────────────────────┐
│    阶段 3：科学审查（Codex）              │
│  • 评估物理正确性                         │
│  • 检查结果质量                           │
│  • 识别差距和改进空间                     │
│  • 决策：                                │
│    ├─ CONTINUE → 细化、扩展              │
│    ├─ REVISE → 新方法/假设               │
│    ├─ VALIDATE → 准备进入验证             │
│    └─ STOP → 需要人工介入                │
└──────────────┬───────────────────────────┘
               ▼
         ┌── 决策 ──┐
         │          │
    CONTINUE/     VALIDATE
    REVISE           │
         │           ▼
    回到          ┌──────────────────┐
    阶段 2        │ 阶段 4：验证     │
                  │ • 完备性检查     │
                  │ • 收敛性         │
                  │ • 文献对比       │
                  └────────┬─────────┘
                           ▼
                  ┌──────────────────┐
                  │ 阶段 5：报告     │
                  │ • 编写 LaTeX     │
                  │ • 编译 PDF       │
                  │ • 最终文档       │
                  │ • 标记 COMPLETE  │
                  └──────────────────┘
```

#### 交接机制

智能体间的通信完全**基于文件** — 不共享聊天历史：

| 方向 | 文件 | 内容 |
|------|------|------|
| 科学主管 → 执行工程师 | `SCIENCE_DIRECTIVE.md` | 待执行的实验、参数、成功标准 |
| 执行工程师 → 科学主管 | `EXECUTION_SUMMARY.md` | 结果摘要、异常情况、计算备注 |
| 双方 → 状态 | `study_state.json` | 机器可读的唯一事实来源 |
| 双方 → 历史 | `PROGRESS_LOG.md` | 仅追加的执行记录 |

#### 错误检测与纠正

| 错误类型 | 检测方式 | 纠正方式 |
|---------|---------|---------|
| **代码错误** | 脚本执行的回溯信息 | 执行工程师自行调试（3 次尝试） |
| **物理错误** | 非物理值、单位错误、完备性检查失败 | 科学主管在审查中标记，发出 REVISE 指令 |
| **缺失对照** | 实验覆盖不完整 | 科学主管在审查中识别，添加到下一指令 |
| **局部最小值** | 平坦代价景观、多起点对比 | 多个初始种子、替代优化器 |
| **未收敛结果** | 参数加倍时保真度变化 > 阈值 | 增加分辨率/截断重新运行 |

#### 收敛判定

循环仅在**所有**研究质量停止标准通过时终止：

- [ ] 科学问题已有证据回答
- [ ] 物理一致性已验证（极限情况、守恒律）
- [ ] 诊断与对照完成（收敛性、完备性检查）
- [ ] 结果对小参数扰动具有鲁棒性
- [ ] 文档完整（report.tex、IMPROVEMENTS.md、所有图表）
- [ ] 开放问题已记录

### 3.3 提示策略

#### 结构化上下文而非散文

智能体接收结构化状态文件而非对话上下文：

```
SCIENCE_DIRECTIVE.md         ← 含参数的有序行动项
  ┌─ ## 课题目标
  ├─ ## 问题分类
  ├─ ## 假设
  ├─ ## 实验设计    ← 每个实验：目的、方法、参数、预期结果
  ├─ ## 执行计划    ← 编号的 [IMPLEMENT]/[RUN]/[ANALYZE]/[DOCUMENT] 任务
  ├─ ## 假设条件
  └─ ## 停止标准
```

#### 令牌效率规则

| 原则 | 实现方式 |
|------|---------|
| Codex 看摘要，不看原始数据 | EXECUTION_SUMMARY.md 限制为约 500 行 |
| Opus 接收指令，不是散文 | 含具体文件路径和参数的行动列表 |
| 状态持久化在文件中，不在上下文中 | study_state.json 是唯一事实来源 |
| 大文件从不完整传递 | 仅函数签名 + 关键结果 |
| 每次迭代自包含 | 不依赖先前聊天历史 |
| 图表用描述，不嵌入 | 文本描述 + 文件路径，非 base64 |

#### 各智能体上下文策略

| 智能体 | 获取 | 不获取 |
|--------|------|--------|
| 科学主管 | study_state.json、EXECUTION_SUMMARY.md、图表路径、结果摘要 | 完整脚本、原始数据数组、编译日志 |
| 执行工程师 | SCIENCE_DIRECTIVE.md、完全文件访问、API 参考 | 先前聊天历史、其他课题数据 |
| 研究循环 | 全部（根据当前角色切换上下文） | 无 — 自行管理上下文 |

### 3.4 AGENTS.md 系统

#### 规则如何定义

`AGENTS.md` 是所有智能体阅读的主规则文档。它定义了：

1. **5 条不可违反的规则** — 始终使用 cqed_sim、不跳过步骤、按需安装包、无虚拟环境、无笔记本
2. **问题分类** — OPT/REP/DES/ANA 分类法及典型交付物
3. **工作流步骤** — 初始化 → 规划 → 实施 → 验证 → 报告
4. **报告模板** — 含强制章节的完整 LaTeX 模板
5. **图表标准** — 300 dpi PNG + 矢量 PDF、色盲友好配色
6. **验证清单** — 三项必需检查（完备性、收敛性、文献对比）
7. **决策树** — 能否使用独立代码？是否应安装软件包？课题是否完成？

#### 智能体如何遵循约束

智能体通过 VS Code 的 `.github/` 目录配置：

```
.github/
├── agents/              ← 智能体定义文件（.agent.md）
│   ├── science-director.agent.md
│   ├── execution-engineer.agent.md
│   ├── research-loop.agent.md
│   ├── autonomous-planner.agent.md
│   ├── autonomous-implementer.agent.md
│   └── autonomous-resume.agent.md
├── prompts/             ← 斜杠命令提示文件（.prompt.md）
│   ├── autonomous-plan.prompt.md
│   ├── autonomous-implement.prompt.md
│   └── autonomous-resume.prompt.md
├── instructions/        ← 上下文相关指令（.instructions.md）
│   └── task-run-state.instructions.md
└── skills/              ← 可复用能力定义
    ├── cqed-sim-lookup/SKILL.md
    ├── latex-report/SKILL.md
    ├── publication-figures/SKILL.md
    ├── study-init/SKILL.md
    └── validate-results/SKILL.md
```

每个智能体文件指定：
- **描述**：何时使用该智能体
- **工具**：智能体可访问哪些工具（读取、搜索、编辑、执行、待办）
- **参数提示**：预期调用格式
- **系统提示**：详细的行为指令

#### 更新如何传播

当 AGENTS.md 更新时，所有智能体自动接收新规则，因为：
1. AGENTS.md 附加到工作区上下文
2. 智能体定义文件引用它（"阅读 AGENTS.md 快速参考部分"）
3. 技能引用它（"AGENTS.md 步骤 5"、"AGENTS.md 验证清单"）

---

## 第四部分 — VS Code 搭建指南（可复现）

### 4.1 所需工具

#### 必备软件

| 工具 | 版本 | 用途 | 安装方式 |
|------|------|------|---------|
| **VS Code** | 最新稳定版 | IDE 和智能体宿主 | [code.visualstudio.com](https://code.visualstudio.com) |
| **GitHub Copilot** | 最新 | 智能体基础设施（聊天、智能体、技能） | VS Code 扩展市场 |
| **Python** | 3.12.x | 仿真运行时 | [python.org](https://python.org) — 系统安装，**不使用虚拟环境** |
| **Git** | 最新 | 版本控制 | [git-scm.com](https://git-scm.com) |

#### Python 依赖

```bash
pip install --user numpy scipy matplotlib qutip pandas seaborn lmfit
```

| 软件包 | 版本要求 | 用途 |
|--------|---------|------|
| `numpy` | ≥ 1.24 | 数组运算 |
| `scipy` | ≥ 1.10 | 优化、线性代数 |
| `qutip` | ≥ 5.0 | 量子仿真后端 |
| `matplotlib` | ≥ 3.8 | 绘图 |
| `pandas` | ≥ 2.0 | 数据处理 |
| `seaborn` | （可选） | 统计可视化 |
| `lmfit` | （可选） | 曲线拟合 |

#### cqed_sim 安装

```bash
# 以可编辑模式从本地源安装
pip install --user -e /path/to/cQED_simulation
```

> **Windows 注意**：QuTiP 5.x 导入可能因 NumPy 导入期间 `platform._wmi_query` 阻塞而挂起。解决方法：在导入 qutip 前将 `platform._wmi_query` 补丁为抛出 `OSError`。仓库中包含自动应用此补丁的 `runtime_compat.py` 辅助模块。

#### 可选工具

| 工具 | 用途 |
|------|------|
| **LaTeX**（MiKTeX 或 TeX Live） | 编译 PDF 报告 |
| **JAX** | GPU 加速的 GRAPE 最优控制 |

### 4.2 模型使用策略

本工作流使用**三种不同的 AI 模型**，各自针对研究循环的特定部分优化：

#### Copilot Chat（Codex-medium） — 快速迭代器

```
┌─────────────────────────────────────────────┐
│  COPILOT CHAT（Codex-medium）               │
│                                             │
│  最适合：                                    │
│  • 快速代码编辑和修复                        │
│  • 交互式运行实验                            │
│  • 临时调试                                 │
│  • 代码补全和 IntelliSense                  │
│                                             │
│  权衡：                                      │
│  ✓ 快速（< 5 秒响应）                       │
│  ✓ 成本低（低令牌消耗）                      │
│  ✓ 适合常规编码任务                          │
│  ✗ 复杂推理可靠性较低                        │
│  ✗ 可能出现物理公式幻觉                      │
│  ✗ 多步规划能力有限                          │
│                                             │
│  在 VS Code 中使用：                         │
│  • Ctrl+I 内联编辑                          │
│  • Copilot Chat 侧边栏                      │
│  • 代码补全（Tab）                          │
└─────────────────────────────────────────────┘
```

#### OpenAI Codex（高级模式） — 物理验证器

```
┌─────────────────────────────────────────────┐
│  OPENAI CODEX（高级模式）                    │
│                                             │
│  最适合：                                    │
│  • 物理关键实现                              │
│  • 数值正确性验证                            │
│  • 实验设计与假设检验                        │
│  • 符号约定和单位检查                        │
│                                             │
│  权衡：                                      │
│  ✓ 强大的物理推理                            │
│  ✓ 善于发现数值问题                          │
│  ✓ 高风险正确性可靠                          │
│  ✗ 响应时间较慢                              │
│  ✗ 每次查询成本较高                          │
│  ✗ 某些模式下工具访问受限                    │
│                                             │
│  用作 @science-director：                    │
│  • 规划和审查阶段                            │
│  • 评估物理正确性                            │
│  • 设计数值实验                              │
└─────────────────────────────────────────────┘
```

#### Claude Opus 4.6 — 实现与文档专家

```
┌─────────────────────────────────────────────┐
│  CLAUDE OPUS 4.6                            │
│                                             │
│  最适合：                                    │
│  • 编写完整仿真脚本                          │
│  • 多步实现任务                              │
│  • LaTeX 报告生成                           │
│  • 结构化文档编写                            │
│  • 调试复杂故障                              │
│  • 文件组织与状态管理                        │
│                                             │
│  权衡：                                      │
│  ✓ 出色的长上下文推理                        │
│  ✓ 最佳可读性和结构                         │
│  ✓ 可靠的多步执行                            │
│  ✓ 完全工具访问（读/写/执行）               │
│  ✗ 琐碎编辑比 Copilot 慢                    │
│  ✗ 每会话成本较高                            │
│                                             │
│  用作 @execution-engineer：                  │
│  • 引导、实施、验证、报告                    │
│  • 所有动手研究工程                          │
└─────────────────────────────────────────────┘
```

#### 模型选择决策树

```
这是物理推理或实验设计任务吗？
├── 是 → 使用 Codex（高级模式）作为 @science-director
└── 否 → 这是多步实现或报告编写任务吗？
         ├── 是 → 使用 Opus 4.6 作为 @execution-engineer
         └── 否 → 这是快速编辑、补全或临时问题吗？
                  ├── 是 → 使用 Copilot Chat（Codex-medium）
                  └── 都是？ → 使用 @research-loop（在模式间切换）
```

### 4.3 VS Code 推荐工作流

#### 分步指南：运行首个智能体辅助课题

**1. 在 VS Code 中打开仓库**

```
文件 → 打开文件夹 → 选择 cqed_based_study 目录
```

**2. 验证工作区已被识别**

检查 VS Code 在 Copilot Chat 侧边栏中显示智能体定义。应能看到 `@science-director`、`@execution-engineer`、`@research-loop` 等可用智能体。

**3. 方案 A：全自主循环（推荐）**

在 Copilot Chat 面板中：

```
@research-loop study=studies/my_new_study goal='优化色散频移以达到 99.5% 读出保真度'
```

研究循环智能体将：
1. 创建课题文件夹结构
2. 规划实验（科学主管角色）
3. 编写并运行仿真代码（执行工程师角色）
4. 审查结果（科学主管角色）
5. 迭代直到验证通过
6. 编写最终报告

**4. 方案 B：逐步控制**

如需更多控制权，手动驱动每个阶段：

```bash
# 步骤 1：初始化课题
# 在 VS Code 中：终端 → 运行任务 → "Research: New Study"
# 或在聊天中：
@execution-engineer study=studies/my_study run=task_runs/my_study phase=bootstrap
```

```bash
# 步骤 2：科学主管规划
@science-director study=studies/my_study run=task_runs/my_study phase=plan
```

```bash
# 步骤 3：执行工程师实施
@execution-engineer study=studies/my_study run=task_runs/my_study phase=implement
```

```bash
# 步骤 4：科学主管审查
@science-director study=studies/my_study run=task_runs/my_study phase=review
```

```bash
# 步骤 5：重复步骤 3-4 直到 VALIDATE 决策，然后：
@execution-engineer study=studies/my_study run=task_runs/my_study phase=validate
@execution-engineer study=studies/my_study run=task_runs/my_study phase=report
```

**5. 检查课题状态**

```bash
# 通过 VS Code 任务：
终端 → 运行任务 → "Research: Study Status"

# 或通过 PowerShell：
.\tools\research_loop.ps1 -Action status -StudyName "my_study"
```

**6. 恢复中断的课题**

```
@research-loop study=studies/my_study resume
```

或：
```bash
.\tools\research_loop.ps1 -Action resume -StudyName "my_study"
```

### 4.4 自动化循环（高级）

#### PowerShell 研究循环脚本

`tools/research_loop.ps1` 脚本提供命令行编排：

```powershell
# 初始化新课题
.\tools\research_loop.ps1 -Action init -StudyName "chi_optimization" -StudyGoal "为读出优化 chi"

# 检查状态
.\tools\research_loop.ps1 -Action status -StudyName "chi_optimization"

# 恢复中断的课题
.\tools\research_loop.ps1 -Action resume -StudyName "chi_optimization"

# 运行特定阶段
.\tools\research_loop.ps1 -Action execute -StudyName "chi_optimization"
.\tools\research_loop.ps1 -Action validate -StudyName "chi_optimization"
.\tools\research_loop.ps1 -Action report -StudyName "chi_optimization"
```

#### VS Code 任务

通过**终端 → 运行任务**可用：

| 任务 | 功能 |
|------|------|
| **Research: New Study** | 初始化课题文件夹 + 状态文件 |
| **Research: Study Status** | 显示当前循环状态和下一步操作 |
| **Research: Resume Study** | 检测阶段并继续 |
| **Research: Run Loop Action** | 选择任意阶段运行 |
| **Copilot: Init Task Run** | 为可恢复任务引导状态文件 |
| **Copilot: Show Task Run Status** | 显示任务运行状态 |

#### 失败重试

自调试协议自动处理大多数故障。对于持续性阻塞：

1. 检查任务运行目录中的 `BLOCKERS.md`
2. 手动解决阻塞或提供指导
3. 恢复：`@research-loop study=studies/<name> resume`

#### 长时间课题执行

对于需要多次会话的课题：

1. 基于文件的状态系统（`study_state.json`、`TASK_CHECKLIST.md`、`PROGRESS_LOG.md`）在 VS Code 重启后存续
2. `autonomous-resume` 智能体从文件重建上下文，而非聊天历史
3. 每次迭代自包含 — 不依赖先前对话

### 4.5 终端 / 执行集成

#### 从 VS Code 终端运行实验

所有仿真脚本设计为可从终端运行：

```powershell
# 导航到课题脚本目录
cd studies/gray_box_adaptive_control/scripts

# 运行仿真阶段
python study_phase4.py

# 运行验证
python validate_results.py
```

#### 监控进度

- **study_state.json** — 机器可读的当前状态
- **PROGRESS_LOG.md** — 人类可读的时间顺序记录
- **TASK_CHECKLIST.md** — 基于复选框的完成度跟踪
- **控制台输出** — 实时仿真进度

#### 恢复中断的任务

```powershell
# 检查待处理项
Get-Content task_runs/my_study/TASK_CHECKLIST.md | Select-String "^\- \[ \]"

# 检查阻塞
Get-Content task_runs/my_study/BLOCKERS.md

# 通过智能体恢复
# @autonomous-resume task=... run=task_runs/my_study
```

---

## 第五部分 — 向其他领域的推广

智能体辅助研究工作流在设计上是**领域无关的** — 可以替换 cQED 特定组件，同时保留完整的编排基础设施。

### 保持不变的部分

以下组件可直接转移到任何计算物理领域：

| 组件 | 描述 | 文件 |
|------|------|------|
| **智能体基础设施** | VS Code 智能体、提示、技能、指令 | `.github/agents/`、`.github/prompts/`、`.github/skills/` |
| **课题文件夹结构** | 标准化的 `studies/<name>/` 布局 | AGENTS.md 模板 |
| **状态管理** | study_state.json、TASK_CHECKLIST.md、PROGRESS_LOG.md、BLOCKERS.md | `task_runs/<name>/` |
| **双模型循环** | 科学主管（物理推理）+ 执行工程师（实现） | RESEARCH_LOOP.md |
| **自调试协议** | 4 级升级：检查 → 修复 → 记录 → 停止 | 智能体定义 |
| **报告模板** | 含强制章节的 LaTeX（摘要到附录） | `.github/skills/latex-report/` |
| **验证框架** | 3 项检查门控：完备性、收敛性、文献 | `.github/skills/validate-results/` |
| **IMPROVEMENTS.md 模式** | 含 P1/P2/P3 优先级标签的动态改进日志 | AGENTS.md 规格 |
| **图表标准** | PNG + PDF、色盲友好、带单位的标注坐标轴 | `.github/skills/publication-figures/` |
| **自动化脚本** | PowerShell 编排、VS Code 任务 | `tools/`、`.vscode/tasks.json` |

### 因领域而变的部分

| 组件 | cQED 版本 | 所需改编 |
|------|----------|---------|
| **物理仿真器** | `cqed_sim`（基于 QuTiP） | 替换为领域特定求解器 |
| **API 参考** | cqed_sim 的 `API_REFERENCE.md` | 为新仿真器编写等效文档 |
| **技能：sim-lookup** | `cqed-sim-lookup/SKILL.md` | 使用新仿真器的 API 重写 |
| **物理约定** | 色散哈密顿量、χ、Kerr 等 | 领域特定的哈密顿量和约定 |
| **验证检查** | cQED 特定完备性检查（幺正性、色散极限） | 领域特定的极限情况 |
| **默认参数** | Transmon/腔体频率、耦合强度 | 领域特定的参数表 |
| **问题分类** | OPT/REP/DES/ANA（面向 cQED） | 可能需要领域特定分类 |
| **智能体物理知识** | 科学主管提示中的 cQED 专业知识 | 更新系统提示以包含领域专业知识 |

### 领域特定改编指南

#### 等离子体物理

| cQED 组件 | 等离子体物理替代 |
|-----------|----------------|
| `cqed_sim` | PlasmaPy、GENE 或自定义 MHD 求解器 |
| 哈密顿量构造 | MHD 方程、Vlasov-Poisson 系统 |
| 噪声建模（T1、T2、κ） | 碰撞输运、电阻率、辐射损耗 |
| 色散频移验证 | Alfvén 波色散、MHD 稳定性判据 |
| 保真度度量 | 能量守恒、动量守恒、增长率 |

#### 核理论

| cQED 组件 | 核理论替代 |
|-----------|-----------|
| `cqed_sim` | LAMMPS、Geant4、核壳模型程序 |
| QuTiP 求解器 | 多体薛定谔方程求解器、密度泛函理论 |
