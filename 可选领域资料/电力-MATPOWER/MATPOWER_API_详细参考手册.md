# MATPOWER 7.1 API 详细参考手册

> 本文档针对"基于MATPOWER + IPOPT实现多时段AC OPF（含自定义非线性约束）"的需求，
> 系统整理MATPOWER的OPF框架、opt_model API、用户扩展接口等关键函数的签名、参数说明和使用范例。
> 所有内容基于对matpower7.1源码的实际阅读。

---

## 目录

1. [总体架构与OPF调用流程](#1-总体架构与opf调用流程)
2. [高层入口函数](#2-高层入口函数)
3. [opt_model 核心API](#3-opt_model-核心api)
4. [用户函数扩展机制 (userfcn)](#4-用户函数扩展机制-userfcn)
5. [OPF标准变量、约束、成本的名称与结构](#5-opf标准变量约束成本的名称与结构)
6. [求解器调用链与IPOPT接口](#6-求解器调用链与ipopt接口)
7. [辅助函数](#7-辅助函数)
8. [结果数据结构](#8-结果数据结构)
9. [完整使用示例：通过userfcn扩展OPF](#9-完整使用示例通过userfcn扩展opf)

---

## 1. 总体架构与OPF调用流程

### 1.1 调用链总览

```
runopf(casedata, mpopt)
  └─► opf(mpc, mpopt)
        ├─ opf_args(): 解析输入参数
        ├─ ext2int(mpc): 外部→内部编号转换 + 调用userfcn('ext2int')
        ├─ opf_setup(mpc, mpopt): ★ 构建优化模型
        │    ├─ om = opf_model(mpc)                   创建模型对象
        │    ├─ om.add_var('Va'/'Vm'/'Pg'/'Qg'...)     添加标准OPF变量
        │    ├─ om.add_nln_constraint('Pmis'/'Qmis'...)  添加AC潮流等式约束
        │    ├─ om.add_nln_constraint('Sf'/'St'...)      添加支路潮流不等式约束
        │    ├─ om.add_lin_constraint(...)               添加线性约束
        │    ├─ om.add_quad_cost / add_nln_cost(...)    添加发电成本
        │    └─ run_userfcn('formulation', om, mpopt)   ★ 调用用户自定义扩展
        │
        ├─ opf_execute(om, mpopt): 调用求解器
        │    └─ nlpopf_solver(om, mpopt)
        │         └─ om.solve(opt)
        │              └─ nlps_master() → ipopt() MEX
        │
        └─ int2ext(results): 内部→外部编号 + 调用userfcn('int2ext')
```

### 1.2 关键设计理念

- **opt_model对象**（mp-opt-model库）是整个框架的核心，它封装了一个NLP问题：变量、约束、成本均以"命名块"方式组织
- **用户扩展完全通过userfcn + opt_model的API实现**，无需修改MATPOWER内部代码
- **所有命名块最终被拼接为一个大NLP**，由IPOPT（或其他求解器）一次性求解
- 变量通过名称引用，opt_model自动管理全局索引映射

---

## 2. 高层入口函数

### 2.1 `runopf` — 最外层运行入口

**文件**：`lib/runopf.m`

```matlab
results = runopf(casedata)
results = runopf(casedata, mpopt)
results = runopf(casedata, mpopt, fname)
results = runopf(casedata, mpopt, fname, solvedcase)
[results, success] = runopf(...)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `casedata` | string/struct | 算例文件名或mpc结构体，默认'case9' |
| `mpopt` | struct | MATPOWER选项，由`mpoption()`创建 |
| `fname` | string | 可选，输出到文件 |
| `solvedcase` | string | 可选，保存求解后的算例 |

**返回值**：results结构体，包含bus/gen/branch/gencost/f/success等所有标准字段加求解结果。

### 2.2 `opf` — 核心OPF函数

**文件**：`lib/opf.m`（第166-244行为核心流程）

```matlab
results = opf(mpc, mpopt)
[results, success] = opf(mpc, mpopt)
```

**核心流程**（第181-244行）：
```matlab
[mpc, mpopt] = opf_args(varargin{:});    % 第181行：解析参数
mpc = ext2int(mpc, mpopt);                % 第222行：外部→内部编号
om = opf_setup(mpc, mpopt);               % 第225行：★ 构建模型
[results, success, raw] = opf_execute(om, mpopt);  % 第232行：求解
results = int2ext(results);               % 第244行：恢复外部编号
```

### 2.3 `opf_setup` — 构建OPF模型（最关键的函数）

**文件**：`lib/opf_setup.m`（561行）

```matlab
om = opf_setup(mpc, mpopt)
```

此函数做了以下工作（按行号标注关键位置）：

1. **创建模型对象**（~第380行）：`om = opf_model(mpc);`
2. **添加标准变量**（~第391-417行）：Va, Vm, Pg, Qg
3. **添加非线性约束**（~第420-435行）：功率平衡Pmis/Qmis、支路潮流Sf/St
4. **添加线性约束**（~第438-443行）：PQ能力曲线、电压约束等
5. **添加成本**（~第445-472行）：多项式成本或分段线性成本
6. **★ 调用用户函数**（~第514行）：
   ```matlab
   om = run_userfcn(userfcn, 'formulation', om, mpopt);
   ```
   这是你的POFC约束/成本注入的入口点。

### 2.4 `opf_execute` — 执行求解

**文件**：`lib/opf_execute.m`

```matlab
[results, success, raw] = opf_execute(om, mpopt)
```

**求解器路由逻辑**（第49-136行）：
```
DC OPF    → dcopf_solver()
AC OPF:
  legacy  → tspopf_solver() (PDIPM/TRALM)
  modern  → nlpopf_solver() (MIPS/FMINCON/IPOPT/KNITRO)  ← 你会用到这个
```

### 2.5 `nlpopf_solver` — NLP求解器封装

**文件**：`lib/nlpopf_solver.m`

```matlab
[results, success, raw] = nlpopf_solver(om, mpopt)
```

**核心逻辑**（第72-295行）：
```matlab
mpc = om.get_mpc();                        % 获取算例数据
[vv, ll, nne, nni] = om.get_idx();         % 获取变量/约束索引
opt = mpopt2nlpopt(mpopt, model);          % 转换选项
[x, f, eflag, output, lambda] = om.solve(opt);  % ★ 调用求解器

% 提取解
Va = x(vv.i1.Va:vv.iN.Va);
Vm = x(vv.i1.Vm:vv.iN.Vm);
Pg = x(vv.i1.Pg:vv.iN.Pg);
Qg = x(vv.i1.Qg:vv.iN.Qg);
```

---

## 3. opt_model 核心API

opt_model是mp-opt-model库中定义的优化模型类，位于：
`mp-opt-model/lib/@opt_model/`

### 3.1 构造函数

**文件**：`mp-opt-model/lib/@opt_model/opt_model.m`

```matlab
om = opt_model()       % 空模型
om = opt_model(S)      % 从结构体S初始化
```

MATPOWER的OPF使用子类：
```matlab
om = opf_model(mpc)    % lib/@opf_model/opf_model.m，保存mpc
```

### 3.2 `add_var` — 添加优化变量

**文件**：`mp-opt-model/lib/@opt_model/add_var.m`

```matlab
om.add_var(NAME, N)
om.add_var(NAME, N, V0)
om.add_var(NAME, N, V0, VL)
om.add_var(NAME, N, V0, VL, VU)
om.add_var(NAME, N, V0, VL, VU, VT)
```

| 参数 | 类型 | 说明 | 默认值 |
|------|------|------|--------|
| `NAME` | string | 变量集名称（如'Pg', 'my_var'）| 必填 |
| `N` | integer | 变量个数 | 必填 |
| `V0` | N×1 double | 初始值 | 0 |
| `VL` | N×1 double | 下界 | -Inf |
| `VU` | N×1 double | 上界 | +Inf |
| `VT` | char/char vector | 变量类型：'C'连续, 'I'整数, 'B'二进制 | 'C' |

**注意事项**：
- 标量V0/VL/VU会自动扩展为N×1向量
- 变量名在模型中必须唯一
- 变量的添加顺序决定了它们在全局向量x中的排列顺序
- 支持索引命名集（先调用`init_indexed_name`），用于分组管理

**使用示例**：
```matlab
om.add_var('Va', nb, Va0, -pi*ones(nb,1), pi*ones(nb,1));
om.add_var('Pg', ng, Pg0, Pmin, Pmax);
om.add_var('H2', n_pofc, H2_0, H2_min, H2_max);  % 你的氢产量变量
om.add_var('SOC', n_batt*nt, SOC0, SOC_min, SOC_max); % 储能SOC变量
```

### 3.3 `add_lin_constraint` — 添加线性约束

**文件**：`mp-opt-model/lib/@opt_model/add_lin_constraint.m`

```matlab
om.add_lin_constraint(NAME, A, L, U)
om.add_lin_constraint(NAME, A, L, U, VARSETS)
```

**约束形式**：`L ≤ A * x ≤ U`

| 参数 | 类型 | 说明 |
|------|------|------|
| `NAME` | string | 约束集名称 |
| `A` | M×N sparse | 约束矩阵，N等于VARSETS中变量总维度 |
| `L` | M×1 double | 下界，空或-Inf表示无下界 |
| `U` | M×1 double | 上界，空或+Inf表示无上界 |
| `VARSETS` | cell array of strings | 涉及的变量名列表，如`{'Pg','Qg'}` |

**注意事项**：
- 若省略VARSETS，A的列数必须等于全局变量总维度
- 若指定VARSETS，A的列数 = 各变量维度之和，列按VARSETS中变量名的顺序排列
- 纯等式约束：L = U
- 纯上界约束：L = []
- 纯下界约束：U = []
- 标量L/U会自动扩展

**使用示例**：
```matlab
% 爬坡约束：-Ramp ≤ Pg(t+1) - Pg(t) ≤ Ramp
A_ramp = sparse([I, -I]);  % I是单位矩阵
om.add_lin_constraint('ramp', A_ramp, -ramp_limit, ramp_limit, {'Pg'});

% 储能SOC更新（等式）：SOC(t+1) = SOC(t) - P_discharge(t)*dt + P_charge(t)*dt
% 变形为 SOC(t+1) - SOC(t) + P_dis(t)*dt - P_chg(t)*dt = 0
A_soc = [...];  % 构建稀疏矩阵
om.add_lin_constraint('SOC_update', A_soc, zeros(M,1), zeros(M,1), {'SOC','P_dis','P_chg'});
```

### 3.4 `add_nln_constraint` — 添加非线性约束 ★最关键的API

**文件**：`mp-opt-model/lib/@opt_model/add_nln_constraint.m`

```matlab
om.add_nln_constraint(NAME, N, ISEQ, FCN, HESS)
om.add_nln_constraint(NAME, N, ISEQ, FCN, HESS, VARSETS)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `NAME` | string 或 cell | 约束名称，cell时表示多个共享FCN/HESS的约束组 |
| `N` | integer 或 vector | 约束个数（cell NAME时为各组个数的向量） |
| `ISEQ` | 0 或 1 | 1=等式约束 g(x)=0，0=不等式约束 h(x)≤0 |
| `FCN` | function_handle | 约束值及Jacobian的计算函数 |
| `HESS` | function_handle | 约束的Lagrangian Hessian计算函数 |
| `VARSETS` | cell of strings | 涉及的变量名列表 |

**★ FCN的签名**（当指定VARSETS时）：
```matlab
function G = my_fcn(X)
function [G, DG] = my_fcn(X)
```
- `X`：**cell array**，每个元素对应VARSETS中一个变量的向量
  - 例如VARSETS = {'fgo', 'p'}，则 X = {fgo_vector, p_vector}
  - 通过 `[fgo, p] = deal(X{:})` 解包
- `G`：N×1约束值向量
- `DG`：N×M Jacobian矩阵，M = VARSETS中所有变量的总维度
  - DG(i,j) = ∂G(i)/∂x(j)
  - 列按VARSETS中变量的顺序排列

**★ HESS的签名**（当指定VARSETS时）：
```matlab
function D2G = my_hess(X, LAMBDA)
```
- `X`：同上，cell array
- `LAMBDA`：N×1乘子向量
- `D2G`：M×M 矩阵，= Σ_i LAMBDA(i) * ∇²G_i(x)
  - 这是约束的加权Hessian
  - 行列都按VARSETS中变量的总维度排列
  - 必须返回对称矩阵（或至少下三角部分正确）

**注意事项**：
- 不等式约束ISEQ=0对应 h(x) ≤ 0（不是 ≥ 0）
- 等式约束ISEQ=1对应 g(x) = 0
- 如果NAME是cell（如`{'Sf','St'}`），FCN需一次返回所有组的约束值拼接向量，DG也是拼接的Jacobian
- opt_model会自动将局部Jacobian/Hessian映射到全局矩阵中的正确位置

**★ 使用示例（这是你最需要参考的模式）**：

```matlab
% 定义POFC可行域非线性约束 g(P_net, H2) ≤ 0
fcn_pofc = @(x) pofc_feasibility_fcn(x, pofc_params);
hess_pofc = @(x, lambda) pofc_feasibility_hess(x, lambda, pofc_params);
om.add_nln_constraint('pofc_region', n_constraints, 0, fcn_pofc, hess_pofc, {'Pg', 'H2'});

% 其中：
function [G, DG] = pofc_feasibility_fcn(X, params)
    [Pg, H2] = deal(X{:});   % 解包变量
    % ... 计算约束值 G(N×1) 和 Jacobian DG(N×(ng+n_h2))
end

function D2G = pofc_feasibility_hess(X, lambda, params)
    [Pg, H2] = deal(X{:});
    % ... 计算加权Hessian D2G((ng+n_h2)×(ng+n_h2))
end
```

### 3.5 `add_quad_cost` — 添加二次成本

**文件**：`mp-opt-model/lib/@opt_model/add_quad_cost.m`

```matlab
om.add_quad_cost(NAME, Q, C, K, VARSETS)
```

**成本形式**：

- 当Q是N×N矩阵时：`f(x) = ½ x'Qx + c'x + k`
- 当Q是N×1向量或空时：`f(x) = ½ Q.*x² + c.*x + k`（逐元素）
- 当Q为空（[]）时：`f(x) = c'x + k`（纯线性成本）

| 参数 | 类型 | 说明 |
|------|------|------|
| `NAME` | string | 成本名称 |
| `Q` | N×N sparse 或 N×1 或 [] | 二次系数 |
| `C` | N×1 | 线性系数 |
| `K` | scalar | 常数项，默认0 |
| `VARSETS` | cell of strings | 涉及的变量名 |

**使用示例**：
```matlab
% 线性气井开采成本：cost = gcost' * g
om.add_quad_cost('gas_cost', [], gcost, 0, {'g'});

% 二次发电成本：cost = ½ Pg' * Q * Pg + c' * Pg + k
om.add_quad_cost('gen_cost', Q_gen, c_gen, k_gen, {'Pg'});

% 氢收益（负成本）：cost = -price_H2' * H2
om.add_quad_cost('H2_revenue', [], -price_H2, 0, {'H2'});

% 碳交易成本：cost = carbon_price * emission_factor' * Pg
om.add_quad_cost('carbon_cost', [], carbon_price * emission_factor, 0, {'Pg'});
```

### 3.6 `add_nln_cost` — 添加非线性成本

**文件**：`mp-opt-model/lib/@opt_model/add_nln_cost.m`

```matlab
om.add_nln_cost(NAME, N, FCN, VARSETS)
```

| 参数 | 说明 |
|------|------|
| `NAME` | 成本名称 |
| `N` | 必须等于1（标量成本） |
| `FCN` | 成本函数句柄 |
| `VARSETS` | 涉及的变量名 |

**FCN签名**：
```matlab
function F = cost_fcn(X)            % 仅值
function [F, DF] = cost_fcn(X)      % 值 + 梯度
function [F, DF, D2F] = cost_fcn(X) % 值 + 梯度 + Hessian
```
- `F`：标量，成本值
- `DF`：NX×1 梯度向量
- `D2F`：NX×NX Hessian矩阵

**使用示例**：
```matlab
% 三次及以上多项式发电成本（MATPOWER内部对高阶多项式使用此方式）
cost_Pg = @(x) opf_gen_cost_fcn(x, baseMVA, pcost, ig, mpopt);
om.add_nln_cost('polPg', 1, cost_Pg, {'Pg'});
```

### 3.7 `get_idx` — 获取索引信息

**文件**：`mp-opt-model/lib/@opt_model/get_idx.m`

```matlab
vv = om.get_idx()                              % 仅变量索引
[vv, ll] = om.get_idx()                         % + 线性约束索引
[vv, ll, nne] = om.get_idx()                    % + 非线性等式约束索引
[vv, ll, nne, nni] = om.get_idx()               % + 非线性不等式约束索引
[vv, ll, nne, nni, qq] = om.get_idx()           % + 二次成本索引
[vv, ll, nne, nni, qq, nnc] = om.get_idx()      % + 非线性成本索引
```

**每个索引结构体包含**：
```
vv.i1.NAME  — 变量NAME在全局向量x中的起始索引
vv.iN.NAME  — 变量NAME在全局向量x中的结束索引
vv.N.NAME   — 变量NAME的维度
```

**使用示例**（在求解后提取结果）：
```matlab
[vv, ll, nne, nni] = om.get_idx();
Pg = x(vv.i1.Pg:vv.iN.Pg);        % 提取发电出力
H2 = x(vv.i1.H2:vv.iN.H2);        % 提取氢产量
SOC = x(vv.i1.SOC:vv.iN.SOC);      % 提取储能SOC

% 提取线性约束乘子
mu_ramp = lambda.mu_l(ll.i1.ramp:ll.iN.ramp);
```

### 3.8 `params_var` — 获取变量参数

**文件**：`mp-opt-model/lib/@opt_model/params_var.m`

```matlab
[v0, vl, vu] = om.params_var()             % 全部变量
[v0, vl, vu] = om.params_var(NAME)          % 指定变量集
[v0, vl, vu, vt] = om.params_var(...)       % 含变量类型
```

### 3.9 `get_mpc` — 获取算例数据

**文件**：`lib/@opf_model/opf_model.m`

```matlab
mpc = om.get_mpc()
```

仅在opf_model（opt_model的子类）中可用，返回创建模型时保存的mpc结构体。

### 3.10 `solve` — 求解优化问题

**文件**：`mp-opt-model/lib/@opt_model/solve.m`

```matlab
[x, f, exitflag, output, lambda] = om.solve()
[x, f, exitflag, output, lambda] = om.solve(opt)
```

**opt选项结构体**：
```matlab
opt.alg       % 求解器：'DEFAULT', 'MIPS', 'IPOPT', 'FMINCON', 'KNITRO'
opt.verbose   % 0=静默, 1=摘要, 2=详细
opt.x0        % 初始点（覆盖add_var中的V0）
opt.ipopt_opt % IPOPT特定选项（struct）
opt.mips_opt  % MIPS特定选项
```

**返回值**：

| 返回值 | 说明 |
|--------|------|
| `x` | NX×1 最优解向量 |
| `f` | 最优目标函数值 |
| `exitflag` | 1=收敛, ≤0=失败 |
| `output` | 求解器信息，含`.alg`字段 |
| `lambda` | 乘子结构体 |

**lambda结构体**：
```matlab
lambda.eqnonlin    % 非线性等式约束乘子
lambda.ineqnonlin  % 非线性不等式约束乘子
lambda.mu_l        % 线性约束下界乘子
lambda.mu_u        % 线性约束上界乘子
lambda.lower       % 变量下界乘子
lambda.upper       % 变量上界乘子
```

### 3.11 `display` — 显示模型信息

```matlab
om          % 在MATLAB命令行直接输入模型对象名
display(om) % 显示所有变量块、约束块、成本块的摘要
```

---

## 4. 用户函数扩展机制 (userfcn)

### 4.1 `add_userfcn` — 注册回调函数

**文件**：`lib/add_userfcn.m`

```matlab
mpc = add_userfcn(mpc, stage, fcn)
mpc = add_userfcn(mpc, stage, fcn, args)
```

| 参数 | 说明 |
|------|------|
| `mpc` | MATPOWER算例结构体 |
| `stage` | 回调阶段名称（string） |
| `fcn` | 回调函数句柄 |
| `args` | 可选，传递给回调的额外参数 |

**5个回调阶段及其签名**：

| 阶段 | 调用时机 | 回调签名 |
|------|----------|----------|
| `'ext2int'` | 外部→内部编号转换后 | `mpc = fcn(mpc, mpopt, args)` |
| `'formulation'` | opf_setup中，标准OPF模型构建完成后 | `om = fcn(om, mpopt, args)` |
| `'int2ext'` | 求解完成，内部→外部编号转换前 | `results = fcn(results, mpopt, args)` |
| `'printpf'` | 标准输出打印完成后 | `results = fcn(results, fd, mpopt, args)` |
| `'savecase'` | 保存算例文件时 | `mpc = fcn(mpc, fd, prefix, args)` |

**注册顺序**：同一阶段可注册多个回调，按注册顺序依次调用。

**使用示例**：
```matlab
% 注册POFC扩展的三个核心回调
mpc = add_userfcn(mpc, 'ext2int', @userfcn_pofc_ext2int, pofc_args);
mpc = add_userfcn(mpc, 'formulation', @userfcn_pofc_formulation, pofc_args);
mpc = add_userfcn(mpc, 'int2ext', @userfcn_pofc_int2ext);
```

### 4.2 `run_userfcn` — 执行回调

**文件**：`lib/run_userfcn.m`

```matlab
rv = run_userfcn(userfcn, stage, varargin)
```

内部遍历 `userfcn.(stage)` 数组，依次调用每个注册的回调。通常不需要用户直接调用。

### 4.3 formulation回调的完整模式

```matlab
function om = userfcn_pofc_formulation(om, mpopt, args)
    % 1. 从om中获取电力系统数据
    mpc = om.get_mpc();

    % 2. 获取问题维度
    nb = size(mpc.bus, 1);
    ng = size(mpc.gen, 1);

    % 3. 读取自定义数据（存储在mpc的自定义字段中）
    pofc_data = mpc.pofc;

    % 4. 添加新的优化变量
    om.add_var('H2', n_h2, H2_0, H2_min, H2_max);
    om.add_var('SOC', n_soc, SOC_0, SOC_min, SOC_max);

    % 5. 添加线性约束
    om.add_lin_constraint('soc_balance', A_soc, L_soc, U_soc, {'SOC', 'P_batt'});

    % 6. 添加非线性约束
    fcn = @(x) pofc_constraint_fcn(x, pofc_params);
    hess = @(x, lam) pofc_constraint_hess(x, lam, pofc_params);
    om.add_nln_constraint('pofc_feasible', n_con, 0, fcn, hess, {'Pg', 'H2'});

    % 7. 添加成本
    om.add_quad_cost('H2_revenue', [], -H2_price, 0, {'H2'});
    om.add_quad_cost('carbon_cost', [], carbon_c, 0, {'Pg'});
end
```

---

## 5. OPF标准变量、约束、成本的名称与结构

### 5.1 标准AC OPF变量

在opf_setup.m中添加的标准变量及其在全局向量x中的排列顺序：

| 顺序 | 变量名 | 维度 | 含义 | 单位（标幺值） |
|------|--------|------|------|----------------|
| 1 | `Va` | nb | 节点电压相角 | rad |
| 2 | `Vm` | nb | 节点电压幅值 | p.u. |
| 3 | `Pg` | ng | 发电机有功出力 | p.u. (baseMVA) |
| 4 | `Qg` | ng | 发电机无功出力 | p.u. (baseMVA) |
| (5) | `y` | ny | 分段线性成本辅助变量 | - |

用户变量（通过userfcn添加）排在标准变量之后。

**重要：MATPOWER内部所有功率值都是标幺值**，以baseMVA为基准。
- `Pg = Pg_MW / baseMVA`
- 在userfcn中引用`{'Pg'}`时，得到的是标幺值
- gencost中的成本系数已经考虑了这个标幺化

### 5.2 标准AC OPF约束

| 约束名 | 类型 | 维度 | 含义 |
|--------|------|------|------|
| `Pmis` | nle (等式) | nb | 有功功率平衡 |
| `Qmis` | nle (等式) | nb | 无功功率平衡 |
| `Sf` | nli (不等式) | nl2 | 支路首端潮流限制 |
| `St` | nli (不等式) | nl2 | 支路末端潮流限制 |
| `PQh` | lin | - | PQ能力曲线上限 |
| `PQl` | lin | - | PQ能力曲线下限 |
| `vl` | lin | - | 恒功率因数负荷 |
| `ang` | lin | nang | 电压角差限制 |
| `ycon` | lin | ny | PWL成本基底约束 |

### 5.3 标准成本

| 成本名 | 类型 | 含义 |
|--------|------|------|
| `polPg` | quad_cost 或 nln_cost | 有功发电成本（二次或高阶多项式） |
| `polQg` | quad_cost 或 nln_cost | 无功发电成本 |
| `pwl` | quad_cost | 分段线性成本 |

---

## 6. 求解器调用链与IPOPT接口

### 6.1 完整调用链

```
om.solve(opt)
  └─ opt_model/solve.m
       ├─ 检测问题类型（NLP/QP/LP/...）
       ├─ 组装目标函数和约束的回调
       └─ nlps_master(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt)
            └─ 根据opt.alg路由到具体求解器：
                 ├─ nlps_ipopt()  → ipopt() MEX接口  ← IPOPT
                 ├─ nlps_fmincon() → fmincon()
                 ├─ nlps_knitro() → ktrlink()
                 └─ mips()        → MATPOWER内置Interior Point Solver
```

### 6.2 nlps_master的NLP问题格式

**文件**：`mp-opt-model/lib/nlps_master.m`

```
min   F(X)
 X
s.t.  G(X) = 0          (非线性等式约束)
      H(X) ≤ 0          (非线性不等式约束)
      L ≤ A*X ≤ U       (线性约束)
      XMIN ≤ X ≤ XMAX   (变量界)
```

**传给nlps_master的回调签名**：

```matlab
% 目标函数
[F, DF, D2F] = f_fcn(X)

% 约束（等式+不等式拼接）
[H, G, DH, DG] = gh_fcn(X)

% Lagrangian的Hessian
LXX = hess_fcn(X, LAM)
% 其中LAM = struct('eqnonlin', λ_eq, 'ineqnonlin', λ_ineq)
```

opt_model内部自动将所有注册的变量/约束/成本块组装成这些回调函数。

### 6.3 IPOPT选项设置

```matlab
mpopt = mpoption();
mpopt.opf.ac.solver = 'IPOPT';
mpopt.ipopt.opts.max_iter = 100000;
mpopt.ipopt.opts.tol = 1e-6;
mpopt.ipopt.opts.acceptable_tol = 1e-4;
% 其他IPOPT选项可通过mpopt.ipopt.opts.xxx设置
```

---

## 7. 辅助函数

### 7.1 `ext2int` / `int2ext` — 编号转换

**文件**：`lib/ext2int.m`, `lib/int2ext.m`

```matlab
mpc = ext2int(mpc)         % 外部→内部编号（连续1:nb）
mpc = ext2int(mpc, mpopt)  % 含选项，会调用userfcn('ext2int')

results = int2ext(results)  % 内部→外部编号，会调用userfcn('int2ext')
```

ext2int做的事：
- 将非连续的bus编号重映射为1:nb
- 删除离线设备和孤立节点
- 保存映射信息在`mpc.order`中
- 调用所有注册的`ext2int`阶段回调

### 7.2 `modcost` — 修改发电成本

**文件**：`lib/modcost.m`

```matlab
gencost = modcost(gencost, alpha)            % 默认SCALE_F
gencost = modcost(gencost, alpha, modtype)
```

| modtype | 效果 |
|---------|------|
| `'SCALE_F'`(默认) | F_α(X) = α * F(X)，成本乘以α |
| `'SCALE_X'` | F_α(αX) = F(X)，等效于X缩放α倍 |
| `'SHIFT_F'` | F_α(X) = F(X) + α |
| `'SHIFT_X'` | F_α(X+α) = F(X) |

**在MPNG中的用途**：`multi_period.m`中用`modcost(gencost, time(j))`将发电成本乘以时段持续时间（小时），使目标函数的单位变为$/day。

### 7.3 `extract_islands` — 提取电网岛

**文件**：`lib/extract_islands.m`

```matlab
mpc_array = extract_islands(mpc)     % 返回cell array，每个元素是一个岛
mpc_k = extract_islands(mpc, k)      % 提取第k个岛
```

MPNG在int2ext阶段用此函数将堆叠的多时段结果拆分回各时段。

### 7.4 `mpoption` — 创建/修改选项

```matlab
mpopt = mpoption()                          % 默认选项
mpopt = mpoption('opf.ac.solver', 'IPOPT')  % 设置AC OPF求解器
mpopt = mpoption(mpopt, 'verbose', 2)       % 修改已有选项
```

**常用选项**：

| 选项 | 说明 | 常用值 |
|------|------|--------|
| `opf.ac.solver` | AC OPF求解器 | 'IPOPT', 'MIPS', 'FMINCON' |
| `verbose` | 输出详细程度 | 0, 1, 2, 3 |
| `ipopt.opts.max_iter` | IPOPT最大迭代 | 100000 |
| `ipopt.opts.tol` | IPOPT收敛容差 | 1e-6 |
| `opf.start` | 初始点策略 | 0(默认), 2(使用算例值) |
| `out.all` | 输出控制 | -1(抑制), 0(默认), 1(全部) |

---

## 8. 结果数据结构

runopf/opf返回的results结构体：

```matlab
results.success     % 1=收敛, 0=失败
results.f           % 最优目标函数值
results.et          % 求解时间(秒)

% 标准MATPOWER数据（含求解结果）
results.bus         % bus矩阵，含求解后的VM, VA, LAM_P, LAM_Q等
results.gen         % gen矩阵，含求解后的PG, QG, MU_PMAX等
results.branch      % branch矩阵，含求解后的PF, QF, PT, QT等
results.gencost     % 发电成本矩阵

% 按命名块组织的变量结果
results.var.val.Va  % 电压角（标幺rad）
results.var.val.Vm  % 电压幅值（标幺）
results.var.val.Pg  % 有功出力（标幺，需乘baseMVA得MW）
results.var.val.Qg  % 无功出力（标幺）
results.var.val.XXX % 用户自定义变量XXX的最优值
results.var.mu.l.XXX % 变量XXX下界乘子
results.var.mu.u.XXX % 变量XXX上界乘子

% 按命名块组织的约束乘子
results.lin.mu.l.XXX  % 线性约束XXX的下界乘子
results.lin.mu.u.XXX  % 线性约束XXX的上界乘子
results.nle.lambda.Pmis % 功率平衡约束乘子（节点边际电价）
results.nle.lambda.Qmis
results.nli.mu.Sf    % 支路潮流约束乘子
results.nli.mu.St

% 按命名块组织的成本结果
results.qdc.XXX     % 二次成本XXX的值

% 原始求解器输出
results.raw.xr      % 原始求解器的解向量
results.raw.pimul   % 原始乘子向量
results.raw.info    % 求解器退出代码
results.raw.output  % 求解器详细输出
```

---

## 9. 完整使用示例：通过userfcn扩展OPF

以下示例展示如何在MATPOWER的AC OPF基础上，通过userfcn机制添加自定义变量、非线性约束和成本——这正是你需要做的事情的最小化完整模板。

```matlab
%% ==================== 主入口文件 ====================
function results = my_extended_opf(mpc, my_data, mpopt)
    % 设置默认选项
    if nargin < 3
        mpopt = mpoption();
        mpopt.opf.ac.solver = 'IPOPT';
        mpopt.ipopt.opts.max_iter = 1e5;
    end

    % 将自定义数据附加到mpc结构体
    mpc.my_data = my_data;

    % 注册userfcn回调
    mpc = add_userfcn(mpc, 'ext2int', @my_ext2int);
    mpc = add_userfcn(mpc, 'formulation', @my_formulation);
    mpc = add_userfcn(mpc, 'int2ext', @my_int2ext);

    % 调用标准OPF
    results = runopf(mpc, mpopt);
end

%% ==================== ext2int回调 ====================
function mpc = my_ext2int(mpc, mpopt, args)
    % 数据验证、标幺化等
    % mpc.my_data仍然可访问
    my_data = mpc.my_data;

    % 示例：验证维度
    ng = size(mpc.gen, 1);
    if size(my_data.param, 1) ~= ng
        error('维度不匹配');
    end
end

%% ==================== formulation回调（核心）====================
function om = my_formulation(om, mpopt, args)
    mpc = om.get_mpc();
    my_data = mpc.my_data;

    ng = size(mpc.gen, 1);
    baseMVA = mpc.baseMVA;

    %% --- 添加自定义变量 ---
    n_h2 = my_data.n_h2;
    om.add_var('H2', n_h2, zeros(n_h2,1), zeros(n_h2,1), my_data.H2_max);

    %% --- 添加线性约束 ---
    % 示例：H2总产量上限
    A_h2_total = ones(1, n_h2);
    om.add_lin_constraint('H2_total', A_h2_total, [], my_data.H2_total_max, {'H2'});

    %% --- 添加非线性约束 ---
    % POFC可行域约束：涉及Pg和H2
    params.coeff_a = my_data.coeff_a;
    params.coeff_b = my_data.coeff_b;
    n_nln = 2 * n_h2;  % 例如每个POFC单元2个约束

    fcn = @(x) my_nln_constraint_fcn(x, params);
    hess = @(x, lam) my_nln_constraint_hess(x, lam, params);
    om.add_nln_constraint('pofc', n_nln, 0, fcn, hess, {'Pg', 'H2'});
    % 0表示不等式约束 h(x) ≤ 0

    %% --- 添加成本 ---
    % 氢收益（目标函数中为负成本）
    om.add_quad_cost('H2_revenue', [], -my_data.H2_price * ones(n_h2,1), 0, {'H2'});
end

%% ==================== 非线性约束函数 ====================
function [G, DG] = my_nln_constraint_fcn(X, params)
    [Pg, H2] = deal(X{:});  % 解包：X{1}=Pg向量, X{2}=H2向量

    n_pg = length(Pg);
    n_h2 = length(H2);
    n_con = 2 * n_h2;  % 约束总数

    % 计算约束值 G (n_con × 1)
    % 示例：H2(i)^2 + a*Pg(i) - b ≤ 0 和 -H2(i) ≤ 0
    G = zeros(n_con, 1);
    G(1:n_h2) = H2.^2 + params.coeff_a .* Pg(1:n_h2) - params.coeff_b;
    G(n_h2+1:end) = -H2;

    if nargout > 1
        % 计算Jacobian DG (n_con × (n_pg + n_h2))
        % DG(i,j) = dG(i)/dx(j)
        dG_dPg = sparse(1:n_h2, 1:n_h2, params.coeff_a, n_con, n_pg);
        dG_dH2 = sparse([1:n_h2, n_h2+1:n_con], [1:n_h2, 1:n_h2], ...
                        [2*H2; -ones(n_h2,1)], n_con, n_h2);
        DG = [dG_dPg, dG_dH2];
    end
end

%% ==================== 非线性约束Hessian ====================
function D2G = my_nln_constraint_hess(X, lambda, params)
    [Pg, H2] = deal(X{:});

    n_pg = length(Pg);
    n_h2 = length(H2);
    N = n_pg + n_h2;  % 总变量维度

    % D2G = Σ_i lambda(i) * ∇²G_i
    % 在此例中只有 d²G/dH2² = 2*lambda(1:n_h2) 是非零的
    hess_H2 = sparse(1:n_h2, 1:n_h2, 2*lambda(1:n_h2), n_h2, n_h2);

    D2G = [sparse(n_pg, n_pg),   sparse(n_pg, n_h2);
           sparse(n_h2, n_pg),   hess_H2           ];
end

%% ==================== int2ext回调 ====================
function results = my_int2ext(results, mpopt, args)
    % 从results.var.val中提取自定义变量结果
    H2 = results.var.val.H2;

    % 可以将结果存储到results的自定义字段
    results.my_results.H2 = H2;
    results.my_results.H2_revenue = results.qdc.H2_revenue;
end
```

---

## 附录A：MATPOWER标准数据矩阵列定义速查

### bus矩阵（常用列）
```
BUS_I=1, BUS_TYPE=2, PD=3, QD=4, GS=5, BS=6, BUS_AREA=7,
VM=8, VA=9, BASE_KV=10, ZONE=11, VMAX=12, VMIN=13
```

### gen矩阵（常用列）
```
GEN_BUS=1, PG=2, QG=3, QMAX=4, QMIN=5, VG=6, MBASE=7,
GEN_STATUS=8, PMAX=9, PMIN=10, RAMP_30=23
```

### branch矩阵（常用列）
```
F_BUS=1, T_BUS=2, BR_R=3, BR_X=4, BR_B=5, RATE_A=6,
BR_STATUS=11, PF=14, QF=15, PT=16, QT=17
```

### gencost矩阵
```
MODEL=1 (1=PW_LINEAR, 2=POLYNOMIAL), STARTUP=2, SHUTDOWN=3,
NCOST=4, COST=5:end
```

---

## 附录B：关键注意事项总结

1. **标幺值**：MATPOWER内部所有功率值都除以baseMVA，在userfcn中操作Pg/Qg等变量时要注意单位
2. **Jacobian方向**：`add_nln_constraint`的FCN返回的DG是 N_con × N_var 的矩阵（行=约束，列=变量）
3. **不等式方向**：`add_nln_constraint` ISEQ=0 对应 h(x) ≤ 0
4. **Hessian对称性**：HESS返回的D2G必须是对称的（至少下三角正确），IPOPT要求这一点
5. **稀疏矩阵**：所有矩阵应使用sparse格式，对大规模问题至关重要
6. **VARSETS的列顺序**：当指定VARSETS = {'A','B'}时，Jacobian和Hessian的列（行）顺序为 [A的变量, B的变量]
7. **多时段堆叠后的Pg**：如果使用岛堆叠实现多时段，Pg = [Pg_t1; Pg_t2; ...; Pg_nt]，每个时段ng个发电机
8. **成本函数的时间修正**：multi_period中用modcost将gencost乘以时段小时数，确保目标函数单位一致
