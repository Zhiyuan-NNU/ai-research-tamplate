# MPNG (MATPOWER-Natural Gas) 代码结构完整分析报告

> **分析目的**：为在MATPOWER框架上实现"多时段AC OPF + POFC氢电联产可行域约束"提供技术参考
> **分析对象**：mpng-master (v0.99beta) + matpower7.1
> **分析方法**：逐文件源码阅读，以代码实际实现为准

---

## 目录

1. [MPNG代码文件清单及各文件功能说明](#1-mpng代码文件清单及各文件功能说明)
2. [核心数据流：从输入到求解器](#2-核心数据流从输入到求解器)
3. [问题1：多时段堆叠机制详解](#3-问题1多时段堆叠机制详解)
4. [问题2：与MATPOWER的接口方式详解](#4-问题2与matpower的接口方式详解)
5. [问题3：目标函数的扩展方式详解](#5-问题3目标函数的扩展方式详解)
6. [问题4：数据结构详解](#6-问题4数据结构详解)
7. [问题5：求解器调用方式详解](#7-问题5求解器调用方式详解)
8. [路径评估："在MPNG上修改" vs "在MATPOWER上重写"](#8-路径评估在mpng上修改-vs-在matpower上重写)
9. [对POFC项目的具体参考价值与局限性分析](#9-对pofc项目的具体参考价值与局限性分析)

---

## 1. MPNG代码文件清单及各文件功能说明

### 1.1 Functions/ 目录（核心代码，共13个.m文件）

| 文件名 | 行数 | 功能说明 |
|--------|------|----------|
| **`mpng.m`** | 1646 | **主入口+核心逻辑**。前125行为主函数（注册userfcn → 数据预处理 → runopf）；第126-292行为`userfcn_mpng_ext2int`（数据验证+标幺化）；第293-878行为`userfcn_mpng_formulation`（★核心：添加天然气变量/约束/成本到opt_model）；第879-1023行为`userfcn_mpng_int2ext`（结果提取+反标幺化）；第1024-1394行为`userfcn_mpng_printpf`+`savecase`（结果打印）；第1395-1646行为**6个非线性约束的Jacobian/Hessian函数**（fpipe_fcn/hess, wcompgas_fcn/hess, compgas_fcn/hess, wcomppower_fcn/hess） |
| **`multi_period.m`** | 171 | **多时段堆叠**：将单时段mpc复制nt份，创建互不相连的"岛"——bus编号偏移nb*i, gen/branch编号对应偏移，负荷替换为各时段的pd/qd，发电成本乘以时段持续小时数。**这是你最可能直接复用的文件**。 |
| **`mpc2gas_prep.m`** | 112 | **数据预处理总调度**：依次调用 compressor2gen() → multi_period() → nsd2gen()，还处理了Unit Commitment对发电上下限的修改。 |
| **`compressor2gen.m`** | 92 | 将电力驱动的压缩机建模为mpc.gen中的"负发电机"（PMIN=-inf, PMAX=0），代表压缩机的电力消耗。同时创建零成本的gencost行。压缩机的ID记录在`mpc.genid.comp`中。 |
| **`nsd2gen.m`** | 109 | 将所有正负荷（PD>0的bus）转为"可调度负荷"——在gen矩阵中添加负发电机，PG=-PD, PMIN=-PD, PMAX=0。成本为切负荷惩罚（默认5000$/MWh）。原始bus的PD/QD清零。这样切负荷就成为优化变量。记录在`mpc.nsd`和`mpc.genid.dl`中。 |
| **`wey_approx.m`** | 82 | **Weymouth方程及其导数**：管道气流 f = sgn(Δπ)·sqrt(K·|Δπ|)，在零点附近（|Δπ|<π*）用5次多项式平滑替代以避免不可微。返回函数值f、一阶导df、二阶导d2f。**对你无直接用途**（你没有管道流方程），但其平滑近似思路可能对POFC可行域边界的平滑处理有参考价值。 |
| **`mgc_PU.m`** | 118 | 天然气数据标幺化：压力→平方除以pbase²，流量除以fbase，成本乘以fbase等。**注意MPNG对压力变量使用的是"压力的平方"的标幺值**（π = p²/p_base²），这使得Weymouth方程在标幺空间中更简洁。 |
| **`mgc_REAL.m`** | （对称于mgc_PU） | 标幺值反变换为实际值，在int2ext阶段使用。 |
| **`idx_node.m`** | ~60 | 定义天然气节点数据矩阵的列索引常量：DEM=1, WELL=2, NODE_I=3, NODE_TYPE=4, PR=5, PRMAX=6, PRMIN=7, OVP=8, UNP=9, COST_OVP=10, COST_UNP=11, GD=12, NGD=13 |
| **`idx_pipe.m`** | ~40 | 管道列索引：F_NODE=1, T_NODE=2, FG_O=3, K_O=4, DIAM=5, LNG=6, FMAX_O=7, FMIN_O=8, COST_O=9 |
| **`idx_comp.m`** | ~50 | 压缩机列索引：COMP_G=1, COMP_P=2, F_NODE=3, T_NODE=4, TYPE_C=5, RATIO_C=6, B_C=7, Z_C=8, ALPHA=9, BETA=10, GAMMA=11, FMAX_C=12, COST_C=13, FG_C=14, GC_C=15, PC_C=16 |
| **`idx_well.m`** | ~35 | 气井列索引：WELL_NODE=1, G=2, PW=3, GMAX=4, GMIN=5, WELL_STATUS=6, COST_G=7 |
| **`idx_sto.m`** | ~60 | 储气列索引：STO_NODE=1, STO=2, STO_0=3, STOMAX=4, STOMIN=5, FSTO=6, FSTO_OUT=7, FSTO_IN=8, FSTOMAX=9, FSTOMIN=10, S_STATUS=11, COST_STO=12, COST_OUT=13, COST_IN=14 |
| **`idx_connect.m`** | ~50 | 互联数据列索引：GEN_ID=1, MAX_ENER=2, COMP_ID=1, BUS_ID=2, NODE_ID=2, EFF=3 |
| **`define_constants_gas.m`** | 57 | 便捷脚本，一次性调用所有idx_xxx函数定义常量 |

### 1.2 Cases/ 目录（算例文件）

| 文件名 | 说明 |
|--------|------|
| `case2.m` | 2节点电力系统（用于纯天然气分析时的占位电网） |
| `case9_new.m` | 修改版9节点电力系统（IEEE 9 bus的MPNG版本） |
| `ng_case8.m` | 8节点天然气系统（1个气井、6条管道、2个压缩机） |
| `ng_case48.m` | 48节点天然气系统（大规模，用于Example2和Example4） |
| `connect_pg_case2.m` | 2-bus电力与天然气的互联数据 |
| `connect_pg_case9.m` | 9-bus电力与8-node天然气的互联数据 |
| `connect_pg_case118.m` | 118-bus电力与48-node天然气的互联数据 |

### 1.3 Examples/ 目录

| 文件名 | 场景 | 说明 |
|--------|------|------|
| `Example1.m` | 9-bus + 8-node, 4时段 | 基础示例：展示电气耦合OPF的标准用法 |
| `Example2.m` | 2-bus + 48-node, 单时段 | 纯天然气优化：管道N-1分析（绘制对比图） |
| `Example3.m` | 9-bus + 8-node, 4时段 | 涡轮膨胀机（压缩比<1）：3个场景对比 |
| `Example4.m` | 118-bus + 48-node, 24时段 | ★大规模：24小时机组组合+天然气耦合优化（883行，含UC） |

### 1.4 MPNG_User's_Manual/ 目录

LaTeX源文件，含以下章节：
- Introduction.tex
- Getting_started.tex
- Natural_Gas_Flow.tex
- **Formulation.tex** ← 数学公式推导
- Examples.tex
- Appendixes.tex

---

## 2. 核心数据流：从输入到求解器

### 2.1 完整流程图

```
用户构造输入数据
  mpgc = mpc;                    % 标准MATPOWER电力算例
  mpgc.mgc = mgc;                % 天然气系统数据
  mpgc.connect = connect;         % 电-气互联数据
         │
         ▼
mpng(mpgc, mpopt)                 [mpng.m 第1-125行]
  │
  ├─ 默认设置 mpopt.opf.ac.solver = 'IPOPT'    [第80行]
  │
  ├─ add_userfcn() × 5                          [第111-115行]
  │    注册 ext2int / formulation / int2ext / printpf / savecase
  │
  ├─ mpc2gas_prep(mpgc, mpopt)                   [第122行]
  │    │
  │    ├─ compressor2gen(mpgc)                    [mpc2gas_prep.m 第72行]
  │    │    └─ 电驱压缩机 → gen矩阵中的负发电机
  │    │       mpgc.gen = [mpgc.gen; gen_comp]    (PMIN=-inf, PMAX=0)
  │    │       mpgc.genid.comp = [ng+1 : ng+nc_p]
  │    │
  │    ├─ multi_period(mpgc, time, pd, qd)        [mpc2gas_prep.m 第94行]
  │    │    └─ ★ 电网多时段堆叠
  │    │       bus: nb → nb*nt (复制nt份，每份用不同PD/QD)
  │    │       gen: ng → ng*nt (bus编号偏移)
  │    │       branch: nl → nl*nt (bus编号偏移)
  │    │       gencost: 每时段成本 × time(j)小时
  │    │       mpgc.multi_period.time = time
  │    │       mpgc.multi_period.status = 1
  │    │
  │    └─ nsd2gen(mpgc, cost_nsd)                 [mpc2gas_prep.m 第101行]
  │         └─ 所有正负荷 → gen矩阵中的可调度负荷（负发电机）
  │            mpgc.gen = [mpgc.gen; gen_dem]      (PG=-PD, PMAX=0)
  │            mpgc.nsd.N = ndl (可调度负荷总数)
  │            mpgc.bus(:,PD) = 0 (原始PD清零)
  │
  └─ results = runopf(mpgc, mpopt)                [第124行]
       │
       └─ MATPOWER标准OPF流程启动
            │
            ├─ ext2int(mpc)
            │    └─ userfcn_mpng_ext2int()         [mpng.m 第126-292行]
            │         ├─ 数据维度验证
            │         ├─ 天然气数据标幺化 mgc_PU()
            │         └─ 检查储气、压缩机比、Weymouth近似等
            │
            ├─ opf_setup(mpc, mpopt)
            │    ├─ om = opf_model(mpc)
            │    ├─ 添加标准OPF变量:
            │    │    Va(nb*nt个), Vm(nb*nt个), Pg(ng*nt个), Qg(ng*nt个)
            │    │    注意：此时ng已包含原始gen+压缩机gen+可调度负荷gen
            │    ├─ 添加标准AC潮流约束:
            │    │    Pmis/Qmis (nb*nt个等式), Sf/St (nl2*nt个不等式)
            │    ├─ 添加标准发电成本
            │    │
            │    └─ ★ run_userfcn('formulation', om, mpopt)
            │         └─ userfcn_mpng_formulation(om)  [mpng.m 第293-878行]
            │              │
            │              ├─ 添加天然气优化变量 (om.add_var):
            │              │    g(nw), p(nn), ovp(nn), unp(nn),
            │              │    sto_diff(ns), sto_out(ns), sto_in(ns),
            │              │    fgo(no), fgopos(no), fgoneg(no),
            │              │    gamma(ngamma),
            │              │    fgc_g(nc_g), psi_g(nc_g), phi(nc_g),
            │              │    fgc_p(nc_p), psi_p(nc_p)
            │              │
            │              ├─ 添加天然气线性约束 (om.add_lin_constraint):
            │              │    p_pos, p_neg, sto, p_comp1, p_comp2, fgo,
            │              │    Nodes(节点气平衡),
            │              │    power_gencomp, power_gencomp_1,
            │              │    SR(旋转备用), hydro_energy(水电日能量),
            │              │    UCpg/UCqg(机组组合), ramp(爬坡)
            │              │
            │              ├─ 添加天然气非线性约束 (om.add_nln_constraint):
            │              │    f_pipe: Weymouth管道方程   (no个不等式)
            │              │    w_comp_g: 气驱压缩机功率方程 (nc_g个等式)
            │              │    g_comp_g: 气驱压缩机耗气方程 (nc_g个等式)
            │              │    w_comp_p: 电驱压缩机功率方程 (nc_p个等式)
            │              │
            │              └─ 添加天然气成本 (om.add_quad_cost):
            │                   gcost, ovpcost, unpcost, gammacost,
            │                   sto_cost, sto_out_cost, sto_in_cost,
            │                   fgoposcost, fgonegcost,
            │                   fgccost_g, fgccost_p
            │
            ├─ opf_execute(om, mpopt)
            │    └─ nlpopf_solver(om) → om.solve()
            │         └─ nlps_master() → ipopt() MEX
            │              IPOPT一次性求解整个大规模NLP:
            │              变量x = [Va; Vm; Pg; Qg; y; g; p; ovp; unp;
            │                       sto_diff; sto_out; sto_in; fgo;
            │                       fgopos; fgoneg; gamma; ...]
            │
            └─ int2ext(results)
                 └─ userfcn_mpng_int2ext()         [mpng.m 第879-1023行]
                      ├─ 提取天然气变量最优值
                      ├─ 反标幺化 mgc_REAL()
                      ├─ 更新mgc数据结构
                      └─ extract_islands() 将堆叠的多时段电网拆回各时段
```

### 2.2 gen矩阵的最终组成

经过mpc2gas_prep处理后，gen矩阵的行按以下顺序排列：

```
gen矩阵 (ng_total × cg):
┌─────────────────────────────────────────────────────┐
│ 时段1：原始发电机 (ng_ori个)                          │  ← genid.original(1:ng_ori)
│ 时段1：压缩机负发电机 (nc_p个)                        │  ← genid.comp(1:nc_p)
├─────────────────────────────────────────────────────┤
│ 时段2：原始发电机 (ng_ori个)                          │  ← genid.original(ng_ori+1:2*ng_ori)
│ 时段2：压缩机负发电机 (nc_p个)                        │  ← genid.comp(nc_p+1:2*nc_p)
├─────────────────────────────────────────────────────┤
│ ...重复nt次...                                        │
├─────────────────────────────────────────────────────┤
│ 时段nt：原始发电机 + 压缩机                           │
├═════════════════════════════════════════════════════┤
│ 可调度负荷 (ndl个, 来自所有nb*nt个bus中PD>0的)        │  ← genid.dl
└─────────────────────────────────────────────────────┘
ng_total = (ng_ori + nc_p) * nt + ndl
```

因此Pg变量的排列是：
```
Pg = [Pg_t1_gen1, ..., Pg_t1_genN, Pg_t1_comp1, ..., Pg_t1_compM,
      Pg_t2_gen1, ..., Pg_t2_genN, Pg_t2_comp1, ..., Pg_t2_compM,
      ...,
      Pg_nt_gen1, ..., Pg_nt_genN, Pg_nt_comp1, ..., Pg_nt_compM,
      Pg_dl_1, ..., Pg_dl_ndl]
```

---

## 3. 问题1：多时段堆叠机制详解

### 3.1 核心思想

MPNG的多时段实现采用"**电网岛堆叠 + 天然气日均平衡**"策略：

- **电力系统**：将单时段电网复制nt份，每份是一个独立的"岛"（bus编号不同，彼此不相连），各岛有不同的负荷
- **天然气系统**：不复制，保持单份。天然气的时间尺度被处理为**一天的累计平衡**
- **耦合方式**：各时段电网的燃气机组出力（Pg×时段小时数×效率）累计贡献到天然气节点的日平衡方程中

### 3.2 multi_period.m 逐行解析

```matlab
%% 核心参数
nb = size(mpc.bus,1);    % 原始bus数
ng = size(mpc.gen,1);    % 原始gen数(含压缩机负gen)
nl = size(mpc.branch,1); % 原始branch数
nt = length(time);       % 时段数

%% 预分配堆叠后的矩阵
bus    = zeros(nb*nt, cb);       % 总bus数 = nb × nt
gen    = zeros(ng*nt, cg);       % 总gen数 = ng × nt
branch = zeros(nl*nt, cl);      % 总branch数 = nl × nt
gencost= zeros(ngen*nt, cgen);  % 总gencost行 = ngen × nt

%% 逐时段复制
for i = 0:(nt-1)
    j = i+1;

    % Bus：复制整个bus矩阵，替换PD和QD为该时段的值
    bus(((i*nb)+1):j*nb, :) = mpc.bus;
    bus(((i*nb)+1):j*nb, PD) = pd(:,j);  % 各bus在时段j的有功负荷
    bus(((i*nb)+1):j*nb, QD) = qd(:,j);  % 各bus在时段j的无功负荷

    % Gen：复制gen矩阵，GEN_BUS编号偏移 nb*i
    gen(((i*ng)+1):j*ng, :) = mpc.gen;
    gen(((i*ng)+1):j*ng, GEN_BUS) = gb_temp + nb*i;
    % 即时段2的发电机连接到bus (nb+1)~(2*nb)

    % Branch：复制branch矩阵，F_BUS和T_BUS偏移 nb*i
    branch(((i*nl)+1):j*nl, F_BUS) = fb_temp + nb*i;
    branch(((i*nl)+1):j*nl, T_BUS) = tb_temp + nb*i;

    % Gencost：复制并乘以时段持续小时数
    gencost(((i*ngen)+1):j*ngen, :) = mpc.gencost;
    gc_temp = modcost(gencost(...), time(j));  % 成本 × 小时数
    % 即如果原成本为$/MWh，乘以6小时变为$/6h → 目标函数单位=$/day
end

%% 重新编号
bus(:, BUS_I) = 1:nb*nt;  % 连续编号 1, 2, ..., nb*nt

%% 记录多时段信息
mpc_out.multi_period.time = time;
mpc_out.multi_period.status = 1;  % 标记已堆叠
```

### 3.3 为什么"岛堆叠"有效？

MATPOWER的opf_setup在构建Ybus和功率平衡方程时，是基于bus/branch矩阵的。由于不同时段的bus编号不同且branch不跨岛连接，Ybus自然是块对角的：

```
Ybus = diag(Ybus_t1, Ybus_t2, ..., Ybus_nt)
```

功率平衡方程 Pmis(i) = 0 对 i = 1 到 nb×nt，自动变成nt组独立的AC潮流方程。**MATPOWER不需要知道这是"多时段"问题**——它只是看到一个有nb×nt个bus的大电网（其中岛之间不相连）。

### 3.4 天然气变量的排列（不堆叠）

天然气变量在formulation回调中添加，每种变量只有**一份**（不按时段复制）：

```
天然气变量 = [g(nw) | p(nn) | ovp(nn) | unp(nn) |
              sto_diff(ns) | sto_out(ns) | sto_in(ns) |
              fgo(no) | fgopos(no) | fgoneg(no) |
              gamma(ngamma) |
              fgc_g(nc_g) | psi_g(nc_g) | phi(nc_g) |
              fgc_p(nc_p) | psi_p(nc_p)]
```

其中 nw=气井数, nn=气节点数, ns=储气站数, no=管道数, nc_g/nc_p=气驱/电驱压缩机数, ngamma=ngd×ndn(切气变量数)。

### 3.5 时段间耦合约束的添加方式

所有时段间耦合都在`userfcn_mpng_formulation`中以**线性约束**方式添加：

#### (a) 燃气机组的电-气耦合（第691-709行）

这是**电力时段 → 天然气日平衡**的关键耦合。

```matlab
% 构建 A_pggas 矩阵 (nn × ng_total)
% 对于连接到气节点k的燃气机组i，在所有时段j：
%   节点k的气平衡 += Σ_j [-time(j) × eff_pg(i)] × Pg(i,j)
% 即各时段的Pg乘以小时数和效率，累计消耗天然气

val_pggas = -(time' * eff_pg');  % nt × n_gasfired 矩阵
val_pggas = val_pggas(:);        % 拉成列向量
A_pggas = sparse(row_pggas, col_pggas, val_pggas, nn, ng);
```

最终通过节点气平衡约束连接：
```matlab
om.add_lin_constraint('Nodes', A_node, L_node, U_node, {'fgo',...,'Pg','gamma'});
```

**注意**：这里`'Pg'`直接引用了MATPOWER标准OPF的Pg变量（包含所有时段的所有发电机），MPNG通过精确的行列索引将不同时段的Pg映射到对应的气节点。

#### (b) 爬坡约束（第848-876行）

```matlab
% 构建 Pg(t+1) - Pg(t) 的差分约束矩阵
idg_ori = mpc.genid.original;  % 各时段原始发电机在gen矩阵中的行索引
ng_pg_ori = ng_ori/nt;          % 单时段原始发电机数

col1 = idg_ori(ng_pg_ori+1:end);  % Pg(t+1)的列索引
col2 = idg_ori(1:end-ng_pg_ori);  % Pg(t)的列索引
val = [ones(...); -ones(...)];     % Pg(t+1) - Pg(t)

A_r_t = sparse(row, col, val, ng_pg_ori*(nt-1), ng);
L_r_t = -2 * ramp30 .* ramp_time;  % 下坡限制
U_r_t =  2 * ramp30 .* ramp_time;  % 上坡限制

om.add_lin_constraint('ramp', A_r_t, L_r_t, U_r_t, {'Pg'});
```

**关键技巧**：利用`genid.original`记录的各时段发电机在堆叠后gen矩阵中的位置，直接用稀疏矩阵的行列索引构建跨时段差分约束。

#### (c) 日能量上限约束（第811-831行）

```matlab
% Σ_j Pg(i,j) × time(j) ≤ MaxEnergy(i)
A_ener = sparse(row, col, val, ngenr, ng);  % val包含time向量
om.add_lin_constraint('hydro_energy', A_ener, L_ener, U_ener, {'Pg'});
```

#### (d) 机组组合约束（第833-846行）

```matlab
% 当UC(i,t)=0时，强制Pg(i,t)=0
om.add_lin_constraint('UCpg', A_UC, 0, 0, {'Pg'});
om.add_lin_constraint('UCqg', A_UC, 0, 0, {'Qg'});
```

#### (e) 压缩机功耗时段一致约束（第745-773行）

```matlab
% 电驱压缩机在各时段的功耗必须相等：Pg_comp(t+1) = Pg_comp(t)
% 且第一个时段的Pg_comp = psi_p (天然气侧的压缩机功率变量)
om.add_lin_constraint('power_gencomp', Agencomp, 0, 0, {'Pg'});
om.add_lin_constraint('power_gencomp_1', Agencomp_1, 0, 0, {'Pg','psi_p'});
```

#### (f) 储气约束（第445-469行）

储气变量只有一份（日均），约束为：
```
sto_diff = sto_out - sto_in    (线性约束)
sto_min ≤ sto_0 - sto_diff ≤ sto_max  (通过变量界实现)
```

**⚠️ 重要发现**：MPNG的储气是**单步**的（整天一个sto_diff），不是多时段时序耦合。这与你需要的多时段储氢SOC（每个时段一个SOC值，时段间有递推关系）完全不同。你需要自行设计多时段储氢约束。

### 3.6 旋转备用约束（第777-808行）

```matlab
% 分区域分时段的旋转备用：
% Σ_{i∈zone_k} Pg(i,t) ≤ Σ_{i∈zone_k} Pg_max(i) - SR(k,t)
om.add_lin_constraint('SR', A_sr, L_sr, U_sr, {'Pg'});
```

---

## 4. 问题2：与MATPOWER的接口方式详解

### 4.1 接口机制总结

**MPNG完全通过MATPOWER的标准扩展接口实现，未修改MATPOWER的任何内部文件。**

具体使用的接口：

| MATPOWER API | MPNG使用情况 |
|-------------|-------------|
| `add_userfcn()` | 注册5个回调（ext2int, formulation, int2ext, printpf, savecase） |
| `om.add_var()` | 添加17种天然气变量 |
| `om.add_lin_constraint()` | 添加12种线性约束 |
| `om.add_nln_constraint()` | 添加4种非线性约束（管道流、压缩机功率×2、压缩机耗气） |
| `om.add_quad_cost()` | 添加11种线性/二次成本 |
| `om.get_mpc()` | 在formulation回调中获取mpc数据 |
| `runopf()` | 最终调用标准OPF求解 |

**这意味着你完全可以用相同的方式添加POFC约束，不需要修改MATPOWER的任何文件。**

### 4.2 变量添加详情

在`userfcn_mpng_formulation`中（第572-592行），天然气变量的添加：

```matlab
% 气井产量
om.add_var('g', nw, g0, gmin, gmax);

% 节点压力（标幺，实际是π = p²/p_base²）
om.add_var('p', nn, p0, pmin0, []);  % 无上界，通过ovp松弛

% 超压/欠压松弛变量
om.add_var('ovp', nn, ovp0, ovpmin, []);
om.add_var('unp', nn, unp0, unpmin, []);

% 储气相关（3个变量）
om.add_var('sto_diff', ns, sto_diff0, sto_diffmin, sto_diffmax);
om.add_var('sto_out', ns, sto_out0, sto_outmin, sto_outmax);
om.add_var('sto_in', ns, sto_in0, sto_inmin, sto_inmax);

% 管道流量（分正负方向建模）
om.add_var('fgo', no, fgo0, fgomin, fgomax);
om.add_var('fgopos', no, fgopos0, fgoposmin, fgoposmax);
om.add_var('fgoneg', no, fgoneg0, fgonegmin, fgonegmax);

% 切气变量
om.add_var('gamma', ngamma, gamma0, gammamin, gammamax);

% 气驱压缩机
om.add_var('fgc_g', nc_g, fgc0_g, fgcmin_g, fgcmax_g);
om.add_var('psi_g', nc_g, psi0_g, psimin_g, []);
om.add_var('phi', nc_g, phi0, phimin, []);

% 电驱压缩机
om.add_var('fgc_p', nc_p, fgc0_p, fgcmin_p, fgcmax_p);
om.add_var('psi_p', nc_p, psi0_p, psimin_p, []);
```

### 4.3 Jacobian和Hessian的处理方式

**MPNG不需要手动将天然气的Jacobian/Hessian拼接到电力系统的矩阵中。**

opt_model的工作原理：
1. 每个`add_nln_constraint`调用都指定了`VARSETS`——该约束涉及哪些变量
2. 回调函数只需返回**关于VARSETS中变量的局部Jacobian/Hessian**
3. opt_model内部维护了`变量名 → 全局索引`的映射表（通过get_idx获取）
4. 在组装传给IPOPT的全局Jacobian/Hessian时，opt_model自动将局部矩阵放到正确的位置

**具体例子**：`fpipe_fcn`（管道Weymouth方程约束）

```matlab
% 注册时指定 VARSETS = {'fgo', 'p'}
om.add_nln_constraint('f_pipe', no, 'true', fcn_fpipe, hess_fpipe, {'fgo','p'});
```

```matlab
function [g, dg] = fpipe_fcn(x, parpipe)
    [fgo, p] = deal(x{:});  % x是cell: x{1}=fgo, x{2}=p

    % g(no×1) = fgo - Weymouth(p_from, p_to)
    g = fgo - f;

    % dg(no × (no+nn)) = [∂g/∂fgo | ∂g/∂p]
    jac1 = speye(no);                    % no×no: ∂g/∂fgo = I
    jac2 = sparse(..., npipe, nn);       % no×nn: ∂g/∂p
    dg = [jac1 jac2];
end
```

opt_model内部会将这个 no×(no+nn) 的局部Jacobian 映射到全局 Jacobian 矩阵中 fgo 和 p 对应的列位置。

**Hessian同理**：
```matlab
function d2L = fpipe_hess(x, lambda, parpipe)
    % 返回 (no+nn) × (no+nn) 的局部Hessian
    d2L = [hess_fgo_fgo    hess_fgo_p;     % no×no  no×nn
           hess_p_fgo      hess_p_p    ];   % nn×no  nn×nn
end
```

opt_model将这个矩阵嵌入全局Hessian中对应fgo和p的行列位置。

### 4.4 线性约束引用MATPOWER标准变量

MPNG在添加线性约束时直接引用MATPOWER标准变量名：

```matlab
% 节点气平衡——同时涉及天然气变量'fgo','g','sto_diff'和电力变量'Pg'
om.add_lin_constraint('Nodes', A_node, L_node, U_node,
    {'fgo', 'fgc_p', 'fgc_g', 'phi', 'g', 'sto_diff', 'Pg', 'gamma'});

% 爬坡约束——仅涉及电力变量'Pg'
om.add_lin_constraint('ramp', A_r_t, L_r_t, U_r_t, {'Pg'});

% 旋转备用——仅涉及'Pg'
om.add_lin_constraint('SR', A_sr, L_sr, U_sr, {'Pg'});
```

**这证明了：你可以在userfcn中自由混合引用MATPOWER的标准变量（Va/Vm/Pg/Qg）和你自己添加的变量。**

---

## 5. 问题3：目标函数的扩展方式详解

### 5.1 目标函数的组成

MPNG的总目标函数 = MATPOWER标准发电成本 + 天然气相关成本

**MATPOWER标准发电成本**：在opf_setup中根据gencost矩阵自动添加（通过add_quad_cost或add_nln_cost），用户无需处理。注意gencost已经被multi_period中的modcost乘以了时段小时数。

**天然气成本**：在formulation回调中添加（第627-641行）。

### 5.2 天然气成本详解

所有天然气成本都使用`om.add_quad_cost(name, Q, c, k, varsets)`添加，格式为`f = ½x'Qx + c'x + k`。

由于MPNG的天然气成本全部是线性的（Q=[]），实际形式为 `f = c'x + k`：

```matlab
% (1) 气井开采成本：gcost(i) × g(i)，单位 $/day
om.add_quad_cost('gcost', [], gcost, 0, {'g'});
% gcost已乘以fbase标幺化：gcost_pu = gcost_real * fbase

% (2) 超压惩罚：ovpcost(i) × ovp(i)
om.add_quad_cost('ovpcost', [], ovpcost, 0, {'ovp'});

% (3) 欠压惩罚
om.add_quad_cost('unpcost', [], unpcost, 0, {'unp'});

% (4) 切气惩罚：gammacost(k) × gamma(k)
om.add_quad_cost('gammacost', [], gammacost, 0, {'gamma'});

% (5) 储气成本（含初始储气的常数项）
% f = -sto_cost' × sto_diff + sto_cost' × sto_0
% 即最终储量越高（sto_diff越小），成本越低
om.add_quad_cost('sto_cost', [], -sto_cost, k_sto_cost, {'sto_diff'});
% k_sto_cost = sto_cost .* sto0

% (6) 储气流出/流入成本
om.add_quad_cost('sto_out_cost', [], sto_out_cost, 0, {'sto_out'});
om.add_quad_cost('sto_in_cost', [], -sto_in_cost, 0, {'sto_in'});
% 注意sto_in的系数为负——表示注入储气是一种"收益"

% (7) 管道运输成本（分正负方向）
om.add_quad_cost('fgoposcost', [], fgocost, 0, {'fgopos'});
om.add_quad_cost('fgonegcost', [], -fgocost, 0, {'fgoneg'});

% (8) 压缩机运行成本
om.add_quad_cost('fgccost_g', [], fgccost_g, 0, {'fgc_g'});
om.add_quad_cost('fgccost_p', [], fgccost_p, 0, {'fgc_p'});
```

### 5.3 对你的参考价值

你可以用完全相同的方式添加：
- **碳交易成本**：`om.add_quad_cost('carbon', [], c_carbon .* emission_factor, 0, {'Pg'})`
- **氢销售收益**（负成本）：`om.add_quad_cost('H2_rev', [], -H2_price, 0, {'H2'})`
- 如果碳成本或氢收益是非线性的，改用`om.add_nln_cost()`

---

## 6. 问题4：数据结构详解

### 6.1 顶层结构

MPNG的输入数据是一个扩展的MATPOWER mpc结构体：

```
mpgc (顶层结构体)
│
├── .baseMVA          标准MATPOWER字段
├── .bus              标准MATPOWER bus矩阵 [nb × 13+]
├── .gen              标准MATPOWER gen矩阵 [ng × 25]
├── .branch           标准MATPOWER branch矩阵 [nl × 21]
├── .gencost          标准MATPOWER gencost矩阵 [ng × 7+]
│
├── .mgc              ★ 天然气系统数据（MPNG新增）
│   ├── .version      版本号 '1'
│   ├── .fbase        流量基准 (MMSCFD)，如50
│   ├── .pbase        压力基准 (PSI)，如500
│   ├── .wbase        功率基准 (MW)，如100（建议与baseMVA一致）
│   │
│   ├── .node
│   │   ├── .info     [nn × (11+NGD)] 节点信息矩阵
│   │   │              列: NODE_I, NODE_TYPE, PR, PRMAX, PRMIN,
│   │   │                   OVP, UNP, COST_OVP, COST_UNP, GD, NGD
│   │   ├── .dem      [nn × ngd] 各类型气需求矩阵
│   │   └── .demcost  [nn × ngd] 各类型切气成本矩阵
│   │
│   ├── .well          [nw × 7] 气井数据矩阵
│   │                  列: WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G
│   │
│   ├── .pipe          [no × 9] 管道数据矩阵
│   │                  列: F_NODE, T_NODE, FG_O, K_O, DIAM, LNG,
│   │                       FMAX_O, FMIN_O, COST_O
│   │
│   ├── .comp          [nc × 16] 压缩机数据矩阵
│   │                  列: COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C,
│   │                       RATIO_C, B_C, Z_C, ALPHA, BETA, GAMMA,
│   │                       FMAX_C, COST_C, FG_C, GC_C, PC_C
│   │
│   ├── .sto           [ns × 14] 储气数据矩阵
│   │                  列: STO_NODE, STO, STO_0, STOMAX, STOMIN,
│   │                       FSTO, FSTO_OUT, FSTO_IN, FSTOMAX, FSTOMIN,
│   │                       S_STATUS, COST_STO, COST_OUT, COST_IN
│   │
│   └── .wey_percent   (可选) Weymouth近似区间百分比，标量或no×1向量
│
└── .connect           ★ 电-气互联数据（MPNG新增）
    ├── .version       版本号 '1'
    │
    ├── .power         电力系统时序数据
    │   ├── .time      [1 × nt] 各时段持续小时数，如 [6 6 6 6]，sum=24
    │   ├── .demands
    │   │   ├── .pd    [nb × nt] 各bus各时段有功负荷 (MW)
    │   │   └── .qd    [nb × nt] 各bus各时段无功负荷 (MVar)
    │   ├── .cost      切负荷成本 ($/MWh)，标量
    │   ├── .sr        [nareas × nt] 分区域分时段旋转备用 (MW)，可为[]
    │   ├── .energy    [ngenr × 2] 日最大发电能量约束
    │   │               列: GEN_ID, MAX_ENER (MWh/day)
    │   ├── .UC        [ng_ori × nt] 机组组合矩阵，0/1，可为[]
    │   └── .ramp_time [ng_ori × nt] 爬坡时间矩阵 (小时)，可为[]
    │
    └── .interc        互联拓扑
        ├── .comp      [nc_p × 2] 电驱压缩机-电网连接
        │               列: COMP_ID, BUS_ID
        └── .term      [n_gasfired × 3] 燃气机组-气网连接
                        列: GEN_ID, NODE_ID, EFF (MMSCFD/MW)
```

### 6.2 数据在处理流程中的变化

经过mpc2gas_prep后，mpgc结构体新增以下字段：

```
mpgc（经mpc2gas_prep处理后）
│
├── .genid             发电机分类索引（★关键）
│   ├── .original      [ng_ori×nt × 1] 原始发电机在堆叠后gen矩阵中的行索引
│   ├── .comp          [nc_p×nt × 1] 压缩机负发电机的行索引
│   └── .dl            [ndl × 1] 可调度负荷的行索引
│
├── .gencomp           压缩机负发电机信息
│   ├── .N             电驱压缩机数量
│   └── .id            在gen矩阵中的行索引
│
├── .nsd               可调度负荷信息
│   ├── .N             可调度负荷数量
│   ├── .id_dem        在bus矩阵中原始PD>0的bus索引
│   └── .original
│       ├── .PD        原始有功负荷向量
│       └── .QD        原始无功负荷向量
│
├── .multi_period      多时段信息
│   ├── .time          各时段持续小时数
│   └── .status        1 = 已堆叠
│
└── bus(:,PD) = 0      ★ 原始PD/QD已清零（负荷转为可调度负荷）
```

### 6.3 数据与mpc的关联方式

MPNG将所有自定义数据作为mpc结构体的**附加字段**：
- `mpc.mgc` → 天然气系统
- `mpc.connect` → 互联数据
- `mpc.genid` → 发电机分类索引
- `mpc.nsd` → 可调度负荷信息

MATPOWER会忽略它不认识的字段，在ext2int/int2ext转换时这些字段会被保留。在formulation回调中通过`om.get_mpc()`即可访问。

---

## 7. 问题5：求解器调用方式详解

### 7.1 调用链

```
mpng.m 第124行:
  results = runopf(mpgc, mpopt)
    → opf(mpgc, mpopt)                    [lib/opf.m]
      → opf_setup(mpc, mpopt)             [lib/opf_setup.m]
      → opf_execute(om, mpopt)            [lib/opf_execute.m]
        → nlpopf_solver(om, mpopt)        [lib/nlpopf_solver.m]
          → om.solve(opt)                 [mp-opt-model/@opt_model/solve.m]
            → nlps_master(...)            [mp-opt-model/nlps_master.m]
              → nlps_ipopt(...)           IPOPT MEX接口
                → ipopt(prob)             底层C++求解器
```

### 7.2 MPNG设置的默认IPOPT选项

```matlab
mpopt.opf.ac.solver = 'IPOPT';     % 使用IPOPT求解AC OPF
mpopt.ipopt.opts.max_iter = 1e5;   % 最大迭代次数100000
```

### 7.3 传给IPOPT的NLP问题格式

```
min   f(x) = Σ(发电成本) + Σ(天然气成本)
 x
s.t.
  [标准AC潮流等式约束]
    Pmis(x) = 0     (nb×nt个)
    Qmis(x) = 0     (nb×nt个)

  [标准支路潮流不等式约束]
    Sf(x) ≤ 0       (nl2×nt个)
    St(x) ≤ 0       (nl2×nt个)

  [天然气非线性约束]
    f_pipe(fgo,p) ≤ 0  (no个, Weymouth管道方程)
    w_comp_g(...) = 0   (nc_g个, 气驱压缩机功率)
    g_comp_g(...) = 0   (nc_g个, 气驱压缩机耗气)
    w_comp_p(...) = 0   (nc_p个, 电驱压缩机功率)

  [线性约束]
    所有add_lin_constraint添加的约束
    （节点气平衡、压力限制、储气、爬坡、旋转备用、能量、UC等）

  [变量界]
    xmin ≤ x ≤ xmax
```

变量x的完整组成（按添加顺序）：
```
x = [Va(nb×nt) | Vm(nb×nt) | Pg(ng_total×nt+ndl) | Qg(ng_total×nt+ndl) |
     y(ny) |
     g(nw) | p(nn) | ovp(nn) | unp(nn) |
     sto_diff(ns) | sto_out(ns) | sto_in(ns) |
     fgo(no) | fgopos(no) | fgoneg(no) | gamma(ngamma) |
     fgc_g(nc_g) | psi_g(nc_g) | phi(nc_g) |
     fgc_p(nc_p) | psi_p(nc_p)]
```

### 7.4 问题规模估算

以Example4（118-bus + 48-node, 24时段）为例：

| 量 | 估算 |
|----|------|
| 电力变量 | Va: 118×24=2832, Vm: 2832, Pg: ~54×24+ndl≈1500, Qg: ~1500 |
| 天然气变量 | g+p+ovp+unp+sto+fgo+...≈ 300-500 |
| **总变量** | **~9000-10000** |
| AC潮流约束 | 2×118×24 = 5664等式 |
| 支路约束 | ~186×24×2 = ~8928不等式 |
| 天然气约束 | ~200-300 |
| 线性约束 | ~2000-3000 |
| **总约束** | **~17000-18000** |

这是IPOPT可以轻松处理的规模。

---

## 8. 路径评估："在MPNG上修改" vs "在MATPOWER上重写"

### 8.1 路径A：在MPNG代码上修改

#### 需要做的工作

| 步骤 | 工作量 | 说明 |
|------|--------|------|
| 删除天然气网络代码 | 中 | 删除管道、压缩机、气井、节点气平衡、Weymouth等（约800行） |
| 删除compressor2gen | 小 | 你没有压缩机 |
| 修改multi_period | 小 | 基本可复用，可能需要微调 |
| 保留nsd2gen | 可复用 | 切负荷建模仍然需要 |
| 重写formulation回调 | 大 | 完全替换为POFC变量/约束/成本 |
| 重写ext2int回调 | 中 | 替换数据验证逻辑 |
| 重写int2ext回调 | 中 | 替换结果提取逻辑 |
| 重写printpf回调 | 中 | 输出格式完全不同 |
| 添加储能多时段SOC约束 | 大 | MPNG没有这个，需从零设计 |
| 添加储氢多时段约束 | 大 | 需从零设计 |
| 添加POFC可行域约束 | 大 | 核心工作，需编写fcn+Jacobian+Hessian |

**总评**：需要删除和重写的代码量 > 可复用的代码量

#### 风险
- mpng.m是单一大文件（1646行），修改过程容易出错
- 残留的天然气相关逻辑可能导致隐蔽bug
- 代码结构不够清晰，维护困难

### 8.2 路径B：参考MPNG，在MATPOWER上重新搭建

#### 需要做的工作

| 步骤 | 工作量 | 说明 |
|------|--------|------|
| 搭建主框架 | 小 | 参考mpng.m前125行，几乎可以照搬模式 |
| 复用multi_period | 可直接用 | 直接复制或微调 |
| 复用nsd2gen | 可直接用 | 同上 |
| 编写ext2int回调 | 中 | 数据验证+标幺化（如果需要） |
| 编写formulation回调 | 大 | ★核心工作 |
| 编写int2ext回调 | 中 | 结果提取 |
| 编写POFC可行域约束 | 大 | fcn + Jacobian + Hessian |
| 编写储能SOC约束 | 中 | 线性约束即可 |
| 编写储氢约束 | 中 | 线性约束即可 |
| 编写碳交易+氢收益成本 | 小 | add_quad_cost即可 |

**总评**：代码更清晰，工作量与路径A接近但可维护性更好

### 8.3 综合评估

| 维度 | 路径A（修改MPNG） | 路径B（重新搭建） |
|------|-------------------|-------------------|
| 初始工作量 | ★★★（需要先理解再删改） | ★★★（需要搭建框架） |
| 代码清晰度 | ★★（单文件残留混合逻辑） | ★★★★★（模块化，专为你的问题设计） |
| 调试难度 | ★★（可能有隐蔽的残留问题） | ★★★★（干净起步，问题定位明确） |
| 可维护性 | ★★ | ★★★★★ |
| 可扩展性 | ★★ | ★★★★★ |
| 论文展示 | ★★★ | ★★★★（代码结构更专业） |

### 8.4 最终建议

**推荐路径B**，具体策略为"**复用核心模式，重写业务逻辑**"：

**直接复用**：
- `multi_period.m`（电网岛堆叠机制，99%可直接用）
- `nsd2gen.m`（切负荷建模，100%可直接用）
- userfcn注册框架模式（mpng.m的第77-124行可作为模板）
- 爬坡约束的稀疏矩阵构建方式（mpng.m第848-876行可作为参考）

**参考思路但需重写**：
- 非线性约束的Jacobian/Hessian编写模式（fpipe_fcn/hess是很好的模板）
- 跨时段线性约束的索引构建技巧
- 标幺化策略

**完全新写**：
- POFC可行域约束函数
- 多时段储能SOC递推约束
- 多时段储氢SOC递推约束
- 碳交易和氢收益的成本函数
- 新能源出力上限约束

---

## 9. 对POFC项目的具体参考价值与局限性分析

### 9.1 高参考价值的部分

#### (a) userfcn扩展框架（★★★★★）
MPNG证明了**不修改MATPOWER内部代码就能实现复杂的多能系统耦合OPF**的可行性。你可以完全沿用这个框架。

#### (b) 岛堆叠多时段机制（★★★★★）
multi_period.m是经过验证的多时段实现方式，对于你的24时段问题直接适用。

#### (c) 非线性约束的Jacobian/Hessian编写范式（★★★★）
`fpipe_fcn`/`fpipe_hess`、`wcompgas_fcn`/`wcompgas_hess`提供了完整的编写模板：
- 如何从cell array中解包变量
- 如何构建稀疏Jacobian
- 如何构建Lagrangian加权Hessian
- 如何处理分块矩阵结构

#### (d) 跨时段约束的稀疏矩阵索引技巧（★★★★）
爬坡约束、能量约束、UC约束的实现展示了如何利用genid索引在堆叠后的gen矩阵中精确定位各时段的发电机，构建跨时段差分矩阵。

#### (e) 切负荷建模（★★★）
nsd2gen的"负发电机"方法是成熟的切负荷建模方式，可直接复用。

### 9.2 局限性（需要你自行解决的部分）

#### (a) 天然气无多时段时序耦合（★关键局限★）
MPNG的天然气网络是**日均平衡**（单步），储气也只有一个sto_diff变量。但你的：
- **储能（电池）**：需要每个时段一个SOC值，`SOC(t+1) = SOC(t) ± P_batt(t)*dt`
- **储氢**：需要每个时段一个储氢量，`H2_tank(t+1) = H2_tank(t) + H2_prod(t)*dt - H2_demand(t)*dt`

这些都需要你自行设计**多时段时序递推约束**。建议方式：
```matlab
% 储能SOC变量：nt个时段
om.add_var('SOC', nt, SOC0, SOC_min*ones(nt,1), SOC_max*ones(nt,1));
om.add_var('P_charge', nt, zeros(nt,1), zeros(nt,1), P_ch_max*ones(nt,1));
om.add_var('P_discharge', nt, zeros(nt,1), zeros(nt,1), P_dis_max*ones(nt,1));

% SOC递推约束（线性等式）：
% SOC(t+1) - SOC(t) + P_dis(t)*dt/E_cap - P_ch(t)*dt*η/E_cap = 0
% 构建 (nt-1) × (3*nt) 的稀疏矩阵A
om.add_lin_constraint('SOC_dynamics', A_soc, zeros(nt-1,1), zeros(nt-1,1),
    {'SOC', 'P_charge', 'P_discharge'});
```

#### (b) POFC可行域的非线性程度
MPNG的非线性约束（Weymouth方程、压缩机功率方程）是解析函数，导数可以精确计算。你的POFC可行域如果是：
- 解析函数（如椭圆/多项式）→ 可以直接写Jacobian/Hessian
- 查表数据的拟合 → 需要选择合适的拟合函数并推导导数
- 黑箱模型 → 需要用有限差分近似导数（但IPOPT效率会降低）

#### (c) 电-氢耦合的非线性特性
MPNG的电-气耦合是线性的（Pg × 效率 = 气消耗）。你的POFC电-氢耦合如果是非线性的（如可行域边界是曲线），需要作为非线性约束添加。

#### (d) 新能源出力处理
MPNG没有显式处理可再生能源（风/光），你需要自行添加。最简单的方式是将风/光作为特殊的发电机（PMIN=0, PMAX=forecast_value(t)），每个时段的PMAX不同——这可以在multi_period中通过修改gen矩阵实现，或在formulation中通过约束实现。

### 9.3 建议的代码组织结构

```
pofc_opf/
├── pofc_opf.m                   % 主入口（~50行，参考mpng.m前125行）
├── multi_period.m               % 直接复用MPNG的
├── pofc_prep.m                  % 数据预处理（参考mpc2gas_prep.m）
│                                  调用 multi_period → nsd2gen → 新能源处理
├── nsd2gen.m                    % 直接复用MPNG的
│
├── userfcn_pofc_ext2int.m       % 数据验证（~100行）
├── userfcn_pofc_formulation.m   % ★核心：添加变量/约束/成本（~300-500行）
├── userfcn_pofc_int2ext.m       % 结果提取（~100行）
├── userfcn_pofc_printpf.m       % 结果打印（可选）
│
├── constraints/
│   ├── pofc_feasibility_fcn.m   % POFC可行域约束值+Jacobian
│   ├── pofc_feasibility_hess.m  % POFC可行域约束Hessian
│   └── (如有更多非线性约束...)
│
├── idx_pofc.m                   % POFC数据列索引常量
├── idx_hydrogen.m               % 氢相关数据列索引常量
│
├── cases/
│   ├── case_XXbus.m             % 电力系统算例
│   ├── pofc_params.m            % POFC电厂参数（可行域系数等）
│   ├── renewable_data.m         % 风/光出力时序数据
│   └── connect_pofc.m           % 互联和时序数据
│
└── utils/
    ├── plot_results.m           % 结果可视化
    └── analyze_results.m        % 经济性分析
```

### 9.4 formulation回调的骨架设计

```matlab
function om = userfcn_pofc_formulation(om, mpopt, args)
    mpc = om.get_mpc();
    baseMVA = mpc.baseMVA;
    pofc = mpc.pofc;
    nt = length(mpc.connect.time);

    %% ===== 1. POFC电厂变量 =====
    % 每个时段一个H2产量
    om.add_var('H2', nt, H2_0, H2_min, H2_max);
    % POFC净电出力已包含在Pg中（作为特殊发电机）

    %% ===== 2. 储能变量 =====
    om.add_var('SOC', nt, SOC_0, SOC_min, SOC_max);
    om.add_var('P_ch', nt, zeros, zeros, P_ch_max);
    om.add_var('P_dis', nt, zeros, zeros, P_dis_max);

    %% ===== 3. 储氢变量 =====
    om.add_var('H2_tank', nt, H2_tank_0, H2_tank_min, H2_tank_max);

    %% ===== 4. POFC可行域非线性约束（每个时段）=====
    fcn = @(x) pofc_feasibility_fcn(x, pofc_params);
    hess = @(x, lam) pofc_feasibility_hess(x, lam, pofc_params);
    om.add_nln_constraint('pofc_region', n_con, 0, fcn, hess, {'Pg', 'H2'});

    %% ===== 5. 储能SOC递推约束（线性）=====
    A_soc = build_soc_matrix(nt, dt, eta_ch, eta_dis, E_cap);
    om.add_lin_constraint('SOC_dyn', A_soc, zeros, zeros,
        {'SOC', 'P_ch', 'P_dis'});

    %% ===== 6. 储氢递推约束（线性）=====
    A_h2tank = build_h2tank_matrix(nt, dt);
    om.add_lin_constraint('H2tank_dyn', A_h2tank, zeros, zeros,
        {'H2_tank', 'H2'});

    %% ===== 7. 新能源出力上限约束 =====
    % 通过修改gen的PMAX实现，或添加线性约束

    %% ===== 8. 爬坡约束 =====
    % 参考MPNG的ramp约束构建方式

    %% ===== 9. 成本 =====
    om.add_quad_cost('carbon', [], c_carbon, 0, {'Pg'});
    om.add_quad_cost('H2_rev', [], -price_H2, 0, {'H2'});
end
```

---

*报告完成。以上所有分析基于对mpng-master和matpower7.1源码的实际阅读，代码行号均指向实际文件位置。*
