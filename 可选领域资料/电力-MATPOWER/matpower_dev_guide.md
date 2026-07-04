# MATPOWER opt_model Extension Development Guide

> Practical lessons from building a multi-period AC OPF extension using MATPOWER 7.1 + IPOPT.
> This guide covers the undocumented patterns, API usage, and common pitfalls.

---

## 1. Architecture: userfcn Callback Mechanism

MATPOWER's standard way to extend OPF is via **user function callbacks**. Three hooks
are available, called in order during `runopf()`:

```
your_main_function(mpc, mpopt)
  |
  +-- Preprocessing: validate & transform mpc
  +-- add_userfcn(mpc, 'ext2int',     @your_ext2int_cb)
  +-- add_userfcn(mpc, 'formulation', @your_formulation_cb, args)
  +-- add_userfcn(mpc, 'int2ext',     @your_int2ext_cb)
  |
  +-- runopf(mpc, mpopt)    <-- MATPOWER takes over
        |
        +-- opf()
              +-- ext2int(mpc)         --> calls your_ext2int_cb
              +-- opf_setup(mpc)       --> builds om (opt_model)
              |     +-- adds standard vars: Va, Vm, Pg, Qg
              |     +-- run_userfcn('formulation', om)
              |           --> calls your_formulation_cb(om, mpopt, args)
              |                 om.add_var(...)
              |                 om.add_lin_constraint(...)
              |                 om.add_nln_constraint(...)
              |                 om.add_quad_cost(...)
              +-- opf_execute(om)
              |     +-- nlpopf_solver  --> om.solve() --> IPOPT
              |     +-- populates results.var.val.*, results.lin.mu.*, etc.
              +-- int2ext(results)     --> calls your_int2ext_cb
                    +-- extract custom results from results.var.val.*
```

**Key points:**
- `ext2int`: data transformation (external → internal ordering). Runs before formulation.
- `formulation`: add variables, constraints, costs to `om` (opt_model). This is where all the action happens.
- `int2ext`: extract results after solving. Convert back to external ordering.

---

## 2. Core opt_model API

### 2.1 add_var — Add optimization variables

```matlab
om.add_var(NAME, N, V0, VL, VU);
% NAME : string identifier (e.g., 'my_var')
% N    : number of variables (scalar)
% V0   : initial value vector [N x 1]
% VL   : lower bound vector [N x 1]
% VU   : upper bound vector [N x 1]
```

### 2.2 add_lin_constraint — Linear constraints

```matlab
om.add_lin_constraint(NAME, A, L, U, VARSETS);
% Enforces: L <= A * x <= U
% A       : [num_constraints x num_variables_in_VARSETS] sparse matrix
% L, U    : bound vectors [num_constraints x 1]
% VARSETS : cell array of variable names, e.g., {'Pg'} or {'var1', 'var2'}
```

**Multi-VARSETS column layout:** When using `{'var1', 'var2'}`, the A matrix columns are
`[var1(1..n1), var2(1..n2)]` concatenated:
- var1 columns: 1 to n1
- var2 columns: n1+1 to n1+n2

### 2.3 add_nln_constraint — Nonlinear constraints

```matlab
om.add_nln_constraint(NAME, N, ISEQ, FCN, HESS, VARSETS);
% N       : number of constraints
% ISEQ    : 1 for equality (g=0), 0 for inequality (g<=0)
% FCN     : @(x) -> [g, dg]
% HESS    : @(x, lam) -> H
% VARSETS : cell array of variable names
```

**FCN details:**
- When VARSETS is used, `x` is a **cell array**: `x{1}`, `x{2}`, ...
- `g`: [N x 1] constraint values
- `dg`: [N x M] Jacobian matrix (**rows = constraints, cols = variables, NOT transposed**)

**HESS details:**
- `H`: [M x M] weighted Hessian = Σᵢ lam(i) × ∂²gᵢ/∂x²

**ISEQ note:** Both integer (1/0) and string ('true'/'=') work in MATPOWER 7.1. Recommend using integer 0/1 to match MATPOWER's internal convention.

### 2.4 add_quad_cost — Quadratic/linear cost terms

```matlab
om.add_quad_cost(NAME, Q, C, K, VARSETS);
% cost = 0.5 * x' * Q * x + C' * x + K
% Q : [N x N] sparse matrix (use sparse(N,N) or [] for zero)
% C : [N x 1] vector
% K : scalar constant
% VARSETS : cell array of variable names
```

**Tips:**
- Q=0 → pure linear cost: C'*x
- Negative C pushes variable upward (revenue/benefit terms)
- Custom costs via `add_quad_cost` are NOT auto-scaled by period duration; you must include Δt manually

---

## 3. Critical Technical Notes

### 3.1 Per-Unit Conversion (MW ↔ p.u.)

**Inside formulation**, all power variables (Pg, Qg) are in **per-unit** (÷ baseMVA).

```matlab
baseMVA = om.get_mpc().baseMVA;
ramp_pu = ramp_mw / baseMVA;        % convert MW to p.u. for constraints
```

After solving, `nlpopf_solver.m` converts back: `gen(:, PG) = Pg * baseMVA`.

When building constraint A matrices that mix custom variables (in physical units) with Pg (in p.u.), you must account for baseMVA. Example:

```matlab
% Constraint: Pg_pofc * baseMVA = a1 * M_coal - a2 * G_gas
% Rearranged: a1*M - a2*G - baseMVA*Pg = 0
% VARSETS = {'M_coal', 'G_gas', 'Pg'}
% A columns: [M_coal cols | G_gas cols | Pg cols]
A(t, col_M(t))  =  a1;
A(t, col_G(t))  = -a2;
A(t, col_Pg(t)) = -baseMVA;
```

### 3.2 Jacobian Format for NLN Constraints

MATPOWER expects `dg` as **N × M** (not transposed):
- N = number of constraints
- M = total number of variables in VARSETS
- Row i = gradient of constraint i w.r.t. all variables

Confirmed from `opf_power_balance_fcn.m` (lines 83-86 in MATPOWER 7.1).

### 3.3 gencost Scaling

`modcost(gencost, time)` scales the entire cost curve by `time` (hours):
- Polynomial `c2*P² + c1*P + c0` becomes `c2*T*P² + c1*T*P + c0*T`
- Converts $/h rate → $/period total cost
- Must be done **before** stacking for multi-period

### 3.4 Variable Extraction from Results

After solving, `opf_execute.m` (lines ~211-221) automatically populates:

```matlab
results.var.val.<name>      % variable values by block name
results.var.mu.l.<name>     % lower bound shadow prices
results.var.mu.u.<name>     % upper bound shadow prices
results.lin.mu.l.<name>     % linear constraint lower multipliers
results.lin.mu.u.<name>     % linear constraint upper multipliers
results.nli.mu.<name>       % NLN inequality multipliers
results.qdc.<name>          % quadratic cost values
```

Example: `results.var.val.my_var` gives the optimal values of your custom variable.

### 3.5 Useful opt_model Methods

```matlab
mpc = om.get_mpc();                      % retrieve the case struct from om
n   = om.getN('var', 'Pg');              % get size of variable block 'Pg'
[v0, vl, vu] = om.params_var('Pg');      % get initial/bounds of variable
```

---

## 4. Multi-Period OPF via Island Stacking

### 4.1 Concept

Stack the entire network nt times as electrically disconnected "islands":
- Bus numbers offset by `(t-1) * nb` for period t
- gen, branch tables duplicated with updated bus references
- gencost scaled by period duration via `modcost()`
- Load (Pd, Qd) scaled by per-period factors

```matlab
% Pseudo-code for multi-period stacking
for t = 2:nt
    new_bus = bus_orig;
    new_bus(:, BUS_I) = new_bus(:, BUS_I) + (t-1) * nb;
    % ... scale loads by pd_factor(t) ...
    mpc.bus = [mpc.bus; new_bus];

    new_gen = gen_orig;
    new_gen(:, GEN_BUS) = new_gen(:, GEN_BUS) + (t-1) * nb;
    mpc.gen = [mpc.gen; new_gen];

    new_branch = branch_orig;
    new_branch(:, F_BUS) = new_branch(:, F_BUS) + (t-1) * nb;
    new_branch(:, T_BUS) = new_branch(:, T_BUS) + (t-1) * nb;
    mpc.branch = [mpc.branch; new_branch];

    new_gencost = modcost(gencost_orig, time(t));
    mpc.gencost = [mpc.gencost; new_gencost];
end
```

### 4.2 Important Notes

- Works because MATPOWER treats each island independently for power flow
- Bus numbers must be sequential (1..nb) in each period for simple offset logic
- After `ext2int()`, out-of-service generators may be removed — always use `om.getN('var', 'Pg') / nt` for the internal generator count, NOT the original count
- Similarly, in `int2ext`, use `size(results.gen, 1) / nt` for gen count

### 4.3 Cross-Period Constraints

MATPOWER stacks all period copies of Pg into one block: `[Pg_t1; Pg_t2; ...; Pg_tnt]`.

Ramp constraint example (linking consecutive periods):
```
A = [ -I_ng  +I_ng   0      0    ]   % period 1→2
    [  0     -I_ng  +I_ng   0    ]   % period 2→3
    [  0      0     -I_ng  +I_ng ]   % period 3→4

VARSETS = {'Pg'}
L = -ramp_pu * ones(ng*(nt-1), 1);
U =  ramp_pu * ones(ng*(nt-1), 1);
```

For custom variables (e.g., SOC), they are already sized [nt × 1], so one A matrix
expresses all cross-period relationships directly.

---

## 5. Design Patterns

### 5.1 Virtual Generator Pattern

Model non-generation devices (storage, fuel cell, etc.) as virtual generators:
- **Discharge**: gen with `0 ≤ Pg ≤ P_dis_max`, gencost ≈ 0
- **Charge**: gen with `-P_ch_max ≤ Pg ≤ 0`, gencost ≈ 0
- Add virtual gens to `mpc.gen` **before** multi-period stacking
- Link to custom variables via `add_lin_constraint` in formulation (e.g., SOC dynamics)

**Pitfall — Dispatchable load Q limits:** MATPOWER treats gens with PMIN < 0 as
"dispatchable loads". The `makeAvl` function **requires** either QMIN == 0 or QMAX == 0
for these. For charge gens (PMIN < 0, PMAX = 0), set `QMIN = QMAX = 0`.

### 5.2 SOC Recursion Pattern (Battery/Storage)

```matlab
% SOC(t) = SOC(t-1) + eta_ch * P_ch(t) * dt - P_dis(t)/eta_dis * dt
% Express as: A_soc * [SOC; Pg] = b_soc (via VARSETS = {'SOC', 'Pg'})
%
% A matrix structure (for nt periods):
%   SOC(1) = SOC_init + eta_ch * P_ch(1) * dt - P_dis(1)/eta_dis * dt
%   SOC(t) - SOC(t-1) = eta_ch * P_ch(t) * dt - P_dis(t)/eta_dis * dt
```

### 5.3 Piecewise Linear Cost via Epigraph

For convex piecewise linear functions (e.g., carbon trading penalty), use epigraph formulation:
```matlab
% z >= s_k * x + b_k   for each segment k
% Minimize sum(z)
%
% add_var('z', N, ...)
% add_lin_constraint: A * [x_vars; z] >= [b_1; b_2; ...; b_K] (per unit, per segment)
% add_quad_cost: C = [0; ...; 0; 1; ...; 1]  (only on z)
```

No integer variables needed for convex PWL.

### 5.4 Curtailment Penalty for Renewables

For renewables with zero gencost, add a negative linear cost on Pg:
```matlab
% Penalize curtailment by rewarding dispatch
C_wind = -penalty * dt / baseMVA;   % note: Pg is in p.u. inside formulation
om.add_quad_cost('wind_value', sparse(n,n), C_wind * ones(n,1), 0, {'Pg'});
```

### 5.5 Cost Term Design

- `add_quad_cost` with Q=0 → pure linear cost: C'*x
- Negative C → pushes variable upward (revenue term)
- gencost scaling by `modcost(gc, time)` converts $/h to $/period
- Custom costs via `add_quad_cost` are NOT auto-scaled by period, must include Δt manually
- For custom variables in physical units (not p.u.), no baseMVA conversion needed in cost

---

## 6. Common Pitfalls & Debugging

### 6.1 args Are Frozen Snapshots

The `args` passed at `add_userfcn(..., @callback, args)` registration time are **frozen**.
You cannot modify args in the formulation callback and have int2ext see the changes.
Use the `results` struct (or `mpc.pofc`-style fields) for cross-callback data sharing.

### 6.2 define_constants vs Nested Functions

MATLAB's `define_constants` uses `assignin('caller', ...)` which **fails** in functions
containing nested functions (static workspace). Solution: use **local functions** (defined
at end of file with `function ... end`) instead of nested functions in any file that calls
`define_constants`.

### 6.3 Internal vs External Generator Count

After `ext2int()`, out-of-service generators may be removed. Always use:
```matlab
% In formulation callback:
ng_internal = om.getN('var', 'Pg') / nt;   % NOT args.ng_ori

% In int2ext callback:
ng_pp = size(results.gen, 1) / nt;          % NOT results.pofc.ng_ori
```

### 6.4 case9 gen Matrix Column Count

Standard case9 has only 21 columns in gen matrix. `RAMP_30` is column 23 (from `idx_gen`),
which **doesn't exist**. Solution: pass ramp limits via your own struct field
(e.g., `mpc.pofc.ramp_limit`) instead of using `gen(:, RAMP_30)`.

### 6.5 IPOPT "Local Infeasibility"

Common causes:
- Ramp limits too tight: total ramp capacity < inter-period load change
- Variable bounds contradicting equality constraints
- Missing baseMVA conversion making constraints inconsistent

**Debug approach:**
1. Relax bounds/constraints one group at a time
2. Check constraint feasibility with a known feasible point
3. Set `mpopt = mpoption('verbose', 2)` for detailed IPOPT output

### 6.6 Backward Compatibility Pattern

Gate new features with `isfield()` checks:
```matlab
if isfield(args, 'my_new_feature') && ~isempty(args.my_new_feature)
    % new feature path
else
    % legacy path (must still work!)
end
```

---

## 7. Recommended MATPOWER Source Files to Read

| File | Why |
|------|-----|
| `opf_setup.m` | How om is built, standard variables/constraints added |
| `opf_execute.m` | Post-processing, how results.var.val.* is populated |
| `nlpopf_solver.m` | How om.solve() calls IPOPT, p.u. ↔ MW conversion |
| `opf_power_balance_fcn.m` | Reference for NLN constraint callback pattern (Jacobian format) |
| `run_userfcn.m` | How callbacks are dispatched |
| `ext2int.m` / `int2ext.m` | External ↔ internal ordering conversion |
| `modcost.m` | How gencost scaling works |
| `makeAvl.m` | Dispatchable load handling (Q limit requirements) |
| `define_constants.m` | Where PG, PD, GEN_BUS, etc. come from |

---

*This guide is based on MATPOWER 7.1 with IPOPT solver. Some details may differ in other versions.*
