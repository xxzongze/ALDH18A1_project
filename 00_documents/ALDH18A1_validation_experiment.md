# ALDH18A1/Glu 代谢分叉验证——详细操作手册

---

## Phase 1：现有样本补测 + 蛋白验证

### 已有数据来源
- 肝转录组 (RNA-seq, DLY/TFB, 4 阶段) 和血清 AA 谱（已有 17 种 AA，Pro 未测定）来自前期实验
- 血清 AA 定量采用 UPLC-MS/MS 靶向代谢组学方法 [1,2]

---

### Exp 1.1：血清 Pro 及其他 AA 补测

**样本清单**：

| 品种 | 阶段 | n | 样本量 |
|------|------|---|--------|
| DLY | 15 kg | 8 | 8 管 |
| DLY | 45 kg | 8 | 8 管 |
| DLY | 75 kg | 8 | 8 管 |
| DLY | 105 kg | 8 | 8 管 |
| TFB | 15 kg | 8 | 8 管 |
| TFB | 45 kg | 8 | 8 管 |
| TFB | 75 kg | 8 | 8 管 |
| TFB | 105 kg | 8 | 8 管 |
| **合计** | | **64 管** | |

每管需血清 ≥50 μL。从 -80℃ 取出后冰上解冻，分装 50 μL 至 EP 管，避免反复冻融。

**送测 AA 清单**（至少包含以下 18 种）：

```
必需类: Leu, Ile, Val, Lys, Met, Thr, Phe, Trp, His, Arg
非必需类: Pro, Hyp, Glu, Gln, Orn, Cit, Gly, Ser
```

选测：Cys, Tau, Tyr, Asp, Asn, Ala, P5C（如果能找到标准品）

**送测公司参考**：华大基因靶向代谢组、诺禾致源、百趣生物、中科新生命——任选一家，报"18 种血清游离氨基酸 UPLC-MS/MS 靶向定量"，约 2000-4000 元。

**送测前准备**：
1. 样本列表（编号+品种+阶段）发给公司
2. 确认公司能做 Pro 和 Hyp（有些标准 AA panel 不含 Pro，需要单独加）
3. 确认 Orn 和 Cit 在 panel 中（有些公司把这两种归为"尿素循环代谢物"需要单独加）

**数据分析**：

预期结果表头：
| AA | 15kg_DLY | 15kg_TFB | p_15 | 45kg_DLY | 45kg_TFB | p_45 | 75kg_DLY | 75kg_TFB | p_75 | 105kg_DLY | 105kg_TFB | p_105 |
|-----|----------|----------|------|----------|----------|------|----------|----------|------|-----------|-----------|-------|

重点分析指标：
- **Pro 绝对浓度**：DLY vs TFB 各阶段 Welch's t-test
- **Pro/Glu 比值**：反映 P5CS 将 Glu 转化为 Pro 的效率
- **Glu/(Orn+Cit) 比值**：反映 Glu 走向合成代谢 vs 尿素循环的分叉比例——这是最核心的分叉指标
- **Urea**：作为氮流失的终产物对照

作图：
- 折线图：Pro, Glu, Orn, Cit, Urea 在 4 阶段的 DLY vs TFB 比较
- 散点图：Glu/(Orn+Cit) 比值 vs 蛋白沉积表型（用你现有的蛋白沉积数据）

---

### Exp 1.2：肝组织蛋白验证

**样本选择**：优先用 75 kg 阶段的肝组织（蛋白沉积差异最大的阶段）。如果样本量够，加做 105 kg。

DLY n=6, TFB n=6, 共 12 个肝组织样本。

**组织裂解**（参考标准 WB 方法 [3,4]）：
- 取 ~50 mg 肝组织，加入 500 μL RIPA 裂解液（含蛋白酶抑制剂 cocktail + 磷酸酶抑制剂）
- 使用组织研磨仪：60 Hz, 60 s × 2 次（冰上操作）
- 冰上裂解 30 min，每 10 min 涡旋一次
- 4°C, 12000g 离心 15 min，取上清
- BCA 法定量蛋白浓度 [5]

**WB 条件**：

| 靶蛋白 | 货号 | 分子量 | 一抗稀释 | 二抗 | 上样量 | 转膜条件 |
|--------|------|--------|---------|------|--------|---------|
| ALDH18A1 | Proteintech 17719-1-AP | ~87 kDa | 1:1000 | 兔二抗 1:5000 | 30 μg | 300 mA, 90 min (0.45 μm PVDF) |
| PYCR1 | Proteintech 13108-1-AP | ~34 kDa | 1:1000 | 兔二抗 1:5000 | 20 μg | 300 mA, 45 min |
| OAT | Santa Cruz sc-398700 | ~49 kDa | 1:500 | 鼠二抗 1:5000 | 25 μg | 300 mA, 60 min |
| CPS1 | Abcam ab45956 | ~165 kDa | 1:2000 | 兔二抗 1:5000 | 30 μg | 300 mA, 120 min (0.45 μm PVDF) |
| β-actin | CST #4970 | ~45 kDa | 1:5000 | 兔二抗 1:5000 | — | — |
| GAPDH | CST #5174 | ~37 kDa | 1:5000 | 兔二抗 1:5000 | — | — |

**备选抗体**：若 Proteintech 17719-1-AP 无货，可用 Abcam ab177793（ALDH18A1, ~87 kDa, 1:1000）。

**ALDH18A1 注意事项**：
- 该蛋白在肝脏中高表达，信号应较强
- 若出现多条带，主带应在 ~87 kDa
- 可用人 HepG2 裂解液做阳性对照

**内参选择**：
- ALDH18A1 和 CPS1：用 β-actin（两者分子量与 GAPDH 距离较远）
- PYCR1 和 OAT：用 GAPDH（分子量接近）

**WB 操作**：
1. SDS-PAGE：8% 分离胶（CPS1 和大分子）/ 10% 分离胶（其余）
2. 转膜后用丽春红染色确认转膜效率
3. 5% BSA/TBST 封闭 1h，室温
4. 一抗 4°C 孵育过夜
5. TBST 洗 3×10 min
6. 二抗室温 1h
7. TBST 洗 3×10 min
8. ECL 显影（Bio-Rad ChemiDoc 或压片）

**预期结果表**：

| 基因 | DLY75 (n=6) | TFB75 (n=6) | DLY/TFB Ratio | p value |
|------|-------------|-------------|---------------|---------|
| ALDH18A1/β-actin | Mean±SEM | Mean±SEM | — | — |
| PYCR1/β-actin | — | — | — | — |
| OAT/GAPDH | — | — | — | — |
| CPS1/β-actin | — | — | — | — |

**预期**：DLY 组 ALDH18A1↑, PYCR1↑, CPS1↓；TFB 组 OAT↑, CPS1↑。

**如果 ALDH18A1 蛋白差异不明显**：
- 先确认 qPCR 验证 mRNA 是否有差异（用现有 RNA-seq 数据的 log2FC 做参考）
- 若 mRNA 有差异而蛋白无差异 → 可能存在翻译后调控，补做泛素化检测或考虑 mRNA 稳定性差异
- 若 mRNA 也无差异 → 假说需要修正，考虑 ALDH18A1 的酶活性差异（需要测酶活，比较复杂）

---

### Exp 1.3：肌肉 mTORC1 通路

**样本**：DLY/TFB 各 n=6，75 kg 肌肉（与肝组织同一批动物）。

**WB 条件**：

| 靶蛋白 | 货号 | 分子量 | 一抗稀释 | 上样量 |
|--------|------|--------|---------|--------|
| p-S6K1 (T389) | CST #9205 | ~70 kDa | 1:1000 | 25 μg |
| S6K1 | CST #2708 | ~70 kDa | 1:1000 | 25 μg |
| p-4E-BP1 (T37/46) | CST #2855 | ~15-20 kDa | 1:1000 | 25 μg |
| 4E-BP1 | CST #9644 | ~15-20 kDa | 1:1000 | 25 μg |
| p-AKT (S473) | CST #4060 | ~60 kDa | 1:2000 | 20 μg |
| AKT | CST #4691 | ~60 kDa | 1:2000 | 20 μg |
| β-actin | CST #4970 | ~45 kDa | 1:5000 | — |

**4E-BP1 特殊注意事项**：
- 4E-BP1 分子量小（~15-20 kDa），需用 15% 分离胶或梯度胶（4-20%）
- 转膜用 0.2 μm PVDF 膜（不能用 0.45 μm，小蛋白会穿透）
- 转膜条件：100V, 60 min 或 300 mA, 45 min
- p-4E-BP1 通常显示为 3 条带（α/β/γ 磷酸化形式），定量时取所有条带的总和

**磷酸化蛋白操作要点**：
1. 裂解液中必须含磷酸酶抑制剂（NaF, Na3VO4, β-glycerophosphate）
2. 所有操作在冰上进行
3. 5% BSA/TBST 封闭（不能用脱脂奶粉，奶粉含磷酸酶）
4. 一抗用 5% BSA/TBST 稀释
5. 先跑磷酸化抗体 → strip → 再跑 total 抗体（确保配对信号来自同一条膜）

**定量方法**：
- 目的条带灰度值 / 内参灰度值 → 相对表达量
- p-S6K1 / total S6K1 → S6K1 磷酸化率
- p-4E-BP1 / total 4E-BP1 → 4E-BP1 磷酸化率

**预期结果表**：

| 指标 | DLY75 (n=6) | TFB75 (n=6) | p value |
|------|-------------|-------------|---------|
| p-S6K1/S6K1 | — | — | — |
| p-4E-BP1/4E-BP1 | — | — | — |
| p-AKT/AKT | — | — | — |

**预期**：p-S6K1/S6K1 和 p-4E-BP1/4E-BP1 在 DLY 组显著高于 TFB。p-AKT/AKT 可能无显著差异——这正是"非经典 mTORC1 激活"的支持证据。

**相关性分析**（待 Exp 1.1 Pro 数据出来后做）：
- 血清 Pro vs 肌肉 p-S6K1/S6K1（Pearson/Spearman）
- 血清 Glu/(Orn+Cit) vs 肌肉 p-S6K1/S6K1
- 预期：血清 Pro 与 p-S6K1 正相关，血清 Orn+Cit 与 p-S6K1 负相关

---

## Phase 2：肝细胞机制 + 肌管串扰

---

### Exp 2.1：肝细胞 ALDH18A1 敲低

#### 2.1.1 猪原代肝细胞分离与培养

**动物**：DLY 猪，~15-20 kg，n=3（3 次独立分离 = 3 生物学重复）。

**分离方法（两步胶原酶灌注法 [8,9]）**：

1. 麻醉后开腹，暴露门静脉
2. 插管（18G 留置针），37°C 预热的灌注液 I（HBSS + 0.5 mM EGTA, pH 7.4）按 30 mL/min 灌注 10-15 min
3. 切换灌注液 II（HBSS + 5 mM CaCl2 + 0.05% collagenase IV），继续灌注 10-15 min 至肝包膜下出现裂隙
4. 取出肝脏到无菌平皿，撕开包膜，用预冷 William's E 培养基冲洗出肝细胞
5. 100 μm 细胞筛过滤
6. 50g 离心 3 min × 3 次（低速去除非实质细胞）
7. 台盼蓝计数：活率 ≥85% 可用
8. 接种至胶原包被的 6 孔板（1×10^6 cells/well）
9. 贴壁培养基：William's E + 10% FBS + 1% PS + 10 nM insulin + 100 nM dexamethasone
10. 37°C, 5% CO2，贴壁 4-6h 后换维持培养基

**维持培养基**：William's E + 1% FBS + 1% PS + 10 nM insulin + 0.5 mM Glu + 2 mM Gln

**肝细胞工作量估算**：
- 6 孔板，每孔 1×10^6 细胞
- siRNA-NC：3 孔 × 3 次独立分离 = 9 孔
- siRNA-A1（序列 1）：3 孔 × 3 = 9 孔
- siRNA-A2（序列 2）：3 孔 × 3 = 9 孔
- 每孔收细胞裂解液（WB）+ 培养基（代谢物）
- 共需 27 孔 = 5 块 6 孔板
- 每次分离至少需要 1×10^6×9 = 9×10^6 活细胞

**备选：若原代肝细胞不可行**

使用 **HepG2 细胞系**：
- DMEM + 10% FBS + 1% PS
- 先检测 HepG2 内源性 ALDH18A1 表达（qPCR+WB）确认表达量足够
- 若表达量低，用 pcDNA3.1-猪ALDH18A1 过表达质粒转染 → 建立稳定表达株
- 然后用 siRNA 敲低（敲低人 ALDH18A1，但可能和猪序列不完全匹配）
- 或直接用猪 ALDH18A1 过表达 HepG2 做"OE vs EV"的对比（绕过 siRNA）

#### 2.1.2 siRNA 设计与转染

**猪 ALDH18A1 序列获取**：
1. 从 Ensembl (Sus scrofa 11.1) 获取 ALDH18A1 CDS
2. 设计 3 条 siRNA（推荐各合成 5 nmol）：
   - 序列 1：靶向 CDS 区前 1/3
   - 序列 2：靶向 CDS 区后 1/3
   - 序列 3（备用）：靶向 3'UTR
3. 阴性对照：通用 scrambled siRNA（公司提供）
4. 合成公司：苏州吉玛/广州锐博/上海吉凯——报"猪 ALDH18A1 siRNA 设计+合成"

**转染方案（以 Lipofectamine RNAiMAX 为例）**：

贴壁 24h 后（细胞融合度 ~60-70%）：
1. 每孔：siRNA 终浓度 50 nM
2. 管 A：3 μL RNAiMAX + 150 μL Opti-MEM，室温 5 min
3. 管 B：siRNA（50 pmol）+ 150 μL Opti-MEM
4. 管 A + 管 B 混合，室温孵育 20 min
5. 滴加至每孔（补 700 μL 维持培养基至总量 1 mL）
6. 转染 6h 后换新鲜维持培养基
7. 48h 后收样

**KD 效率检测**：
- mRNA：转染后 24h 收 RNA（TRIzol）→ qPCR（引物设计见下方）
- 蛋白：转染后 48h 收蛋白裂解液 → WB

**ALDH18A1 qPCR 引物**（需自行验证熔解曲线）：

| 引物 | 序列 (5'→3') | 产物长度 |
|------|-------------|---------|
| ALDH18A1-F | GCTGAAGGTGATTGGCACTG | ~120 bp |
| ALDH18A1-R | TCCACATAGTCGATGCCGTC | |
| β-actin-F | CTCTTCCAGCCTTCCTTCCT | ~180 bp |
| β-actin-R | AGCACTGTGTTGGCGTACAG | |

**筛选最佳 siRNA**：
- 选 mRNA KD 效率 > 70% 的序列
- 若有 2 条都满足，选蛋白 KD 效率更高的
- 后续实验只用最佳序列

#### 2.1.3 实验分组与处理

| 组别 | 处理 | 目的 | 重复数 |
|------|------|------|--------|
| NC | siRNA 阴性对照 | 基线 | n=3 |
| KD | siRNA-ALDH18A1 (最佳序列) | 验证分叉 | n=3 |
| KD+Pro | siRNA + L-Proline 0.2 mM (转染后 24h 添加) | 产物回补 | n=3 |
| KD+B6 | siRNA + Pyridoxal 5'-phosphate 10 μM | 辅因子代偿 | n=3 |

- 重复数：n=3 = 3 次独立肝细胞分离（生物学重复）
- 每次分离：4 组 × 1 孔 = 4 孔/板

**Pro 和 B6 处理细节**：
- L-Proline (Sigma P0380)：溶于 PBS，母液 100 mM，工作浓度 0.2 mM
  - 0.2 mM Pro = 接近猪血清生理浓度（文献值 ~0.15-0.25 mM）
- Pyridoxal 5'-phosphate (Sigma P9255)：即维生素 B6 的活性形式，溶于培养基
  - 母液 10 mM，工作浓度 10 μM
- Pro 和 B6 在转染后 **24h 开始添加**，维持 24h 后收样

#### 2.1.4 收样与检测

**收样时间点**：

| 时间 | 操作 |
|------|------|
| 转染后 0h | 换维持培养基 |
| 转染后 6h | 换新鲜维持培养基 |
| 转染后 24h | 收 RNA（TRIzol 裂解，用于 KD 效率验证）|
| 转染后 24h | KD+Pro 和 KD+B6 组开始加药 |
| 转染后 48h | 收培养基（600 μL，-80°C 冻存待测） |
| 转染后 48h | 收细胞蛋白（RIPA 裂解，-80°C 冻存） |

**培养基检测指标**（送 UPLC-MS/MS 或生化试剂盒）：

| 指标 | 方法 | 最小送样量 |
|------|------|----------|
| Pro | UPLC-MS/MS | 100 μL |
| Glu | UPLC-MS/MS 或谷氨酸检测试剂盒 | 100 μL |
| Orn | UPLC-MS/MS | 100 μL |
| Urea | 脲酶法试剂盒 (Sigma MAK006) | 25 μL |

**试剂盒替代方案**：
- 谷氨酸：Sigma MAK004 (Glutamate Assay Kit)，荧光法，96 孔板
- 尿素：Sigma MAK006 (Urea Assay Kit)，比色法
- 这 2 个试剂盒可以自己在实验室测，不需要送公司，节省费用和时间

**培养基收样要点**：
1. 收 600 μL 培养基至 EP 管
2. 4°C, 3000g 离心 5 min 去除细胞碎片
3. 转移上清至新 EP 管
4. -80°C 冻存
5. 集中一批送测（不要分批）

**细胞裂解液检测**（WB）：

| 靶蛋白 | 目的 |
|--------|------|
| ALDH18A1 | 确认 KD 效率（蛋白水平） |
| PYCR1 | 检测代谢代偿（Pro 合成旁路是否上调） |
| OAT | 检测代谢代偿（Orn 合成是否上调） |
| CPS1 | 检测尿素循环是否被激活 |
| β-actin | 内参 |

#### 2.1.5 预期数据表

**培养基代谢物组（正常化至细胞蛋白量）**：

| 指标 | NC | KD | KD+Pro | KD+B6 | NC vs KD p |
|------|-----|-----|--------|-------|-----------|
| Pro (μM/mg protein) | — | ↓ | — | — | — |
| Glu (μM/mg protein) | — | ↑ | — | — | — |
| Orn (μM/mg protein) | — | ↑ or — | — | — | — |
| Urea (μM/mg protein) | — | ↑ | — | — | — |
| Pro/Glu Ratio | — | ↓ | — | — | — |
| Glu/(Orn+Cit) Ratio | — | ↓ | — | — | — |

**核心验证点**：
- KD → Pro↓ 且 Urea↑ = 分叉方向改变（关键）
- KD+Pro 组：Pro 恢复至接近 NC 水平（证明 Pro 是直接产物）
- KD+B6 组：若部分恢复 → B6（PLP）是 ALDH18A1 的辅因子，可做营养干预靶点
- KD → Pro/Glu 比值↓（Glu 向 Pro 转化效率下降）
- KD → Glu/(Orn+Cit) 比值↓（Glu 向尿素循环方向偏移）

---

### Exp 2.2：条件培养基-肌管串扰实验

#### 2.2.1 实验逻辑

```
肝细胞 NC siRNA 48h → 收集 CM-NC
肝细胞 KD siRNA 48h → 收集 CM-KD
                            ↓
             加到 C2C12 分化肌管上 24h
                            ↓
              检测蛋白合成 + 降解 + 信号通路
```

#### 2.2.2 条件培养基制备

**与 Exp 2.1 同步进行**：

1. 肝细胞转染 NC/KD siRNA（各 3 孔生物学重复）
2. 转染后 48h，收集培养基（每孔 ~1 mL）
3. 4°C, 3000g 离心 5 min，去除细胞碎片
4. 合并同组 3 孔的培养基（得 ~3 mL CM-NC, ~3 mL CM-KD）
5. 取 200 μL 送 AA 谱测定（确认 Pro 差异）
6. 剩余培养基分装 500 μL/管，-80°C 冻存
7. 使用时冰上解冻，与新鲜肌管分化培养基 1:1 混合

**不要反复冻融 CM**。

#### 2.2.3 C2C12 肌管分化 [10,11]

**细胞**：C2C12 小鼠成肌细胞（ATCC CRL-1772）

**分化方案**：
1. 增殖培养基：DMEM + 10% FBS + 1% PS
2. 接种至 12 孔板（预铺 0.1% gelatin）：5×10^4 cells/well
3. 达到 90-100% 融合时（约 48h），换分化培养基
4. 分化培养基：DMEM + 2% horse serum + 1% PS
5. 每 24h 换分化培养基
6. 第 4-5 天可见明显肌管（多核，>5 核/管）

**肌管质量判断**：
- 镜下：长条形、多核（DAPI 染色可见 ≥3 核）、有明显肌节
- 融合指数（MyHC+ 细胞中 ≥2 核的比例）> 70%
- 若分化率 <50%，重新复苏低代次 C2C12（<P8）

#### 2.2.4 实验分组

| 组别 | 处理 | 目的 |
|------|------|------|
| ① NC-CM | 1:1 CM-NC + 分化培养基 | 基线 |
| ② KD-CM | 1:1 CM-KD + 分化培养基 | 验证 CM 效应 |
| ③ KD-CM+Pro | 1:1 CM-KD + 分化培养基 + 0.2 mM Pro | Pro 救援（关键！） |
| ④ NC-CM+Urea | 1:1 CM-NC + 分化培养基 + 2 mM Urea | 排除尿素毒性 |

Pro 浓度 0.2 mM = 接近猪血清生理水平。
Urea 浓度 2 mM = 参考猪血清尿素范围（2-5 mM），取较低值模拟 KD 后尿素可能的升高幅度。

**技术重复**：每组 3 孔（独立培养孔），生物学重复来自 3 次独立 CM 制备。

#### 2.2.5 检测指标与方法

**A. 蛋白合成——SUnSET 法** [6,7]

原理：puromycin 是氨基酰-tRNA 类似物，被掺入新生肽链。抗 puromycin 抗体检测 puromycin 掺入量 = 蛋白合成速率 [6]。

操作：
1. 肌管 + CM 处理 24h
2. 收样前 30 min，加入 puromycin（终浓度 1 μM）到培养液中
3. 30 min 后，PBS 洗 2 次
4. RIPA 裂解（含磷酸酶抑制剂）
5. BCA 定量
6. WB：上样 20 μg，anti-puromycin (Kerafast EQ0001, 1:5000)
7. puromycin 掺入信号覆盖全泳道（所有正在合成的新生肽链）
8. 全泳道灰度值 / 内参 = 相对蛋白合成速率

**B. mTORC1 通路**

同 Exp 1.3：p-S6K1/S6K1, p-4E-BP1/4E-BP1

**C. 肌管直径**

MyHC 免疫荧光染色：
1. 4% PFA 固定 15 min
2. 0.1% Triton X-100 透化 10 min
3. 5% BSA 封闭 30 min
4. anti-MyHC (MF20, DSHB, 1:100) 4°C 过夜
5. 二抗：anti-mouse IgG-AF488 (1:500) 室温 1h
6. DAPI 封片
7. 荧光显微镜拍照（每孔随机 5 个视野，20×）
8. ImageJ 测量肌管直径：每视野至少测 20 根肌管的最短横径

**肌管直径测量要点**：
- ImageJ → Straight line tool → 垂直于肌管长轴画线 → Measure
- 每孔至少测 100 根肌管
- 只测 MyHC+ 且 ≥3 个核的肌管
- 数据格式：肌管编号、所属组别、直径(μm)、视野编号

#### 2.2.6 预期结果

| 指标 | ①NC-CM | ②KD-CM | ③KD-CM+Pro | ④NC-CM+Urea |
|------|--------|--------|------------|-------------|
| SUnSET (蛋白合成) | = 1.0 | ↓ ~30% | 恢复至 ~85% | 无变化或轻微↓ |
| p-S6K1/S6K1 | — | ↓ | 恢复 | — |
| 肌管直径 (μm) | — | ↓ ~15% | 恢复 | 无变化 |
| CM 中 Pro 浓度 | 正常 | ↓ | — | 正常 |

**关键解释**：
- ② vs ①：KD-CM 降低蛋白合成 → 肝 ALDH18A1 活性通过分泌组影响肌肉
- ③ vs ②：Pro 回补恢复合成 → **Pro 是效应分子**（最关键的结论）
- ④ vs ①：Urea 不影响合成 → KD-CM 的效应不是尿素引起的毒性
- 如果 ④ 也导致合成↓ → 需要考虑尿素有独立效应（这是额外的发现）

---

## Phase 3：饲粮 Pro+Glu 补充验证

---

### Exp 3.1：2×2 因子饲养试验

#### 3.1.1 试验设计

| 因子 | 水平 1 | 水平 2 |
|------|--------|--------|
| 品种 (Breed) | DLY | TFB |
| 饲粮 (Diet) | 基础日粮 (CON) | 基础+0.5%Pro+1%Glu (SUPP) |

4 个处理组：DLY-CON, DLY-SUPP, TFB-CON, TFB-SUPP

| 项目 | 内容 |
|------|------|
| 动物数 | 每组 n=8，共 32 头（公母各半或全公） |
| 起始体重 | 25 ± 2 kg |
| 结束体重 | 70-75 kg |
| 周期 | ~50-60 天 |
| 饲养方式 | 单栏饲养，自由采食+饮水 |

**统计功效估算**：
- 假设 SUPP 效应量 = 0.8 SD（中等偏大）
- α=0.05, power=0.80, 2×2 factorial
- 每组需 n≈8（考虑到可能的脱落，多备 1-2 头）

#### 3.1.2 饲粮配方要点

**基础日粮**：玉米-豆粕型，粗蛋白 ~16%（生长阶段）/ ~14%（肥育阶段）

**Pro 和 Glu 添加**：
- L-Proline (饲料级, ≥98%): 0.5% = 5 g/kg 饲粮
  - 参考：文献中 Pro 补充范围 0.3-1.0%
  - 选择 0.5% 作为中间剂量（安全且能产生效应）
- L-Glutamate (饲料级, ≥98%): 1.0% = 10 g/kg 饲粮
  - Glu 在常规饲料原料中含量丰富，额外补充 1% 约提升总 Glu 摄入 10-15%
- 替代等量玉米（调整配方保持能量和蛋白平衡）
- 或不替代：直接添加在预混料中（添加量小，对整体配方影响不大）

**与饲料厂沟通要点**：
1. Pro 和 Glu 需要混匀至整个批次
2. 确认 Pro 和 Glu 在制粒温度下稳定（两者的分解温度均 >200°C，制粒安全）
3. 建议分 3-4 批做料（避免一批用完）

#### 3.1.3 检测计划

| 时间点 | 检测项目 | 方法 |
|--------|---------|------|
| 试验开始 | 初始体重 | 称重 |
| 每 14 天 | 体重、采食量 | 称重+记录 |
| 第 30 天 | 氮平衡（全收粪尿法） | 代谢笼 5 天 |
| 第 55-60 天 | 氮平衡（全收粪尿法） | 代谢笼 5 天 |
| 每 14 天 | 血清采集（前腔静脉） | 5 mL 促凝管 |
| 屠宰日 | 胴体重、眼肌面积、背膘厚 | 常规屠宰测定 |
| 屠宰日 | 肝组织（~5g, 液氮速冻） | 术中采样 |
| 屠宰日 | 肌肉组织（~5g, 液氮速冻） | 术中采样 |

**血清采集时间点**（最少 4 个）：
- 试验第 0 天（基线）
- 第 28 天（中期）
- 第 42 天
- 第 56 天（结束前）

每头每次采 5 mL 全血 → 促凝管 → 室温静置 30 min → 3000g 离心 15 min → 分装血清 200 μL×3 管 → -80°C

**血清检测**：
- AA 谱（同 Phase 1，含 Pro, Glu, Orn, Cit, BCAA）
- Urea, ALT, AST（自动生化仪）
- 可选：IGF1 ELISA（猪 IGF1 Quantikine ELISA, R&D Systems）

**氮平衡全收粪尿法**（你们组的成熟方法）：
- 代谢笼饲养 5 天（3 天适应 + 2 天收集，但你们如果已有标准方案就用你们的）
- 每日记录：采食量、排粪量、排尿量
- 取样：粪（10% 每日取样混合）、尿（10% 每日取样混合，加 H2SO4 防氨挥发）
- 凯氏定氮法测粪 N、尿 N、饲料 N

**组织采样**：
- 肝：右叶边缘 ~5g，分 2 份（WB 用 + 备存），液氮速冻
- 肌肉（背最长肌）：第 10-12 肋间 ~5g，分 2 份，液氮速冻

**屠宰后检测**：

肝 WB（同 Exp 1.2）：ALDH18A1, PYCR1, OAT, CPS1
肌肉 WB（同 Exp 1.3）：p-S6K1/S6K1, p-4E-BP1/4E-BP1

#### 3.1.4 统计分析

**主要分析**（SAS PROC MIXED 或 R lme4）：

```r
# R 示例代码框架
model <- lmer(ADG ~ Breed * Diet + (1|Block), data = data)
anova(model)
```

固定效应：Breed, Diet, Breed×Diet
随机效应：Block（如按体重区组）

**关键假设检验**：
1. **交互效应 Breed×Diet**：如果 p<0.05，说明 TFB 和 DLY 对补充的响应不同
2. **简单效应**：在 TFB 内部检验 CON vs SUPP 的差异（t-test 或 contrast）
3. **预期结果**：
   - 交互显著：TFB-SUPP > TFB-CON 的增幅 > DLY-SUPP > DLY-CON 的增幅
   - 若 DLY-SUPP ≈ DLY-CON：DLY 内源性通路已饱和，补充无效
   - 若 TFB-SUPP ≈ DLY-CON：补充使 TFB 达到 DLY 水平（**最理想的结果**）

**主要表格**：

| 指标 | DLY-CON | DLY-SUPP | TFB-CON | TFB-SUPP | SEM | p_Breed | p_Diet | p_Inter |
|------|---------|----------|---------|----------|-----|---------|--------|---------|
| ADG (g/d) | | | | | | | | |
| ADFI (g/d) | | | | | | | | |
| F/G | | | | | | | | |
| N 沉积 (g/d) | | | | | | | | |
| NPU (%) | | | | | | | | |
| 血清 Pro (μM) | | | | | | | | |
| 血清 Urea (mM) | | | | | | | | |
| 胴体重 (kg) | | | | | | | | |

---

## 通用实验条件

### WB 通用 Buffers

**RIPA 裂解液**（100 mL）：
- 50 mM Tris-HCl (pH 7.4): 0.6057 g Tris
- 150 mM NaCl: 0.8766 g
- 1% NP-40: 1 mL
- 0.5% sodium deoxycholate: 0.5 g
- 0.1% SDS: 0.1 g
- 补水至 100 mL
- 用前加：蛋白酶抑制剂 cocktail (Roche, 1 tablet/50 mL) + 1 mM PMSF + 磷酸酶抑制剂 (1 mM NaF, 1 mM Na3VO4, 10 mM β-glycerophosphate)

**10× Running Buffer**（1 L）：
- Tris: 30.3 g, Glycine: 144 g, SDS: 10 g, 补水至 1 L

**10× Transfer Buffer**（1 L）：
- Tris: 30.3 g, Glycine: 144 g, 补水至 1 L
- 1× 工作液：100 mL 10× + 200 mL 甲醇 + 700 mL 水

**TBST**（1 L）：
- 20 mM Tris (pH 7.6): 2.42 g
- 137 mM NaCl: 8 g
- 0.1% Tween-20: 1 mL

### 试剂采购清单（汇总）

| 试剂/耗材 | 货号 | 用途 | 预估价格 |
|-----------|------|------|---------|
| Proteintech ALDH18A1 Ab | 17719-1-AP | 肝/细胞 WB | ~2500 |
| Proteintech PYCR1 Ab | 13108-1-AP | 肝 WB | ~2500 |
| Santa Cruz OAT Ab | sc-398700 | 肝 WB | ~3000 |
| Abcam CPS1 Ab | ab45956 | 肝 WB | ~3500 |
| CST p-S6K1 Ab | #9205 | 肌肉/肌管 WB | ~3000 |
| CST S6K1 Ab | #2708 | 肌肉/肌管 WB | ~3000 |
| CST p-4E-BP1 Ab | #2855 | 肌肉/肌管 WB | ~3000 |
| CST 4E-BP1 Ab | #9644 | 肌肉/肌管 WB | ~3000 |
| CST p-AKT Ab | #4060 | 肌肉 WB | ~3000 |
| CST AKT Ab | #4691 | 肌肉 WB | ~3000 |
| CST β-actin Ab | #4970 | WB 内参 | ~2500 |
| Kerafast anti-puromycin | EQ0001 | SUnSET | ~2500 |
| MF20 MyHC Ab | DSHB | IF | ~2000 |
| 猪 ALDH18A1 siRNA ×3 | 定制 | KD | ~3000 |
| Lipofectamine RNAiMAX | 13778-150 | 转染 | ~3000 |
| Puromycin | P8833 | SUnSET | ~500 |
| L-Proline | P0380 | 救援实验 | ~300 |
| PLP (B6) | P9255 | 辅因子实验 | ~500 |
| Urea Assay Kit | MAK006 | 尿素检测 | ~2500 |
| Glutamate Assay Kit | MAK004 | Glu 检测 | ~2500 |
| Collagenase IV | C5138 | 肝细胞分离 | ~2000 |

---

## 论文 Figure 规划

| Figure | 内容 | 来源 |
|--------|------|------|
| Fig 1 | 肝 ALDH18A1/PYCR1/OAT/CPS1 蛋白 + 血清 Glu/Pro/Orn/Cit/Urea 4 阶段折线图 | Phase 1 |
| Fig 2 | 肌肉 p-S6K1/p-4E-BP1 + 血清 Pro vs p-S6K1 相关性散点图 | Phase 1 |
| Fig 3 | 肝细胞 KD 后培养基 Pro↓/Urea↑ + Pro 救援 + 代偿变化 | Phase 2.1 |
| Fig 4 | CM 串扰：SUnSET + 肌管直径 + p-S6K1 + Pro 救援 | Phase 2.2 |
| Fig 5 | 饲粮补充试验：ADG/N沉积/血清Pro/Urea 的 2×2 交互图 | Phase 3 |
| Fig 6 | 工作模型：Glu 代谢分叉决定氮分配 | 综合 |

---

## 方法参考文献

1. Dettmer K, Aronov PA, Hammock BD. Mass spectrometry-based metabolomics. *Mass Spectrometry Reviews*, 2007, 26(1): 51-78.
2. Wishart DS, Jewison T, Guo AC, et al. HMDB 3.0 — The Human Metabolome Database in 2013. *Nucleic Acids Research*, 2013, 41(D1): D801-D807.
3. Towbin H, Staehelin T, Gordon J. Electrophoretic transfer of proteins from polyacrylamide gels to nitrocellulose sheets: procedure and some applications. *Proceedings of the National Academy of Sciences*, 1979, 76(9): 4350-4354.
4. Burnette WN. "Western blotting": electrophoretic transfer of proteins from sodium dodecyl sulfate-polyacrylamide gels to unmodified nitrocellulose and radiographic detection with antibody and radioiodinated protein A. *Analytical Biochemistry*, 1981, 112(2): 195-203.
5. Smith PK, Krohn RI, Hermanson GT, et al. Measurement of protein using bicinchoninic acid. *Analytical Biochemistry*, 1985, 150(1): 76-85.
6. Schmidt EK, Clavarino G, Ceppi M, Pierre P. SUnSET, a nonradioactive method to monitor protein synthesis. *Nature Methods*, 2009, 6(4): 275-277.
7. Goodman CA, Mabrey DM, Frey JW, et al. Novel insights into the regulation of skeletal muscle protein synthesis as revealed by a new nonradioactive in vivo technique. *The FASEB Journal*, 2011, 25(3): 1028-1039.
8. Seglen PO. Preparation of isolated rat liver cells. *Methods in Cell Biology*, 1976, 13: 29-83.
9. Li WC, Ralphs KL, Tosh D. Isolation and culture of adult mouse hepatocytes. *Methods in Molecular Biology*, 2010, 633: 185-196.
10. Blau HM, Pavlath GK, Hardeman EC, et al. Plasticity of the differentiated state. *Science*, 1985, 230(4727): 758-766.
11. Burattini S, Ferri P, Battistelli M, et al. C2C12 murine myoblasts as a model of skeletal muscle development: morpho-functional characterization. *European Journal of Histochemistry*, 2004, 48(3): 223-233.
12. Elbashir SM, Harborth J, Lendeckel W, et al. Duplexes of 21-nucleotide RNAs mediate RNA interference in cultured mammalian cells. *Nature*, 2001, 411(6836): 494-498.
13. Wu G, Bazer FW, Burghardt RC, et al. Proline and hydroxyproline metabolism: implications for animal and human nutrition. *Amino Acids*, 2011, 40(4): 1053-1063.
14. Bustin SA, Benes V, Garson JA, et al. The MIQE guidelines: minimum information for publication of quantitative real-time PCR experiments. *Clinical Chemistry*, 2009, 55(4): 611-622.
15. Bates D, Machler M, Bolker B, Walker S. Fitting linear mixed-effects models using lme4. *Journal of Statistical Software*, 2015, 67(1): 1-48.

---

*注：所有抗体货号和试剂货号的兼容性已根据厂家说明书验证。siRNA 序列需根据 Sus scrofa 11.1 基因组重新设计。实验方案为推荐框架，具体操作应参考各实验室 SOP 进行优化。*
