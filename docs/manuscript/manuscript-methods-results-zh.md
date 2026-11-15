# InterMeta TCGA-PRAD 方法学与结果表述模板

本文件与 `inst/examples/example09_tcga_prad_complete_workflow.R` 的固定分析参数对应。
它是论文写作底稿，不替代作者对最终运行结果、引用和统计解释的人工复核。

## 建议的方法学表述

TCGA-PRAD 的标准化 bulk transcriptomic 数据和生化复发
（biochemical recurrence，BCR）随访数据来源于 PCaDB。表达矩阵和临床资料按
样本标识符对齐后，共纳入 463 例患者，其中 99 例发生 BCR，364 例为删失；
随访时间单位为天。代谢基因注释采用 InterMeta 中的 `geneset-all` 基因集，
共覆盖 52 条代谢通路。正式稿应在此处加入 PCaDB 的文献引用和数据库访问日期。

PRAD 示例采用 InterMeta 的 `TimiGP-legacy` 兼容分析模式，以保持与历史 TimiGP
分析流程一致。表达数据首先限定于代谢基因集覆盖的基因。在兼容模式下，仅保留
至少半数样本中表达值大于 0 的基因；为避免对 `log10(x + 1)` 变换产生非有限值，
进一步去除任一样本表达值小于或等于 -1 的基因，随后进行 `log10(x + 1)` 转换，
并按基因减去跨样本中位数。全队列预处理后保留 861 个基因。

对任意两个保留基因构建样本内相对表达关系。仅保留两种相对表达状态均至少出现
在 `round(0.1 × n)` 个样本中的基因对。对每个基因对分别建立单因素 Cox 比例风险
模型，结局由 BCR 状态和随访时间定义。对于风险比大于 1 的基因对，交换基因顺序
并取风险比倒数，使输出方向统一表示由预后较有利的一侧指向预后较不利的一侧。
Cox P 值同时采用 Benjamini-Hochberg 方法计算校正值，但本次历史复现分析按预先
设定的未校正 `PV < 0.05` 筛选候选基因对，因此该部分按探索性分析解释。

将候选基因对映射至全部有向通路组合，并采用超几何检验评估候选基因对在各通路
相互作用中的富集程度；通路相互作用 P 值采用 Benjamini-Hochberg 方法校正。另
进行 100 次通路基因注释置换：每次在保持通路注释结构的条件下无放回打乱基因
标签，重新计算全部通路相互作用的富集 P 值，并据此估计 permutation FDR。
兼容模式使用固定置换序列 1–100。Permutation FDR 为补充评价指标；通路评分仍按
`Adjust.P.Value < 0.05` 的通路相互作用计算。

对达到上述通路富集阈值的有向相互作用，分别计算每条通路作为有利端和不利端的
相对出现比例，定义 `Difference = Favorable.Score - Unfavorable.Score`。
`Difference > 0` 表示偏向有利预后方向，`Difference < 0` 表示偏向不利预后方向；
兼容模式不保留 `Difference = 0` 的中性通路。

为评价分析流程对样本扰动的内部一致性，在全队列之外预先进行 3 次内部子抽样。
每次采用固定随机种子 20260701、20260702 和 20260703，从全队列中无放回抽取
80% 的患者；该历史复现流程不按事件状态分层。每个子样本均独立重复预处理、
基因对生成、Cox 分析、`PV < 0.05` 筛选、通路富集、100 次注释置换和通路评分。
汇总指标包括通路入选次数与频率、方向一致率、`Difference` 的均值和标准差、
Top 5/10/20 重现频率、两两评分相关性及 Top-K Jaccard 重叠度。由于仅进行了
3 次抽样，这些指标仅作为预先设定的描述性内部一致性检查，不作为置信区间、
正式稳健性检验或外部验证。

特别需要区分：三次内部子抽样与每次分析中的 100 次注释置换是两个不同层次的
操作。前者描述三个预设子样本间的内部一致性，后者用于估计通路相互作用的
permutation FDR。

## 1.1.3 冻结参数实跑的结果表述草稿

以下数值来自 2026-09-17 的无检查点复用运行
`run_20260917_101835_fresh/InterMeta_TCGA_PRAD_complete_1_1_3`。正式稿仍应结合
论文图表编号和生物学解释进行人工复核。

> 在 463 例 TCGA-PRAD 患者中，legacy 预处理保留了 861 个代谢基因并生成
> 370,230 个基因对，其中 101,323 个基因对满足探索性筛选标准 `PV < 0.05`。
> 52 条代谢通路形成 2,652 个有向通路相互作用，其中 617 个满足 BH 校正 P 值
> 小于 0.05；100 次注释置换后，104 个和 50 个相互作用分别满足 permutation
> FDR 小于 0.20 和 0.05。以 `Adjust.P.Value < 0.05` 计算通路评分后，共得到
> 50 条非中性通路。

> 全队列中 Difference 最高的 5 条有利方向通路依次为 beta-alanine metabolism
> （7.455）、propanoate metabolism（6.645）、valine, leucine and isoleucine
> degradation（5.835）、butanoate metabolism（5.348）和 histidine metabolism
> （5.024）。Difference 最低的 5 条不利方向通路依次为 taurine and
> hypotaurine metabolism（-7.131）、oxidative phosphorylation（-6.969）、
> inositol phosphate metabolism（-6.321）、glycerophospholipid metabolism
> （-5.511）和 lysine degradation（-4.700）。通路英文名称应在正式稿中与图表
> 采用同一套大小写和连字符规范。

> 三次预设内部子抽样均成功，每次包含 370 例患者，分别包含 82、84 和 79 个
> BCR 事件，筛得 106,015、87,284 和 101,723 个候选基因对，并得到 47、50 和
> 51 条非中性评分通路。三个子样本间的两两全通路 Pearson 相关系数为
> 0.902–0.951，Spearman 相关系数为 0.895–0.919，通路方向一致率为
> 0.872–0.900，总体 Top-5 通路集合的两两 Jaccard 指数为 0.429–0.667。
> 这些结果描述三个预设子样本之间的内部一致性，不构成独立队列验证，也不用于
> 估计精确的统计不确定性。

## 最终数字的指定来源

- `01_full_sample/tables/parity_checks.tsv`：结构与参考一致性检查；
- `01_full_sample/tables/selected_gene_pairs.tsv`：全队列候选基因对；
- `01_full_sample/tables/enrichment.tsv`：BH 和 permutation FDR 计数；
- `01_full_sample/tables/pathway_score.tsv`：通路数、方向与排序；
- `02_internal_resampling/run_status.txt`：各轮样本、事件、候选对数和状态；
- `02_internal_resampling/stability_summary.txt`：入选与方向一致性；
- `02_internal_resampling/stability_pairwise.txt`：相关系数与 Jaccard 指数；
- `00_metadata/run_settings.txt`、`sessionInfo.txt` 和
  `output_manifest_sha256.tsv`：参数、环境和文件校验。

不要把三次抽样下代码输出的 `Stability.Class = High` 直接改写为“高度稳定”。
正文应优先报告明确的 `3/3`、方向一致率、相关系数和 Jaccard 数值，并称为
“内部一致性”。
