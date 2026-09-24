# EMP-Web 课程项目报告

会话 ID: xHLj079K91WTfoWpHQbUvH1U
生成时间: 2026-09-24 16:26:50

## 视频测验进度
已通过 5 个步骤测验。

## 科学解读与假设

## 任务反思

- **course_transcriptomics / s1_import**
  基因 ID 呈小鼠基因符号格式（如 0610005C13Rik），导入和差异表达分析时不需要预先转换 ID；后续功能注释可能需要进行 ID 映射。我核对了原始文件：矩阵的 24 个样本列名与 colData 中的 24 个 SampleID 一一对应，顺序也一致

- **course_transcriptomics / s2_prepare**
  DESeq2 用原始整数 count 建立差异分析模型，并在模型内部估计样本间的归一化因子；不能把取对数或其他变换后的数值当作原始 count 输入。低表达基因应按明确的阈值过滤，以减少噪声和多重检验负担。对于潜在批次效应，我会先检查样本信息和 PCA；如果有批次记录且批次与分组不完全重合，就在设计模型中加入批次协变量。若没有可靠的批次信息，则在结果解释中说明这一限制，不声称已经消除了批次效应。

- **course_transcriptomics / s3_analysis**
  本次以 DMSO 为参照组、T4400 为检验组，按 Group 分组进行 DESeq2 差异表达分析。示例样本信息只有 SampleID 和 Group，没有可用于建模的批次字段，因此未加入批次协变量。我采用 padj ≤ 0.05 和 |log2FC| ≥ 1 筛选关注的基因：前者考虑了多重检验，后者表示估计表达变化至少约两倍。这是用于结果筛选和展示的阈值，不等于已经统计证明每个基因的真实变化都超过两倍；后续仍需结合基因功能和实验验证解读。

- **course_transcriptomics / s4_visualization**
  在 DMSO 与 T4400 的差异结果中，我关注 Il6（log2FC = 2.69，FDR = 0.0009）、Ccl3（log2FC = 5.12，FDR 在导出表中显示为 0）和 Cdkn1c（log2FC = −1.25，FDR = 0.0001）。三者均达到 |log2FC| ≥ 1、FDR ≤ 0.05。已知 Il6 参与免疫和炎症反应，Ccl3 与白细胞趋化有关，Cdkn1c 参与细胞周期抑制，因此值得进一步关注。但差异表达不能单独证明细胞来源、通路活性或因果机制，仍需后续验证。

- **course_transcriptomics / s5_interpretation**
  GO-BP 结果中，acute-phase response（11 个差异基因，校正后 p 值约 1.89×10⁻⁹）和 acute inflammatory response（13 个差异基因，校正后 p 值约 4.98×10⁻⁸）显著富集，提示差异基因涉及急性期及炎症相关过程，其中包括 Il6、Il1a 等。人工修正 AI 解读时，我没有把“富集”写成“通路已被激活”或“已证实因果机制”。此外，本次网页 GO 分析采用软件默认注释基因背景，未指定本实验实际检验的基因作为背景，因此显著性需谨慎解释。KEGG 因网络连接失败未得到结果，不能说成“KEGG 无显著通路”。

## Learning Trace 摘要
共 74 条事件。

