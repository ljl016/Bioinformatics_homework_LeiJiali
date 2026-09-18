# AI and Verification Record

**Student:** 雷嘉莉 (`SUAT24000139`)  
**Recorded:** 17 September 2026  
**Scope:** Week 4 synthetic teaching package only  

## Data and execution boundary

- All supplied FASTQ, QC, manifest, variant, coordinate, AF, and ClinVar-style fields are synthetic.
- No real patient data were used, and no synthetic identifier was submitted to a clinical database.
- Q2 includes a local raw-FASTQ integrity/QC inspection. No production alignment, variant calling, mapping rate, BAM, or VCF result is claimed.
- Q4 was executed locally with R 4.6.1 against the supplied 12-row TSV. The original H-drive input was copied temporarily to an ASCII-only path because the local R locale could not open the Chinese path. The source file itself was not modified.

## Q1 AI interaction

**Prompt summary:** Critique an ATAC-first strategy for Gene X upregulation; distinguish what each assay measures from what it cannot prove; convert the assays into a result-dependent decision tree; identify overstatements; finish with a causal validation experiment.

**Correction accepted:** The original statement that accessibility is a prerequisite for TF binding and 3D conformation was too absolute. The final report states that differential accessibility localizes candidate regulatory regions but does not prove TF occupancy, contact, enhancer activity, or causality.

## Q2 plan-first prompt

> I am designing a reproducible workflow for paired-end human whole-genome sequencing (WGS). The available teaching materials include a fictional sample manifest, one tiny synthetic paired-end FASTQ example (S01_CTRL_WGS), and a precomputed FastQC-style summary. A production-scale alignment is not required.
>
> Before suggesting commands or code, produce only a structured analysis plan. For every stage, state: (1) its scientific purpose; (2) the required input; (3) the expected output and file format; (4) the main quality-control checkpoint; (5) one important assumption or limitation; and (6) which official documentation should be checked before implementation.
>
> Cover raw FASTQ inspection and FastQC/MultiQC; at least four supplied QC metrics; adapter and quality trimming followed by repeat QC; reference selection; paired-end alignment; mapped-read processing; small-variant calling and filtering; annotation; visualization; interpretation; and reproducibility records. Treat all inputs as synthetic. Do not claim execution or clinical conclusions. Distinguish required from optional steps, flag parameters requiring verification, and end with a compact workflow and AI-audit checklist. Do not write the final report or implementation code yet.

**Independent check:** I ran `q2_fastq_inspect.py` directly on both supplied gzip-compressed FASTQ files. It counted 120 reads in each mate, found a uniform 80-bp length, and matched all 120 pair identifiers. Exact adapter matches occurred in 18/120 reads in each mate; high-GC reads (GC fraction at least 70%) occurred in 12/120 R1 and 0/120 R2; exact duplicate observations (total reads minus unique sequences) were 34/120 in R1 and 0/120 in R2; the largest identical-sequence group was 18/120 in R1 and 1/120 in R2; and low late-cycle quality (mean Phred at most 15 from cycle 51 onward) occurred in 18/120 reads in each mate. The output is retained as `results/q2_fastq_inspection.tsv`.

**Snapshot comparison:** The direct measurements agree with the planted adapter, R1 high-GC, and R1 duplicate signals stated in `README_data.md` and the supplied snapshot, but expose two simplifications. The snapshot's 15% duplication figure describes the largest repeated template group, not the 28.3% overall R1 exact-duplicate fraction. It also describes late-cycle quality loss only for R1, whereas direct inspection detects the same 15% low-quality subset in R2. These differences are reported rather than silently reconciled. The inspection is not FastQC and does not replace production alignment or variant calling.

## Q3 AI interaction

**Prompt summary:** For a candidate element upstream of Gene Y, classify each ATAC, H3K27ac, methylation, contact, and RNA statement as observation, interpretation, or missing evidence. Generate a rival explanation and a test in which the enhancer hypothesis and rival make different predictions.

**Correction accepted:** No numeric tracks are supplied, so the final answer uses conditional supportive patterns and does not report fictional peaks, methylation percentages, contact strengths, or fold changes.

## Q4 plan-first prompt

> I need to design a reproducible variant-prioritization workflow for a synthetic teaching dataset with the columns CHROM, POS, REF, ALT, FILTER, DP, GQ, AF, GENE, CONSEQUENCE, CLINVAR_SIG, CLINVAR_ID, and NOTE.
>
> Before writing R code, applying filters, counting retained variants, or naming a top variant, produce only a structured analysis plan. Use this pre-specified contract: FILTER must equal PASS; DP must be at least 20; GQ must be at least 30; AF must be at most 0.01; prioritize splice, stop-gained, frameshift, and missense consequences; use clinical annotations only as soft-ranking evidence; repeat with AF at most 0.001; and never treat the synthetic coordinates, IDs, labels, or NOTE field as verified patient evidence. No phenotype, inheritance, segregation, ancestry, or disease-specific gene list is available.
>
> For each step, state whether it is a hard exclusion, soft rank, or post-ranking verification; its rationale; how it could mislead; the authoritative documentation to check; and the reproducibility record to retain. Plan an exclusion log, shortlist, evidence ladder, false-lead audit, prioritization figure, and validation decision. End with pseudocode and an audit checklist. Do not inspect or rank individual rows yet.

**Audit event:** AI initially named TP53 before presenting a reproducible filter. The student rejected that answer and required the rules above to be frozen before row-level ranking.

## Executed Q4 checks

Command concept: `Rscript q4_filter.R variants_q4.tsv results_directory`

- Input rows: 12.
- Candidate rows after locked technical, AF, and consequence filters: 4.
- Excluded rows: 8, with row-specific reasons retained in `results/q4_exclusion_log.tsv`.
- Ranked genes: TP53, BRCA2, KRAS, LDLR.
- AF sensitivity: both 0.01 and 0.001 retained 4 candidates (TP53, BRCA2, LDLR, KRAS).
- Script self-checks passed: 12 input rows, 4 candidates, TP53 ranked first, and identical sensitivity counts.
- Important limitation: the additive score is an explainable classroom ranking aid, not a validated clinical classifier or ACMG implementation.

## Source verification ledger

| ID | Source check | Identifier/status | Claim supported |
|---|---|---|---|
| E01 | PubMed record opened | PMID 24097267; DOI 10.1038/nmeth.2688 | ATAC-seq profiles accessible chromatin/nucleosome organization |
| E02 | PubMed record opened | PMID 31036827; DOI 10.1038/s41467-019-09982-5 | CUT&Tag profiles antibody-targeted chromatin features |
| E03 | PubMed record opened | PMID 27708057; DOI 10.1126/science.aag2445 | CRISPRi can test enhancer-promoter connections |
| E04 | NCBI GRC page opened | GRCh38.p14 official assembly record | Assembly identity and coordinate-version requirement |
| E05 | Babraham official page opened | FastQC official documentation | Modular QC for sequencing data |
| E06 | Official GitHub README opened | BWA-MEM2 paired-end usage | Reference plus separate R1/R2 alignment input |
| E07 | GATK official documentation opened | Germline short-variant workflow | Per-sample GVCF and cohort-genotyping workflow |
| E08 | Ensembl official documentation opened | VEP release 116 | Versioned consequence annotation for GRCh38/VCF inputs |
| E09 | NCBI ClinVar documentation opened | Review-status documentation | Classifications are submissions with differing review levels |
| E10 | PubMed record opened | PMID 25741868; DOI 10.1038/gim.2015.30 | Multiple evidence categories are required for variant interpretation |
| E11 | PubMed record opened | PMID 32461654; DOI 10.1038/s41586-020-2308-7 | Population frequency and annotation QC context from gnomAD |
| E12 | arXiv record opened | arXiv:2609.00065 | Procedural scientific-agent skills used in this workflow |

All sources were checked on 17 September 2026. Search snippets were used only for discovery; identifiers and supported claims were matched to PubMed, arXiv, or official documentation records. The student should still open the final links before submission because human verification remains the final accountability step.

## Limitations and next experiments

- Q1 and Q3 are designs; they contain no collected biological results.
- The sample size for future wet-lab work cannot be fixed without variance/effect-size or precision assumptions. Biological replication, randomization, blocking, and batch balance remain required.
- Q2 validates the supplied pair structure and selected raw-read QC features, but not an alignment or variant call.
- Q4 lacks phenotype, inheritance, segregation, ancestry, zygosity, transcript selection, and real database records.
- The Q4 top variant would require orthogonal DNA confirmation followed by RNA-splicing testing before any mechanistic claim.
