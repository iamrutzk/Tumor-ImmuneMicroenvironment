# Cancer Immunology Sample Extraction Guide

## Overview
This guide provides comprehensive information for extracting and accessing samples related to tumor-infiltrating immune cells across human lung adenocarcinoma, pancreatic cancer, and breast cancer datasets.

---

## 1. HUMAN LUNG ADENOCARCINOMA - Tumor-Infiltrating Immune Cells

### Primary Databases and Accessions

#### 1.1 GEO GSE131907 - Comprehensive scRNA-seq Atlas
- **Title:** Single cell RNA sequencing of lung adenocarcinoma
- **Database:** Gene Expression Omnibus (GEO)
- **Samples:** 58 lung adenocarcinoma samples from 44 patients
- **Total Cells:** 208,506 cells
- **Platform:** Illumina NextSeq 500
- **Coverage:**
  - 11 primary tumors
  - 11 distant normal lung tissues
  - 10 normal lymph nodes
  - 10 brain metastases
  - 7 metastatic lymph nodes
  - 4 lung tumor tissues from advanced stage patients
  - Pleural effusion samples
- **Immune Cells Characterized:** CD8+ T cells, CD4+ T cells, B cells, macrophages, dendritic cells, mast cells, neutrophils
- **Access:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907
- **Download:** RAW data and processed matrices available

#### 1.2 GEO GSE117570 - Early-Stage NSCLC Immune Reprogramming
- **Title:** Reprogramming of Tumor-infiltrating Immune Cells in Early Stage of NSCLC
- **Database:** Gene Expression Omnibus (GEO)
- **Samples:** 8 samples (4 primary tumors + 4 adjacent normal tissues)
- **Total Cells:** 11,485 cells analyzed
- **Patients:** 4 treatment-naive NSCLC patients
- **Platform:** Illumina NextSeq 500
- **Key Features:** 
  - Paired tumor-normal comparison
  - Early-stage disease only
  - Comprehensive immune cell mapping using scRNA-seq
  - Computational deconvolution of 44 TCGA NSCLC + adjacent normal samples
- **Immune Profiling:** CD8+ T cells, NK cells, myeloid cells with differentiation paths
- **GSM Accessions:** GSM3304007-GSM3304014
- **Access:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117570

#### 1.3 GEO GSE123902 - LUAD Heterogeneity
- **Title:** scRNA-seq of primary and metastatic lung adenocarcinoma
- **Database:** Gene Expression Omnibus (GEO)
- **Samples:** 14 samples (primary, metastatic, normal lung tissues)
- **Type:** Single-cell RNA sequencing
- **Coverage:** Inter-patient and intra-tumoral heterogeneity analysis

#### 1.4 TCGA-LUAD - Large-Scale Bulk Data
- **Title:** The Cancer Genome Atlas - Lung Adenocarcinoma
- **Database:** The Cancer Genome Atlas
- **Samples:** 
  - 522 total tumor samples
  - 515 with HiSeq RNA-seq (Gene-level, normalized log2 RPKM)
  - 59 healthy lung tissue samples
- **Data Types Available:**
  - RNA-seq (Gene level and Isoform level)
  - Whole Exome Sequencing (WES)
  - miRNA expression
  - Methylation (HM27 and HM450K)
  - Copy number variations (CNV)
  - Reverse Phase Protein Arrays (RPPA)
  - Clinical data and outcomes
- **Access:** https://www.linkedomics.org/data_download/TCGA-LUAD/
- **Samples:** 533 with somatic mutation calls

#### 1.5 Validation Cohorts
- **GSE50081:** 127 lung carcinoma specimens (bulk RNA-seq)
- **GSE26939, GSE72094:** Additional bulk RNA-seq cohorts for tumor stem cell signatures
- **IMvigor210, GSE78220:** PD-1/PD-L1 blockade immunotherapy response cohorts

### Key Immune Cell Populations Identified in LUAD
- CD8+ T cells (including exhausted phenotypes)
- CD4+ T cells (regulatory and effector subsets)
- Natural Killer (NK) cells
- M1/M2 macrophages with differentiation trajectories
- CD14+ monocytes
- Myeloid-Derived Suppressor Cells (MDSCs)
- Dendritic cells
- B cells and plasma cells
- Mast cells

### Extraction Strategy for LUAD Samples
1. **For comprehensive immune cell atlas:** Use GSE131907 (208k cells from 58 samples)
2. **For tumor-normal paired analysis:** Use GSE117570 (11.5k cells, 4 patients)
3. **For bulk transcriptome validation:** Use TCGA-LUAD (515 samples with RNA-seq)
4. **For immunotherapy prediction:** Use GSE50081, IMvigor210, GSE78220 cohorts
5. **For gene expression associations:** Use TCGA-LUAD with clinical outcomes

---

## 2. HUMAN PANCREATIC CANCER - PDAC Immune Microenvironment

### Primary Databases and Accessions

#### 2.1 GEO GSE212966 - PDAC Immune Landscape
- **Title:** Single-cell RNA-seq reveals immune landscape of pancreatic cancer
- **Database:** Gene Expression Omnibus (GEO)
- **Samples:** 
  - 6 PDAC tumor samples
  - 3 adjacent normal pancreatic tissue samples
- **Type:** Single-cell RNA sequencing
- **Coverage:** Immune cell composition differences between PDAC and normal tissue
- **Access:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212966

#### 2.2 Published Study - Mass Cytometry (Brouwer et al., BMJ JITC 2022)
- **Title:** Local and systemic immune profiles of human pancreatic cancer
- **Patients:** 11 PDAC patients
- **Samples:**
  - PDAC tumor tissues (8 samples)
  - Matched non-malignant pancreatic tissues
  - Regional lymph nodes
  - Spleen samples
  - Portal vein blood samples (pre- and post-surgery)
  - Peripheral blood samples (before and after surgery)
- **Immune Markers:** 41 cell surface and intracellular markers
- **Total Cells Analyzed:** 2+ million immune cells
- **Technology:** Mass cytometry (CyTOF) + flow cytometry
- **Key Finding:** Immunosuppressive landscape with enriched B cells, regulatory T cells, and depletion of cytotoxic CD8+ T cells
- **Clinical Annotation:** Included TMN staging and pretreatment status

#### 2.3 Spatial Transcriptomics Study (GEO GSA)
- **Title:** Spatially resolved multi-omics single-cell analyses inform mechanisms of immune-dysfunction in pancreatic cancer
- **Dataset ID:** STDS0000225
- **Type:** Multi-omics including Spatial Transcriptomics (Visium), scRNA-seq, bulk RNA-seq
- **Tissue Samples:** 8 tissue sections analyzed
- **Coverage:** Spatial immune organization and immune-stromal interactions
- **Access:** https://db.cngb.org/stomics/datasets/STDS0000225

#### 2.4 Published Study - Mass Cytometry (Wang et al., 2022)
- **Title:** A Single-Cell Atlas of Tumor-Infiltrating Immune Cells in Pancreatic Cancer
- **Patients:** 14 patients total
  - 12 PDAC patients (8 with tumor samples for mass cytometry)
  - 2 chronic pancreatitis patients
- **Samples:**
  - 8 treatment-naive PDAC tumor samples
  - 2 chronic pancreatitis (CP) tissues
  - 9 peripheral blood mononuclear cell (PBMC) samples
  - 2 adjacent normal tissues (1 from CP patient, 1 from PDAC patient)
- **Total Cells Analyzed:** 2+ million immune cells from 21 samples (10 patients)
- **Immune Markers:** 33 protein markers
- **Key Cells:** CD8+ T cells, B cells, Tregs, macrophages, dendritic cells, innate lymphoid cells
- **TCGA Validation:** PDAC samples from The Cancer Genome Atlas included

#### 2.5 TCGA-PAAD - Large-Scale Bulk Data
- **Title:** The Cancer Genome Atlas - Pancreatic Adenocarcinoma
- **Database:** The Cancer Genome Atlas / LinkedOmics
- **Samples:**
  - 178 tumor samples with RNA-seq data (HiSeq platform)
  - 185 with complete clinical annotation
  - 178 with miRNA data (Normalized, RPM)
  - 184 with methylation data (CpG and gene level)
  - 126 with somatic mutation data
- **Data Types Available:**
  - RNA-seq (Gene-level, normalized log2 RPKM)
  - miRNA expression (gene and isoform level)
  - Methylation (HM450K platform)
  - Somatic mutations (SNVs)
  - Copy number variations (GISTIC2)
  - Reverse Phase Protein Arrays (RPPA)
  - Clinical data and outcomes
- **Access:** https://www.linkedomics.org/data_download/TCGA-PAAD/
- **Immune Inference:** Multiple deconvolution methods available (CIBERSORT, quanTIseq, etc.)

#### 2.6 Recent Integration Study (Nature Communications, 2025)
- **Title:** Distinct immune cell infiltration patterns in pancreatic ductal adenocarcinoma
- **Authors:** Sivakumar et al.
- **Dataset (PancrImmune):** 
  - 12 treatment-naive PDAC patients
  - Multi-omic profiling: gene expression, ADT-seq, BCR/TCR sequencing
  - Matched tumor tissues and PBMCs
  - Integrated with 2 existing PDAC datasets (Peng and Steele)
- **Key Populations:** B cells, T cells, myeloid cells with detailed sub-classifications

### Key Immune Cell Populations in PDAC
- CD8+ T cells (predominantly exhausted phenotype)
- B cells (both mature and immature)
- Regulatory T cells (Tregs) - INCREASED
- Macrophages (M1 and M2 phenotypes)
- Myeloid-Derived Suppressor Cells (MDSCs)
- Dendritic cells (myeloid and plasmacytoid)
- Innate Lymphoid Cells (ILC1-like populations)
- Natural Killer (NK) cells

### Extraction Strategy for PDAC Samples
1. **For single-cell immune profiling:** Use GSE212966 (9 samples from 6 PDAC + 3 normal)
2. **For comprehensive immune characterization:** Use published mass cytometry studies (Wang et al., Brouwer et al.)
3. **For spatial immune organization:** Use spatial transcriptomics dataset (STDS0000225)
4. **For bulk transcriptome and immune deconvolution:** Use TCGA-PAAD (178+ samples)
5. **For multi-omic integration:** Use recent Nature Communications study (PancrImmune + Peng + Steele datasets)

---

## 3. HUMAN BREAST CANCER - Breast Tumor Immune Landscape

### Primary Databases and Accessions

#### 3.1 GEO GSE161529 - Comprehensive Breast Cancer Single-Cell Atlas
- **Title:** scRNA-seq profiling of breast cancer tumors, BRCA1 mutant pre-neoplastic mammary gland cells and normal mammary gland cells
- **Database:** Gene Expression Omnibus (GEO)
- **Samples:** 69 scRNA-seq profiles from 52 patients
- **Total Cells:** 421,761 cells
- **Platform:** 10X Genomics Chromium
- **Coverage:**
  - 4 Triple Negative Breast Cancer (TNBC)
  - 4 BRCA1-mutant TNBC
  - 6 HER2+ tumors
  - 19 ER+ tumors
  - 6 lymph node metastases from ER+ tumors
  - 13 normal mammary gland samples (hormone-stage dependent)
  - 4 pre-neoplastic BRCA1+/- tissues
- **Immune Cell Diversity:** Characterized across all major immune cell types
- **Access:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529
- **Download:** MTX/TSV formats, 2.2 Gb compressed raw data
- **Note:** Seurat objects available from Chen et al. (Sci Data 2022)

#### 3.2 Integrated Breast Cancer Microenvironment Atlas (bioRxiv 2025)
- **Title:** A Single-Cell Atlas of the Breast Cancer Microenvironment
- **Integration:** 31 single-cell RNA-seq datasets
- **Total Cells:** 1.2 million cells
- **Samples:** 376 samples total
- **Subtypes Covered:**
  - Triple Negative Breast Cancer (TNBC)
  - Luminal A (ER+/PR+, HER2-)
  - Luminal B (ER+/PR+, HER2+)
  - HER2+ subtype
- **Immune Landscape Features:**
  - Subtype-specific immune signatures
  - Stage-specific TME states
  - Checkpoint molecule expression (PD-1, LAG-3, CTLA-4)
  - Distinct tumor-immune interaction programs
- **Clinical Correlation:** Linked to patient survival and therapeutic response

#### 3.3 Published Study - ER+ Breast Cancer Immune Atlas
- **Title:** Immune landscape in estrogen receptor positive breast cancer
- **Type:** Multi-modal analysis
- **Samples:**
  - 50 ER+ invasive ductal carcinoma (IDC) samples
  - 65 ER+ invasive lobular carcinoma (ILC) samples
  - Flow cytometry cohort (subset analyzed)
- **Methods:**
  - Multispectral immunohistochemistry (mIHC)
  - Flow cytometry
  - scRNA-seq
- **Key Findings:** 
  - Macrophages are the predominant immune cells (not T cells)
  - Different immune compositions between IDC and ILC subtypes
  - Macrophage-effector T cell interactions associated with disease-free survival
- **Spatial Analysis:** Multi-region of interest (ROI) based quantification

#### 3.4 TCGA-BRCA - Large-Scale Bulk Data
- **Title:** The Cancer Genome Atlas - Breast Invasive Carcinoma
- **Database:** The Cancer Genome Atlas
- **Samples:** 1,080 breast cancer samples with multi-omics profiles
- **Data Types Available:**
  - RNA-seq (multiple platforms and normalization)
  - Whole Exome Sequencing (WES)
  - Methylation (HM27 and HM450K)
  - Copy number variations (CNV)
  - miRNA expression
  - Reverse Phase Protein Arrays (RPPA)
  - Complete clinical annotation and outcomes
- **Subtypes:** ER+, HER2+, Triple-Negative, and other classifications
- **Immune Profiling:** Six distinct immune subtypes identified with different TME characteristics

#### 3.5 Published Study - Breast Cancer Immunogenomic Landscape
- **Title:** Immunogenomic Landscape in Breast Cancer Reveals Six Immune Subtypes
- **Dataset:** TCGA-BRCA samples (1,080 samples)
- **Analysis:** Integrative immunogenomic analysis
- **Immune Subtypes Identified:**
  1. Immunologically quiet
  2. Chemokine dominant
  3. Lymphocyte depleted
  4. Wounding dominant
  5. Innate immune dominant
  6. IFN-γ dominant
- **Clinical Validation:** 4 cohorts with immune checkpoint inhibitor treatment outcomes

#### 3.6 External Validation Cohort
- **GSE20711:** Microarray-based breast cancer transcriptome validation cohort
- **Use:** External validation for prognostic immune-related signatures

### Key Immune Cell Populations in Breast Cancer
- Macrophages (M1 and M2, with cycling M2 in ER+ tumors)
- CD8+ T cells (highly proliferative in TNBC/HER2+)
- CD4+ T cells (regulatory and effector subsets)
- B cells
- Plasma cells
- Dendritic cells
- Regulatory T cells (Tregs)
- Natural Killer (NK) cells
- Mast cells

### Subtype-Specific Immune Signatures
- **TNBC:** High proliferative T cells, activated inflammatory responses, elevated immune infiltration
- **HER2+:** Similar to TNBC with cycling T cells, good ICI response potential
- **ER+/Luminal A:** Cycling macrophages, immune exclusion, less T cell infiltration, lower immune infiltration
- **ER+/Luminal B:** Mixed immune signatures with higher proliferation

### Extraction Strategy for Breast Cancer Samples
1. **For single-cell immune landscape:** Use GSE161529 (421k cells from 69 profiles, 52 patients)
2. **For integrated multi-dataset analysis:** Use the 31-dataset integration study (1.2M cells)
3. **For ER+ immune microenvironment:** Use the published ER+ breast cancer immune landscape study
4. **For bulk transcriptome and immune deconvolution:** Use TCGA-BRCA (1,080 samples)
5. **For subtype-specific immune signatures:** Use TCGA-BRCA immunogenomic classification study
6. **For external validation:** Use GSE20711 microarray cohort

---

## 4. HUMAN LUNG ADENOCARCINOMA - Bulk Tumor Validation

### Primary Data Sources

#### 4.1 TCGA-LUAD - Comprehensive Bulk RNA-seq
- **Samples:** 515 primary tumors with HiSeq RNA-seq data
- **Data Format:** Gene-level, normalized log2 RPKM
- **Total Genes:** 19,988 genes quantified
- **Associated Data:** 
  - Somatic mutations (533 samples)
  - Copy number variations (516 samples)
  - Methylation data (458 samples)
  - Clinical outcomes and survival data (522 samples)
  - Patient demographics and staging information
- **Access:** https://www.linkedomics.org/data_download/TCGA-LUAD/

#### 4.2 GEO GSE50081 - LUAD Validation Cohort
- **Title:** LUAD validation dataset
- **Samples:** 127 lung carcinoma specimens
- **Type:** Bulk RNA-seq or microarray
- **Use:** Validation of prognostic signatures derived from single-cell data

#### 4.3 Additional Bulk Cohorts
- **GSE26939:** Bulk RNA-seq LUAD cohort (training dataset)
- **GSE72094:** Bulk RNA-seq LUAD cohort (additional validation)
- **Total Coverage:** Multiple cohorts for independent validation

#### 4.4 Immunotherapy Response Cohorts
- **IMvigor210:** Anti-PD-1 (nivolumab) treated NSCLC/LUAD patients with bulk RNA-seq
- **GSE78220:** Additional immunotherapy response cohort with matched RNA-seq and clinical outcomes
- **Key Metric:** Response/progression-free survival following immune checkpoint blockade

### Validation Strategy
1. **Primary Model Development:** TCGA-LUAD (515 samples with outcomes)
2. **Cross-validation:** GSE50081 (127 samples)
3. **Additional Validation:** GSE26939, GSE72094
4. **Clinical Utility Testing:** IMvigor210, GSE78220 (immunotherapy prediction)
5. **Immune Infiltration Prediction:** TIMER2.0 or similar algorithms on bulk RNA-seq

---

## 5. HEALTHY vs NSCLC - Immune Cell Comparison

### Paired Normal-Tumor Datasets

#### 5.1 GEO GSE117570 - Direct Paired Comparison
- **Tumor Samples:** 4 primary NSCLCs (treatment-naive)
- **Normal Samples:** 4 adjacent normal lung tissues (matched pairs)
- **Total Cells:** 11,485 cells analyzed
- **Immune Profiling:** 
  - CD8+ T cells (decreased in tumor)
  - NK cells (decreased in tumor)
  - M1/M2 macrophages (marked reprogramming)
  - Myeloid cell differentiation trajectories
- **Key Finding:** Dynamic immune reprogramming from normal to tumor

#### 5.2 TCGA-LUAD Computational Comparison
- **Tumor Samples:** 44 NSCLC tumors
- **Normal Samples:** 44 adjacent normal lung tissues (matched)
- **Method:** Computational deconvolution of bulk RNA-seq
- **Algorithms:** CIBERSORT, TIMER, and other immune inference methods
- **Outcome:** Heterogeneous patterns of immune cell alterations between tumor and normal

#### 5.3 GEO GSE131907 - LUAD vs Normal Lung
- **Normal Lung:** 10 distant normal lung tissues
- **Normal Lymph Node:** 10 normal lymph nodes
- **Primary Tumors:** 11 primary tumors
- **Metastatic Samples:** Brain, lymph node, and pleural metastases
- **Total Cells:** 208,506 cells
- **Comprehensive Comparison:** Immune composition across tumor microenvironment progression

#### 5.4 Published Study - Immune Cell Characterization
- **Title:** Immune biology of NSCLC revealed by single-cell proteomic analysis
- **Method:** Flow cytometry with 35 immune markers
- **Immune Populations Characterized:**
  - CD8+ T cells
  - CD4+ T cells
  - Regulatory T cells (Tregs)
  - Natural Killer (NK) cells
  - NKT cells
  - B cells
  - Granulocytes
  - Macrophages
  - Dendritic cells
- **Clinical Relevance:** Predictive of immunotherapy response (PD-1 blockade)

### Key Comparisons in Healthy vs NSCLC
1. **T Cell Infiltration:** Generally decreased CD8+ and CD4+ in tumor
2. **NK Cell Abundance:** Reduced in tumor microenvironment
3. **Myeloid Differentiation:** Distinct M1-to-M2 transition in tumors
4. **Immune Checkpoint Expression:** Increased PD-1, TIM-3, LAG-3 in tumor-infiltrating T cells
5. **Functional Status:** Exhausted vs activated immune phenotypes

### Extraction Strategy
1. **For direct single-cell comparison:** Use GSE117570 (4 tumor-normal pairs)
2. **For large-scale computational analysis:** Use TCGA-LUAD (44 tumor-normal pairs)
3. **For comprehensive immune atlas:** Use GSE131907 (normal + tumor tissues)
4. **For immune marker profiling:** Use published proteomic analysis study

---

## 6. Data Access and Download Methods

### GEO Access
1. Navigate to: https://www.ncbi.nlm.nih.gov/geo/
2. Enter Series accession (e.g., GSE131907)
3. Download options:
   - Series matrix files (TXT)
   - RAW data (TAR compressed)
   - Individual sample data (GSM accessions)
   - Supplementary files from authors

### TCGA Access
1. **Via Genomic Data Commons (GDC):** https://portal.gdc.cancer.gov/
2. **Via LinkedOmics:** https://www.linkedomics.org/
3. **Download formats:** Multiple normalized levels available
4. **Clinical data:** Associated with molecular data through patient IDs

### Bioinformatics Tools for Analysis
- **CancerSCEM:** Cancer Single-cell Expression Map database
- **TIMER2.0:** Tumor infiltrating immune estimation
- **CellChat:** Cell-cell interaction analysis
- **Seurat:** R package for single-cell analysis
- **scVelo:** RNA velocity analysis

---

## 7. Key Publications and Citations

### Lung Adenocarcinoma
- GSE131907: scRNA-seq comprehensive atlas of LUAD
- GSE117570: Early-stage immune reprogramming
- TIMER2.0: Immune deconvolution methods

### Pancreatic Cancer
- Wang et al. (2022): Mass cytometry atlas of PDAC immune cells
- Brouwer et al. (BMJ JITC 2022): Local and systemic immune profiles
- Recent studies on spatial transcriptomics and multi-omics

### Breast Cancer
- GSE161529: Comprehensive breast cancer scRNA-seq atlas
- Chen et al. (Scientific Data 2022): Data descriptor and analysis
- Multiple studies on subtype-specific immune signatures

---

## 8. Sample Size Summary Table

| Cancer Type | Study Type | Database | Samples | Cells | Reference |
|---|---|---|---|---|---|
| LUAD | scRNA-seq | GEO GSE131907 | 58 (44 patients) | 208,506 | [1] |
| LUAD | scRNA-seq | GEO GSE117570 | 8 (4 patients) | 11,485 | [2] |
| LUAD | Bulk RNA-seq | TCGA-LUAD | 515 | - | [3] |
| PDAC | scRNA-seq | GEO GSE212966 | 9 (6 PDAC + 3 normal) | - | [4] |
| PDAC | Mass cytometry | Published | 12 patients | 2M+ cells | [5] |
| PDAC | Bulk RNA-seq | TCGA-PAAD | 178 | - | [6] |
| Breast | scRNA-seq | GEO GSE161529 | 69 (52 patients) | 421,761 | [7] |
| Breast | Integration | Published | 376 samples | 1.2M | [8] |
| Breast | Bulk RNA-seq | TCGA-BRCA | 1,080 | - | [9] |
| NSCLC | Paired scRNA-seq | GEO GSE117570 | 8 | 11,485 | [2] |
| NSCLC | Deconvolution | TCGA-LUAD | 44 tumor + 44 normal | - | [3] |

---

## 9. Quality Considerations for Sample Selection

### Single-Cell RNA-seq Quality Metrics
- Cell number per sample (aim for 1,000+ cells)
- Detection rate (percentage of detected genes per cell)
- Mitochondrial gene percentage (<20% considered good)
- Batch effects between samples
- Cell cycle phase distribution

### Bulk RNA-seq Quality Metrics
- Read depth and coverage uniformity
- RNA integrity number (RIN) where available
- Alignment rates to reference genome
- Gene expression consistency
- Clinical metadata completeness

### Immune Cell Purity
- CD45+ enrichment percentage
- Major immune compartment representation
- Absence of major batch effects
- Platform-specific considerations

---

## 10. Recommended Workflows

### For Immune Cell Atlas Development
1. Primary source: GSE131907 (LUAD) or GSE161529 (Breast)
2. Comprehensive integration with multiple datasets
3. Cell type annotation using established markers
4. Immune state classification and functional assessment
5. Validation in independent cohorts

### For Immune Signature Validation
1. Derive signatures from scRNA-seq data
2. Validate using bulk RNA-seq (TCGA cohorts)
3. Cross-validate across cancer types
4. Test predictive value in immunotherapy cohorts
5. Perform functional validation

### For Tumor-Normal Comparison Studies
1. Use paired samples when available (GSE117570)
2. Computational deconvolution on TCGA (GSE131907)
3. Direct scRNA-seq comparison across conditions
4. Pathway enrichment and functional analysis
5. Clinical outcome correlation

---

## Conclusion

This guide provides a comprehensive resource for extracting, accessing, and utilizing tumor-infiltrating immune cell samples across major human cancers. The datasets span from single-cell level analysis to bulk transcriptomics, offering multiple avenues for immunology research and therapeutic target discovery. Researchers should select datasets based on specific research questions, sample size requirements, and analytical approaches needed.

For latest updates and new datasets, consult CancerSCEM database and GEO updates regularly.
