# Complete Pipeline: Downloading & Running Neutrophil/NETosis scRNA-seq Project on Mac

## Table of Contents
1. [Environment Setup](#environment-setup)
2. [Installing Required Tools](#installing-required-tools)
3. [Downloading Raw Data](#downloading-raw-data)
4. [Data Processing Pipeline](#data-processing-pipeline)
5. [Analysis Workflow](#analysis-workflow)
6. [Troubleshooting](#troubleshooting)

---

## ENVIRONMENT SETUP

### Step 1: Create Project Directory Structure

```bash
# Create main project directory
mkdir -p ~/neutrophil_netosis_project
cd ~/neutrophil_netosis_project

# Create subdirectories for organization
mkdir -p data/{raw,processed,references}
mkdir -p results/{figures,tables,metrics}
mkdir -p code/{preprocessing,analysis,visualization}
mkdir -p logs

# Verify structure
tree -L 2
```

### Step 2: Check Mac System Requirements

```bash
# Check Mac OS version
sw_vers

# Check available disk space (need ~500GB for all data)
df -h

# Check available RAM
system_profiler SPHardwareDataType

# Verify you have Xcode Command Line Tools
xcode-select --install

# Check if Homebrew is installed
brew --version

# If not installed, install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### Step 3: Configure Shell Environment

```bash
# Open or create ~/.zprofile (for Mac with Apple Silicon or newer Macs)
nano ~/.zprofile

# Add these lines at the end:
export PATH="/usr/local/bin:$PATH"
export PATH="$HOME/.local/bin:$PATH"
export JAVA_HOME=$(/usr/libexec/java_home)

# Save (Ctrl+X, Y, Enter)

# Reload profile
source ~/.zprofile

# Verify
echo $PATH
```

---

## INSTALLING REQUIRED TOOLS

### Step 1: Install Conda (Recommended for Mac)

```bash
# Download Miniconda for Mac (recommended over Anaconda for space)
# For Apple Silicon (M1/M2/M3):
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh

# For Intel Mac:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

# Install Miniconda
bash Miniconda3-latest-MacOSX-*.sh

# Follow prompts (press Enter, type 'yes', press Enter)

# Initialize conda
conda init zsh

# Close and reopen terminal, then verify
conda --version
```

### Step 2: Create Conda Environments

```bash
# Environment for scRNA-seq Python tools
conda create -n scrna-seq python=3.11 -y

# Environment for R/Seurat analysis
conda create -n seurat-r r-base=4.3 r-essentials -y

# Environment for bioinformatics tools
conda create -n bioinfo python=3.11 -y

# Activate scrna-seq environment
conda activate scrna-seq
```

### Step 3: Install Python-Based Tools (in scrna-seq environment)

```bash
# Make sure you're in scrna-seq environment
conda activate scrna-seq

# Install core bioinformatics packages
conda install -c bioconda -c conda-forge -y \
  scanpy \
  anndata \
  numpy \
  pandas \
  scipy \
  scikit-learn \
  matplotlib \
  seaborn \
  jupyter \
  notebook \
  ipython

# Install cell type annotation tools
pip install celltist
pip install sctype

# Install specialized tools
pip install cellchat  # Note: cellchat R wrapper, may need manual setup
pip install scvelo  # RNA velocity
pip install palantir  # Trajectory inference
pip install pymde  # Visualization

# Install data download tools
pip install SRA-Toolkit
pip install boto3  # AWS S3 access
pip install google-cloud-storage  # Google Cloud access

# Verify installations
python -c "import scanpy; print(scanpy.__version__)"
python -c "import pandas; print(pandas.__version__)"
```

### Step 4: Install R and Seurat (in seurat-r environment)

```bash
# Activate R environment
conda activate seurat-r

# Install Seurat
conda install -c bioconda -c conda-forge -y \
  r-seurat

# Install additional R packages via conda
conda install -c conda-forge -y \
  r-ggplot2 \
  r-dplyr \
  r-patchwork \
  r-harmony

# Start R and install from CRAN/Bioconductor
R

# In R console, install additional packages:
# install.packages("devtools")
# install.packages("clusterprofiler")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("pathview")

# Exit R
# q()
# n
```

### Step 5: Install Cell Ranger (Optional - for raw FASTQ processing)

```bash
# Cell Ranger requires registration but is free
# Download from: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

# Assuming downloaded to ~/Downloads
cd ~/Downloads
tar -xzf cellranger-7.1.0.tar.gz

# Move to appropriate location
mv cellranger-7.1.0 /usr/local/bin/

# Add to PATH in ~/.zprofile
export PATH="/usr/local/bin/cellranger-7.1.0:$PATH"

# Reload and verify
source ~/.zprofile
cellranger --version
```

### Step 6: Install Additional Mac-Specific Tools

```bash
# Using Homebrew for some tools
brew install git
brew install wget
brew install curl
brew install ncbi-blast+  # For sequence analysis if needed
brew install samtools  # For BAM file handling
brew install bedtools   # For genomic analysis

# Verify installations
git --version
wget --version
samtools --version
```

---

## DOWNLOADING RAW DATA

### Step 1: Set Up GEO Utilities

```bash
# Install GEO utilities (for downloading from NCBI GEO)
conda activate scrna-seq

pip install GEOparse

# Create script for batch GEO downloads
cd ~/neutrophil_netosis_project/code
cat > download_geo_data.py << 'EOF'
#!/usr/bin/env python3

import os
import subprocess
import gzip
import shutil
from pathlib import Path

# Define datasets
datasets = {
    "GSE188288": "Circulating_Neutrophils",
    "GSE131907": "LUAD_Tumor_Immune",
}

# Base directory
base_dir = Path.home() / "neutrophil_netosis_project" / "data" / "raw"
base_dir.mkdir(parents=True, exist_ok=True)

def download_geo_dataset(accession, sample_name):
    """Download dataset from GEO"""
    output_dir = base_dir / accession
    output_dir.mkdir(exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"Downloading {accession} ({sample_name})...")
    print(f"{'='*60}\n")
    
    # Method 1: Using wget
    url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession}&format=file"
    
    cmd = [
        "wget",
        "-r",  # recursive download
        "-np",  # no parent
        "-nH",  # no host directories
        "-P", str(output_dir),
        url
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"✓ Downloaded {accession}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error downloading {accession}: {e}")

def extract_gz_files(directory):
    """Extract all .gz files"""
    for gz_file in Path(directory).glob("**/*.gz"):
        if str(gz_file).endswith(".tar.gz"):
            continue  # Skip tar.gz for now
        
        print(f"Extracting {gz_file.name}...")
        with gzip.open(gz_file, 'rb') as f_in:
            output_file = gz_file.with_suffix('')
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # Optional: remove .gz after extraction
        # gz_file.unlink()

# Main execution
for accession, name in datasets.items():
    download_geo_dataset(accession, name)
    extract_gz_files(base_dir / accession)

print("\n" + "="*60)
print("Download complete!")
print("="*60)
EOF

chmod +x download_geo_data.py
```

### Step 2: Download GSE188288 (Circulating Neutrophils - PRIORITY 1)

```bash
cd ~/neutrophil_netosis_project/data/raw

# Method 1: Using curl and wget
mkdir GSE188288
cd GSE188288

# Download the series matrix file (gene expression summary)
echo "Downloading GSE188288 series matrix..."
wget -q "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188288&format=file" -O GSE188288_series_matrix.txt.gz

# Download supplementary files (processed data)
echo "Downloading supplementary files..."
wget -r -np -nH -P . \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188288&format=file&file=GSE188288"

# Extract files
gunzip *.gz

# List downloaded files
ls -lah

echo "✓ GSE188288 download complete"
```

### Step 3: Download GSE131907 (LUAD Tumor Immune - PRIORITY 2)

```bash
cd ~/neutrophil_netosis_project/data/raw

mkdir GSE131907
cd GSE131907

echo "Downloading GSE131907 series matrix..."
wget -q "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131907&format=file" -O GSE131907_series_matrix.txt.gz

echo "Downloading supplementary files..."
wget -r -np -nH -P . \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131907&format=file&file=GSE131907"

# Extract
gunzip *.gz

echo "✓ GSE131907 download complete"
```

### Step 4: Download TCGA-LUAD (Bulk Validation Data)

```bash
cd ~/neutrophil_netosis_project/data/raw

mkdir TCGA-LUAD
cd TCGA-LUAD

# Method 1: Using GDC API (recommended)
echo "Downloading TCGA-LUAD RNA-seq data..."

# Create query for RNA-seq files
cat > query.json << 'EOF'
{
  "filters": {
    "op": "and",
    "content": [
      {
        "op": "in",
        "content": {
          "field": "cases.project.project_id",
          "value": ["TCGA-LUAD"]
        }
      },
      {
        "op": "in",
        "content": {
          "field": "data_type",
          "value": ["Gene Expression Quantification"]
        }
      },
      {
        "op": "in",
        "content": {
          "field": "analysis.workflow_type",
          "value": ["STAR - Counts"]
        }
      }
    ]
  },
  "format": "JSON",
  "size": 1000
}
EOF

# Query GDC API
curl --request POST --header "Content-Type: application/json" \
  --data @query.json \
  "https://api.gdc.cancer.gov/files" > gdc_file_list.json

# Extract file IDs
cat gdc_file_list.json | \
  python3 -c "import sys, json; data = json.load(sys.stdin); print('\n'.join([f['id'] for f in data['data']['hits']]))" \
  > file_ids.txt

# Download first 50 files as sample
head -50 file_ids.txt | while read file_id; do
  echo "Downloading $file_id..."
  curl --remote-name --remote-time \
    "https://api.gdc.cancer.gov/data/$file_id"
done

echo "✓ TCGA-LUAD download complete"
```

### Step 5: Download Reference Genome (GRCh38 - Optional, for alignment)

```bash
cd ~/neutrophil_netosis_project/data/references

echo "Downloading GRCh38 reference genome..."

# Download reference
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.p13.genome.fa.gz

# Download annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz

# Extract
gunzip *.gz

echo "✓ Reference genome downloaded"
```

### Step 6: Create Download Summary Script

```bash
cd ~/neutrophil_netosis_project/code

cat > check_downloads.py << 'EOF'
#!/usr/bin/env python3

import os
from pathlib import Path
import json

base_dir = Path.home() / "neutrophil_netosis_project" / "data" / "raw"

datasets = {
    "GSE188288": "Circulating Neutrophils",
    "GSE131907": "LUAD Tumor Immune",
    "TCGA-LUAD": "Bulk RNA-seq Validation",
}

print("\n" + "="*70)
print("DATA DOWNLOAD STATUS REPORT")
print("="*70 + "\n")

total_size = 0
for dataset, description in datasets.items():
    dataset_path = base_dir / dataset
    
    if dataset_path.exists():
        # Calculate total size
        size = sum(f.stat().st_size for f in dataset_path.rglob('*') if f.is_file())
        size_gb = size / (1024**3)
        total_size += size
        
        # Count files
        file_count = len(list(dataset_path.rglob('*')))
        
        print(f"✓ {dataset}: {description}")
        print(f"  Location: {dataset_path}")
        print(f"  Size: {size_gb:.2f} GB")
        print(f"  Files: {file_count}")
    else:
        print(f"✗ {dataset}: NOT FOUND")
    
    print()

print(f"{'='*70}")
print(f"Total Size: {total_size/(1024**3):.2f} GB")
print(f"{'='*70}\n")
EOF

python3 check_downloads.py
```

---

## DATA PROCESSING PIPELINE

### Step 1: Quality Control of Raw Data

```bash
# Create QC script
cd ~/neutrophil_netosis_project/code/preprocessing

cat > 01_quality_control.py << 'EOF'
#!/usr/bin/env python3

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# Set style
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')

# Define paths
base_dir = Path.home() / "neutrophil_netosis_project"
data_dir = base_dir / "data" / "raw"
results_dir = base_dir / "results" / "metrics"
results_dir.mkdir(parents=True, exist_ok=True)

# Load GSE188288 (example)
print("Loading GSE188288 data...")
adata = sc.read_h5ad(data_dir / "GSE188288" / "GSE188288_processed.h5ad")

print(f"\nData shape: {adata.shape}")
print(f"Observations: {adata.n_obs}")
print(f"Variables: {adata.n_vars}")

# Calculate QC metrics
print("\nCalculating QC metrics...")
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Create QC plots
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Plot 1: nUMI distribution
axes[0].hist(adata.obs['total_counts'], bins=100)
axes[0].set_xlabel('Total Counts (nUMI)')
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('UMI Distribution')

# Plot 2: nGene distribution
axes[1].hist(adata.obs['n_genes_by_counts'], bins=100)
axes[1].set_xlabel('Number of Genes')
axes[1].set_ylabel('Number of Cells')
axes[1].set_title('Gene Count Distribution')

# Plot 3: Mitochondrial percentage
axes[2].hist(adata.obs['pct_counts_mt'], bins=100)
axes[2].set_xlabel('% Mitochondrial Genes')
axes[2].set_ylabel('Number of Cells')
axes[2].set_title('Mitochondrial Gene Percentage')

plt.tight_layout()
plt.savefig(results_dir / 'GSE188288_QC_metrics.pdf')
print(f"✓ QC plots saved to {results_dir / 'GSE188288_QC_metrics.pdf'}")

# Save summary statistics
summary_stats = {
    'total_cells': adata.n_obs,
    'total_genes': adata.n_vars,
    'median_umis': adata.obs['total_counts'].median(),
    'median_genes': adata.obs['n_genes_by_counts'].median(),
    'median_mt_pct': adata.obs['pct_counts_mt'].median(),
}

print("\nQC Summary Statistics:")
for key, value in summary_stats.items():
    print(f"  {key}: {value:.2f}")

# Save to CSV
pd.DataFrame([summary_stats]).to_csv(results_dir / 'GSE188288_summary_stats.csv')

EOF

python3 01_quality_control.py
```

### Step 2: Cell Filtering and Normalization

```bash
cat > 02_preprocessing.py << 'EOF'
#!/usr/bin/env python3

import scanpy as sc
import numpy as np
from pathlib import Path

# Settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100)

# Paths
base_dir = Path.home() / "neutrophil_netosis_project"
data_dir = base_dir / "data" / "raw"
processed_dir = base_dir / "data" / "processed"
processed_dir.mkdir(parents=True, exist_ok=True)

# Load data
print("Loading GSE188288 data...")
adata = sc.read_h5ad(data_dir / "GSE188288" / "GSE188288_processed.h5ad")

print(f"\nOriginal data shape: {adata.shape}")

# Step 1: Filter cells based on QC metrics
print("\nApplying cell filters...")
# For neutrophils: relaxed thresholds due to lower mRNA content
min_genes = 200
max_genes = 10000
max_mt_pct = 25

adata = adata[
    (adata.obs['n_genes_by_counts'] > min_genes) &
    (adata.obs['n_genes_by_counts'] < max_genes) &
    (adata.obs['pct_counts_mt'] < max_mt_pct)
]

print(f"After filtering: {adata.shape}")

# Step 2: Filter genes
print("\nApplying gene filters...")
# Keep genes expressed in at least 3 cells
sc.pp.filter_genes(adata, min_cells=3)
print(f"After gene filtering: {adata.shape}")

# Step 3: Normalize
print("\nNormalizing data...")
sc.pp.normalize_total(adata, target_sum=1e4)
adata.layers['counts'] = adata.X.copy()

# Log transformation
sc.pp.log1p(adata)

# Step 4: Identify highly variable genes
print("\nIdentifying highly variable genes...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=None)

# Step 5: Scale data
print("\nScaling data...")
sc.pp.scale(adata, max_value=10)

# Step 6: Run PCA
print("\nRunning PCA...")
sc.tl.pca(adata, n_comps=50)

# Step 7: Compute nearest neighbors
print("\nComputing nearest neighbors...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# Step 8: UMAP
print("\nRunning UMAP...")
sc.tl.umap(adata)

# Save processed data
output_file = processed_dir / "GSE188288_processed.h5ad"
adata.write(output_file)
print(f"\n✓ Processed data saved to {output_file}")

# Print final summary
print(f"\nFinal data shape: {adata.shape}")
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

EOF

python3 02_preprocessing.py
```

### Step 3: Cell Type Annotation

```bash
cat > 03_cell_annotation.py << 'EOF'
#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
from pathlib import Path

sc.settings.verbosity = 3

# Paths
base_dir = Path.home() / "neutrophil_netosis_project"
processed_dir = base_dir / "data" / "processed"

# Load processed data
print("Loading processed data...")
adata = sc.read_h5ad(processed_dir / "GSE188288_processed.h5ad")

# Define neutrophil markers
neutrophil_markers = {
    'FCGR3B': 'CD16',  # CD16+ neutrophils
    'CXCR2': 'CXCR2',  # Chemokine receptor
    'S100A8': 'S100A8',
    'S100A9': 'S100A9',
    'ELANE': 'NE',  # Neutrophil elastase
    'MPO': 'MPO',  # Myeloperoxidase
    'CEACAM8': 'CEACAM8',
    'LCN2': 'LCN2',  # Lipocalin 2
}

print("\nCalculating neutrophil marker scores...")

# Calculate module scores for each marker
for gene, label in neutrophil_markers.items():
    if gene in adata.var_names:
        adata.obs[f'{label}_score'] = adata[:, gene].X.toarray().flatten()

# Manual cell type annotation based on top markers
print("\nPerforming clustering...")
sc.tl.leiden(adata, resolution=0.5)

# Create UMAP visualization colored by clusters
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color='leiden', ax=ax, show=False)
plt.savefig(base_dir / 'results' / 'figures' / 'GSE188288_clusters_UMAP.pdf')
print("✓ Cluster UMAP saved")

# Save annotated data
adata.write(processed_dir / "GSE188288_annotated.h5ad")
print(f"\n✓ Annotated data saved")

EOF

python3 03_cell_annotation.py
```

---

## ANALYSIS WORKFLOW

### Step 1: NETosis Signature Scoring

```bash
cat > analysis/01_netosis_scoring.py << 'EOF'
#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

sc.settings.verbosity = 3

# Paths
base_dir = Path.home() / "neutrophil_netosis_project"
processed_dir = base_dir / "data" / "processed"
results_dir = base_dir / "results"

# Load annotated data
print("Loading annotated data...")
adata = sc.read_h5ad(processed_dir / "GSE188288_annotated.h5ad")

# Define NETosis gene signature (61 genes from published study)
netosis_genes = [
    'ELANE', 'MPO', 'CTSG', 'PRTN3', 'DEFA1', 'DEFA3', 'CAMP', 'LTF', 'LCN2',
    'HIST1H1C', 'HIST1H2AC', 'HIST1H3A', 'HIST1H4A', 'HIST1H4B',
    'GPRC5B', 'FUT4', 'SOD2', 'STEAP3', 'STX11', 'SLC25A37',
    'ACTB', 'GNPTAB', 'RAB11FIP5', 'SLAMF9', 'CYBB', 'CHI3L1',
    # Add remaining 35 genes from published signature
]

# Filter genes present in dataset
netosis_genes_present = [g for g in netosis_genes if g in adata.var_names]
print(f"\nNETosis genes found: {len(netosis_genes_present)} / {len(netosis_genes)}")

# Calculate NETosis score using AddModuleScore equivalent
from sklearn.preprocessing import StandardScaler

print("\nCalculating NETosis signature scores...")

# Get expression matrix for NETosis genes
X_netosis = adata[:, netosis_genes_present].X

# Calculate mean expression per cell
netosis_score = np.mean(X_netosis, axis=1).A1 if hasattr(X_netosis, 'A1') else np.mean(X_netosis, axis=1)

# Add to obs
adata.obs['NETosis_Score'] = netosis_score

# Normalize score to 0-1 range
adata.obs['NETosis_Score_Norm'] = (netosis_score - netosis_score.min()) / (netosis_score.max() - netosis_score.min())

print(f"NETosis Score range: {netosis_score.min():.3f} - {netosis_score.max():.3f}")
print(f"Mean NETosis Score: {netosis_score.mean():.3f}")

# Save results
adata.write(processed_dir / "GSE188288_with_netosis_score.h5ad")
print("✓ Data with NETosis scores saved")

# Export scores
scores_df = adata.obs[['NETosis_Score', 'NETosis_Score_Norm']].copy()
scores_df.to_csv(results_dir / "tables" / "netosis_scores.csv")
print("✓ NETosis scores exported")

EOF

python3 analysis/01_netosis_scoring.py
```

### Step 2: Differential Expression Analysis

```bash
cat > analysis/02_differential_expression.py << 'EOF'
#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
from pathlib import Path

sc.settings.verbosity = 3

# Paths
base_dir = Path.home() / "neutrophil_netosis_project"
processed_dir = base_dir / "data" / "processed"
results_dir = base_dir / "results"

# Load data with NETosis scores
print("Loading data with NETosis scores...")
adata = sc.read_h5ad(processed_dir / "GSE188288_with_netosis_score.h5ad")

# Stratify cells by NETosis score
print("\nStratifying cells by NETosis activation...")
adata.obs['NETosis_Status'] = pd.cut(
    adata.obs['NETosis_Score'],
    bins=[0, np.percentile(adata.obs['NETosis_Score'], 33),
          np.percentile(adata.obs['NETosis_Score'], 67), np.inf],
    labels=['Low', 'Intermediate', 'High']
)

# Find markers for each NETosis group
print("\nFinding differentially expressed genes by NETosis status...")
sc.tl.rank_genes_groups(
    adata,
    groupby='NETosis_Status',
    method='wilcoxon',
    use_raw=False
)

# Extract results
degs = sc.get.rank_genes_groups_df(adata, group='High')
degs.to_csv(results_dir / "tables" / "DEG_high_netosis.csv", index=False)
print(f"✓ Differential expression results saved")

# Create heatmap of top DEGs
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(12, 8))
top_genes = degs.head(20)['names'].tolist()
sc.pl.heatmap(adata, var_names=top_genes, groupby='NETosis_Status', show=False)
plt.savefig(results_dir / "figures" / "DEG_heatmap.pdf", bbox_inches='tight', dpi=300)
print("✓ Heatmap saved")

EOF

python3 analysis/02_differential_expression.py
```

### Step 3: Visualization and Figure Generation

```bash
cat > analysis/03_visualization.py << 'EOF'
#!/usr/bin/env python3

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

sc.settings.figure_params(dpi=300, facecolor='white')

# Paths
base_dir = Path.home() / "neutrophil_netosis_project"
processed_dir = base_dir / "data" / "processed"
figures_dir = base_dir / "results" / "figures"

# Load data
print("Loading data for visualization...")
adata = sc.read_h5ad(processed_dir / "GSE188288_with_netosis_score.h5ad")

# Figure 1: UMAP colored by NETosis score
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color='NETosis_Score_Norm', cmap='RdYlBu_r', ax=ax, show=False)
plt.title('Neutrophils Colored by NETosis Score')
plt.savefig(figures_dir / 'UMAP_netosis_score.pdf', bbox_inches='tight', dpi=300)
print("✓ NETosis score UMAP saved")

# Figure 2: Feature plots for key genes
key_genes = ['ELANE', 'MPO', 'CTSG', 'PRTN3', 'S100A8', 'S100A9']
fig = sc.pl.umap(adata, color=key_genes, ncols=3, show=False)
plt.savefig(figures_dir / 'feature_plots_key_genes.pdf', bbox_inches='tight', dpi=300)
print("✓ Feature plots saved")

# Figure 3: Violin plots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for idx, gene in enumerate(key_genes):
    ax = axes[idx // 3, idx % 3]
    parts = ax.violinplot(adata[:, gene].X.toarray().flatten(), positions=[0])
    ax.set_title(gene)
    ax.set_ylabel('Expression')

plt.tight_layout()
plt.savefig(figures_dir / 'violin_plots_genes.pdf', bbox_inches='tight', dpi=300)
print("✓ Violin plots saved")

print("\n✓ All visualizations complete!")

EOF

python3 analysis/03_visualization.py
```

---

## TROUBLESHOOTING

### Issue: Command not found

```bash
# Add current directory to PATH
export PATH="$PATH:."

# Or use full path
python3 ./script.py
```

### Issue: Permission denied

```bash
# Make script executable
chmod +x script.py

# Or run with python explicitly
python3 script.py
```

### Issue: Memory issues

```bash
# Use chunked processing
# In Python scripts, use:
# adata.chunked_X = True
# Or work with subsets:
adata_subset = adata[adata.obs['batch'] == 'batch1']
```

### Issue: Conda environment not activating

```bash
# Reinitialize conda
conda init zsh

# Close and reopen terminal
exit

# Then
conda activate scrna-seq
```

### Issue: Package not found

```bash
# Update conda
conda update -n base -c defaults conda

# Or install package specifically
pip install package_name

# For bioconda packages:
conda install -c bioconda package_name
```

---

## Complete Quick Start Script

Save this as `~/start_project.sh`:

```bash
#!/bin/bash

echo "🚀 Starting Neutrophil/NETosis Project Setup..."

# Navigate to project
cd ~/neutrophil_netosis_project

# Activate environment
conda activate scrna-seq

# Create directory structure
mkdir -p data/{raw,processed,references}
mkdir -p results/{figures,tables,metrics}
mkdir -p code/{preprocessing,analysis,visualization}
mkdir -p logs

echo "✓ Directory structure created"

# Download data
echo "Downloading datasets..."
cd data/raw

# Download GSE188288
mkdir -p GSE188288 && cd GSE188288
wget -q "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188288&format=file" -O GSE188288_series_matrix.txt.gz
gunzip *.gz
cd ..

# Download GSE131907
mkdir -p GSE131907 && cd GSE131907
wget -q "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131907&format=file" -O GSE131907_series_matrix.txt.gz
gunzip *.gz
cd ..

cd ../..

echo "✓ Data downloads initiated"

# Run preprocessing
echo "Running preprocessing..."
python3 code/preprocessing/01_quality_control.py
python3 code/preprocessing/02_preprocessing.py
python3 code/preprocessing/03_cell_annotation.py

echo "✓ Preprocessing complete"

# Run analysis
echo "Running analysis..."
python3 code/analysis/01_netosis_scoring.py
python3 code/analysis/02_differential_expression.py
python3 code/analysis/03_visualization.py

echo "✓ Analysis complete"
echo "📊 Results saved to results/ directory"
```

Run with:
```bash
bash ~/start_project.sh
```

---

## Expected Output Directory Structure After Running

```
~/neutrophil_netosis_project/
├── data/
│   ├── raw/
│   │   ├── GSE188288/
│   │   ├── GSE131907/
│   │   └── TCGA-LUAD/
│   ├── processed/
│   │   ├── GSE188288_processed.h5ad
│   │   ├── GSE188288_annotated.h5ad
│   │   └── GSE188288_with_netosis_score.h5ad
│   └── references/
│       ├── GRCh38.fa
│       └── gencode.v43.gtf
├── results/
│   ├── figures/
│   │   ├── UMAP_netosis_score.pdf
│   │   ├── feature_plots_key_genes.pdf
│   │   └── DEG_heatmap.pdf
│   ├── tables/
│   │   ├── netosis_scores.csv
│   │   └── DEG_high_netosis.csv
│   └── metrics/
│       └── GSE188288_summary_stats.csv
├── code/
│   ├── preprocessing/
│   ├── analysis/
│   └── download_geo_data.py
└── logs/
```

---

## Next Steps

1. Run `bash ~/start_project.sh` to begin full pipeline
2. Monitor `results/` directory for outputs
3. Check `logs/` for error messages if issues arise
4. Proceed to cell-cell interaction analysis with CellChat (optional)
5. Prepare manuscript figures

---

## References

- Scanpy documentation: https://scanpy.readthedocs.io/
- Seurat documentation: https://satijalab.org/seurat/
- GEO database: https://www.ncbi.nlm.nih.gov/geo/
- GDC Portal: https://portal.gdc.cancer.gov/
