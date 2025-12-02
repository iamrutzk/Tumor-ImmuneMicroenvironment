# QUICK START: Mac Terminal Commands for Neutrophil/NETosis Project

## 1. INITIAL SETUP (Run Once)

```bash
# Create project directory
mkdir -p ~/neutrophil_netosis_project
cd ~/neutrophil_netosis_project

# Create subdirectories
mkdir -p data/{raw,processed,references}
mkdir -p results/{figures,tables,metrics}
mkdir -p code/{preprocessing,analysis,visualization}
mkdir -p logs

# Check system
sw_vers
df -h
xcode-select --install
brew --version

# If Homebrew not installed:
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

## 2. INSTALL CONDA (if not already installed)

```bash
# For Apple Silicon Mac (M1/M2/M3):
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh

# For Intel Mac:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh

# Initialize
conda init zsh

# Close and reopen terminal, then verify:
conda --version
```

## 3. CREATE CONDA ENVIRONMENTS

```bash
# Create Python environment for scRNA-seq
conda create -n scrna-seq python=3.11 -y

# Create R environment for Seurat
conda create -n seurat-r r-base=4.3 r-essentials -y

# Activate Python environment
conda activate scrna-seq
```

## 4. INSTALL REQUIRED PACKAGES

```bash
# Activate environment
conda activate scrna-seq

# Install bioinformatics tools
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
  notebook

# Install additional tools via pip
pip install GEOparse
pip install SRA-Toolkit
pip install boto3

# Verify installation
python -c "import scanpy; print(scanpy.__version__)"
```

## 5. DOWNLOAD PRIORITY 1: GSE188288 (Circulating Neutrophils)

```bash
cd ~/neutrophil_netosis_project/data/raw

# Create directory
mkdir -p GSE188288
cd GSE188288

# Download series matrix
echo "Downloading GSE188288..."
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE188288&format=file" -O GSE188288.tar.gz

# Extract
tar -xzf GSE188288.tar.gz

# List files
ls -lah

echo "✓ GSE188288 download complete"
cd ~/neutrophil_netosis_project
```

## 6. DOWNLOAD PRIORITY 2: GSE131907 (LUAD Tumor Immune)

```bash
cd ~/neutrophil_netosis_project/data/raw

# Create directory
mkdir -p GSE131907
cd GSE131907

# Download
echo "Downloading GSE131907..."
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131907&format=file" -O GSE131907.tar.gz

# Extract
tar -xzf GSE131907.tar.gz

# List files
ls -lah

echo "✓ GSE131907 download complete"
cd ~/neutrophil_netosis_project
```

## 7. DOWNLOAD TCGA-LUAD (Bulk Validation Data)

```bash
cd ~/neutrophil_netosis_project/data/raw
mkdir -p TCGA-LUAD
cd TCGA-LUAD

# Create query file for RNA-seq files
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
      }
    ]
  },
  "format": "JSON",
  "size": 100
}
EOF

# Query GDC API
curl --request POST --header "Content-Type: application/json" \
  --data @query.json \
  "https://api.gdc.cancer.gov/files" > gdc_files.json

# Extract first 10 file IDs and download
python3 << 'PYTHON'
import json
with open('gdc_files.json') as f:
    data = json.load(f)
    file_ids = [hit['id'] for hit in data['data']['hits'][:10]]
    for fid in file_ids:
        print(fid)
PYTHON

echo "✓ TCGA-LUAD file list ready"
cd ~/neutrophil_netosis_project
```

## 8. QUICK DATA CHECK

```bash
cd ~/neutrophil_netosis_project

# Check what's been downloaded
du -sh data/raw/GSE188288
du -sh data/raw/GSE131907
du -sh data/raw/TCGA-LUAD

# Total size
du -sh data/raw
```

## 9. CREATE FIRST ANALYSIS SCRIPT

```bash
cd ~/neutrophil_netosis_project/code/preprocessing

cat > 01_load_and_check.py << 'EOF'
#!/usr/bin/env python3

import os
from pathlib import Path
import gzip
import tarfile

base_dir = Path.home() / "neutrophil_netosis_project"
data_dir = base_dir / "data" / "raw"

print("="*60)
print("DATA AVAILABILITY CHECK")
print("="*60)

# Check GSE188288
print("\nGSE188288 (Circulating Neutrophils):")
gse188288_path = data_dir / "GSE188288"
if gse188288_path.exists():
    files = list(gse188288_path.glob("*"))
    print(f"  Files: {len(files)}")
    for f in files[:5]:
        print(f"    - {f.name}")
else:
    print("  NOT FOUND - Please download first")

# Check GSE131907
print("\nGSE131907 (LUAD Tumor Immune):")
gse131907_path = data_dir / "GSE131907"
if gse131907_path.exists():
    files = list(gse131907_path.glob("*"))
    print(f"  Files: {len(files)}")
    for f in files[:5]:
        print(f"    - {f.name}")
else:
    print("  NOT FOUND - Please download first")

# Check TCGA-LUAD
print("\nTCGA-LUAD (Bulk RNA-seq):")
tcga_path = data_dir / "TCGA-LUAD"
if tcga_path.exists():
    files = list(tcga_path.glob("*.json"))
    print(f"  Query files: {len(files)}")
else:
    print("  NOT FOUND - Please download first")

print("\n" + "="*60)
print("Next step: Run 02_process_data.py")
print("="*60)

EOF

python3 01_load_and_check.py
```

## 10. PROCESS DOWNLOADED DATA

```bash
cd ~/neutrophil_netosis_project/code/preprocessing

cat > 02_process_data.py << 'EOF'
#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

sc.settings.verbosity = 3

base_dir = Path.home() / "neutrophil_netosis_project"
raw_dir = base_dir / "data" / "raw"
processed_dir = base_dir / "data" / "processed"
processed_dir.mkdir(exist_ok=True)

print("\n" + "="*60)
print("PROCESSING GSE188288")
print("="*60 + "\n")

# Try to read from various formats
gse_path = raw_dir / "GSE188288"

# Look for h5ad file
h5ad_files = list(gse_path.glob("*.h5ad"))
if h5ad_files:
    print(f"Loading H5AD file: {h5ad_files[0]}")
    adata = sc.read_h5ad(h5ad_files[0])
else:
    # Look for MTX format
    mtx_files = list(gse_path.glob("**/*.mtx*"))
    if mtx_files:
        print(f"Found matrix files, attempting to read...")
        # Use scanpy to read from directory
        adata = sc.read_10x_mtx(gse_path)
    else:
        # Try reading text files
        txt_files = list(gse_path.glob("*.txt*"))
        if txt_files:
            print(f"Found text files, reading first one...")
            import gzip
            file_path = txt_files[0]
            if file_path.suffix == '.gz':
                with gzip.open(file_path, 'rt') as f:
                    data = pd.read_csv(f, sep='\t', index_col=0)
            else:
                data = pd.read_csv(file_path, sep='\t', index_col=0)
            adata = sc.AnnData(X=data.T)
        else:
            print("No compatible data files found!")
            print(f"Files in directory: {list(gse_path.glob('*'))[:10]}")
            exit(1)

print(f"\nData loaded:")
print(f"  Shape: {adata.shape}")
print(f"  Cells: {adata.n_obs}")
print(f"  Genes: {adata.n_vars}")

# Basic QC
if 'n_genes' not in adata.obs.columns:
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
if 'n_counts' not in adata.obs.columns:
    adata.obs['n_counts'] = adata.X.sum(axis=1)

# Calculate MT percentage if gene names available
if hasattr(adata.var_names, 'str'):
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
    if adata.var['mt'].any():
        adata.obs['pct_mt'] = (adata[:, adata.var['mt']].X.sum(axis=1) / adata.obs['n_counts'] * 100)

print(f"\nQC Metrics:")
print(f"  Mean genes per cell: {adata.obs['n_genes'].mean():.0f}")
print(f"  Mean counts per cell: {adata.obs['n_counts'].mean():.0f}")

# Save processed
output_file = processed_dir / "GSE188288_loaded.h5ad"
adata.write(output_file)
print(f"\n✓ Saved to {output_file}")

print("\n" + "="*60)
print("SUCCESS: Data loaded and saved")
print("="*60)

EOF

python3 02_process_data.py
```

## 11. VERIFY EVERYTHING WORKS

```bash
# Check if processed files exist
ls -lah ~/neutrophil_netosis_project/data/processed/

# Check results
ls -lah ~/neutrophil_netosis_project/results/

# View log files if any
ls -lah ~/neutrophil_netosis_project/logs/
```

## 12. JUPYTER NOTEBOOK (Optional - for interactive analysis)

```bash
# Activate environment
conda activate scrna-seq

# Start Jupyter
cd ~/neutrophil_netosis_project
jupyter notebook

# This opens http://localhost:8888 in your browser
# Create new notebook and start analyzing data
```

## 13. RUN FULL PIPELINE (After initial setup)

```bash
# Save this as ~/run_pipeline.sh
cat > ~/run_pipeline.sh << 'BASH'
#!/bin/bash
cd ~/neutrophil_netosis_project
conda activate scrna-seq

echo "Running neutrophil analysis pipeline..."
python3 code/preprocessing/01_load_and_check.py
python3 code/preprocessing/02_process_data.py

echo "✓ Pipeline complete!"
BASH

# Make executable
chmod +x ~/run_pipeline.sh

# Run anytime with:
bash ~/run_pipeline.sh
```

## COMMON COMMANDS FOR DAILY USE

```bash
# Activate your environment
conda activate scrna-seq

# Navigate to project
cd ~/neutrophil_netosis_project

# Check conda environments
conda env list

# Deactivate environment
conda deactivate

# Update packages
conda update -n scrna-seq --all

# See what's in your environment
conda list

# Remove environment (if needed)
conda remove -n scrna-seq --all

# Use Python interactively
python3

# Run a script
python3 code/preprocessing/01_load_and_check.py

# Check file sizes
du -sh data/raw/*
du -sh data/processed/*

# Monitor while downloading
watch -n 5 'du -sh data/raw/*'
```

## TROUBLESHOOTING QUICK FIXES

```bash
# Package not found?
pip install --upgrade package_name

# Conda not working?
conda init zsh
# Close and reopen terminal

# Need more space?
du -sh ~  # Check total home dir size

# Jupyter not opening?
jupyter notebook --ip=localhost --port=8888

# Script permission error?
chmod +x script.py

# Python import error?
python3 -m pip install package_name

# Clear conda cache
conda clean --all
```

## EXPECTED TIMELINE

- **Setup (first time):** 30-45 minutes
- **Conda installation:** 5 minutes
- **Package installation:** 10-15 minutes
- **GSE188288 download:** 30-60 minutes (depends on connection)
- **GSE131907 download:** 30-60 minutes
- **Data processing:** 15-30 minutes
- **First analysis:** 10 minutes

**Total first run:** ~2-3 hours

## WHAT YOU'LL HAVE AFTER RUNNING

```
~/neutrophil_netosis_project/
├── data/processed/
│   └── GSE188288_loaded.h5ad (ready for analysis)
├── results/
│   └── (will populate with figures/tables)
└── code/
    └── (your analysis scripts)
```

## NEXT ANALYSIS STEPS

After downloading and processing, run:

1. **Clustering:** `python3 code/preprocessing/03_clustering.py`
2. **Cell typing:** `python3 code/preprocessing/04_annotation.py`
3. **NETosis scoring:** `python3 code/analysis/01_netosis_scoring.py`
4. **Differential expression:** `python3 code/analysis/02_deg_analysis.py`
5. **Visualization:** `python3 code/analysis/03_make_figures.py`

See full guide for complete scripts: mac_pipeline_guide.md
