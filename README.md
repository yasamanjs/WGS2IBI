# WGS2IBI: A Cloud-Based Modular Workflow for Streamlined Individualized Bayesian WGS Data Analysis

## Overview
This repository provides the CWL code, Docker containers, and scripts for **WGS2IBI** — a **cloud-based, modular workflow** that streamlines **Individualized Bayesian Inference (IBI)** for **whole genome sequencing (WGS)** data analysis.  

Optimized for scalability, reproducibility, and cost-efficiency, WGS2IBI enables both population-level and individual-level genomic discovery.  
The workflow has been tested on large datasets, such as the **TOPMed Freeze 8 and Freeze 9 cohorts**.

[Read the full study (link pending publication)]  

![fig1](https://github.com/user-attachments/assets/26c19d09-ae10-44e0-a802-77193d927d67)

---

## Key Features  
- **Scalability**: Supports large cohorts (thousands of participants, >100M variants).  
- **Cost Efficiency**: Preprocessing and IBI analysis can be completed for under **$25 per cohort** on AWS.  
- **Cloud Readiness**: Deployable on cloud platforms (AWS, Seven Bridges) or local HPC clusters.  
- **Reproducibility**: Fully containerized with Docker and standardized using CWL.  
- **Memory Optimization**: Uses numpy arrays and dynamic batching to minimize RAM usage.  
- **Population and Individual Analysis**: Integrates GWAS (Fisher test, GDSearch) and IBI for comprehensive insights.  

---

## Workflow Components
- **Preprocessing Tool**:  
  `cwl/preprocessing.cwl` — Processes raw WGS VCF files using VCF2SNP, PLINK, and BCFtools.
- **Population Analysis Tools**:  
  `cwl/population_analysis.cwl` — Conducts Fisher's exact test and GDSearch for GWAS.
- **Top-k Pool Extraction Scripts**:  
  `scripts/initial_analysis.py`, `scripts/Confirm_topk_pool.py` — Selects top-ranked variants for downstream analysis.
- **Individual-Level Analysis Tool**:  
  `cwl/individual_analysis.cwl` — Performs IBI on top-k selected variants.

---

## Performance Summary
- Reduced preprocessing time and costs by dynamically optimizing cloud resources.
- Achieved preprocessing costs of **$20.83** and individual-level analysis costs of **under $5** for cohorts like JHS (2,777 samples).
- Demonstrated scalability across datasets such as JHS, WHI, MESA, FHS, and GENOA.

---

## Access Options

### 1. GitHub Repository (Current Page)
The full workflow source code, CWL files, and example scripts are provided here for researchers who wish to set up the workflow manually.  
See the **Getting Started** section below for setup instructions.

### 2. Public Project on Seven Bridges (Recommended for Testing)
A public project hosting **WGS2IBI** is available on the Seven Bridges platform.  
This provides the easiest way to run the workflow **without any installation**.  

**Advantages:**
- Fully set up and ready to use — no need to install CWL runners, Docker, or manage AWS instances.
- Test anonymously: Seven Bridges agreed to support reviewer access without needing personal disclosure.
- Free trial credits available:  
  Simply email **support@sevenbridges.com** with the subject "Request for Trial Access to WGS2IBI Project."

**Instructions:**
- Visit Seven Bridges platform.
- Copy the public WGS2IBI project into your workspace.
- Test with provided sample data or upload your own data.

> Detailed instructions for accessing and using the Seven Bridges project are also provided in the [docs/sevenbridges_usage.md](docs/sevenbridges_usage.md) file.

---

## Getting Started (Manual Setup)

### Prerequisites
- **Docker** (for running containerized environments)
- **CWL Runner** (e.g., `cwltool`)
- **AWS account** (for cloud deployment) — optional if running locally

### Quick Start
```bash
git clone https://github.com/your-repo-name
cd your-repo-name

---

## Applications
- **Precision Medicine**: Enables personalized insights by analyzing the most statistically significant genetic variants.
- **Population Studies**: Supports large-scale analyses across diverse cohorts, providing race-specific and phenotype-specific insights.
- **Genetic Research**: Facilitates efficient exploration of deep genetic variants using top-k variant selection.

---

## Citation
If you use this workflow in your research, please cite:  
[Pending Publicatio]  
by Yasaman J. Soofi, Jin Ren, Md. Asad Rahman, David Roberson, and Jinling Liu.

---

## License
This project is licensed under the MIT License.
You are free to use, modify, and distribute this software under the terms of the license.
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

🔒 **Note:** This software is intended for academic and non-commercial research purposes. Commercial use requires prior written permission from the authors.

