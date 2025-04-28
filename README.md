# WGS2IBI: A Cloud-Based Modular Workflow for Streamlined Individualized Bayesian WGS Data Analysis

## Overview
This repository provides the CWL code, Docker containers, and scripts for **WGS2IBI** — a **cloud-based, modular workflow** that streamlines **Individualized Bayesian Inference (IBI)** for **whole genome sequencing (WGS)** data analysis.  

Optimized for scalability, reproducibility, and cost-efficiency, WGS2IBI enables both population-level and individual-level genomic discovery.  
The workflow has been tested on large datasets, such as the **TOPMed Freeze 8 and Freeze 9 cohorts**.

[Read the full study (link pending publication)]  




![Workflow Overview](https://github.com/user-attachments/assets/26c19d09-ae10-44e0-a802-77193d927d67)
> **Figure:** Visual representation of the WGS2IBI workflow, illustrating the modular tools and the flow of inputs and outputs.  
> *(Image reproduced from the manuscript pending publication.)*
>
> 
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
  `cwl/preprocessing.cwl` — VCF2SNP tool Processes raw WGS VCF files using BCFtools, PLINK, and Numpy.
- **Population-Level Analysis Tools**:  
  `cwl/population_analysis_fisher.cwl` — Conducts Fisher's exact test for GWAS.  
  `cwl/population_analysis_gdsearch.cwl` — Conducts GDSearch for GWAS.
- **Top-k Variant Pool Extraction Scripts**:  
  `scripts/initial_analysis.py` — Performs initial ranking analysis and generates the top-k variants pool.  
  `scripts/Confirm_topk_pool.py` — Confirms and manages selection of top-k variants for downstream analysis.
- **Individual-Level Analysis Tool (IBI)**:  
  `cwl/individual_analysis_ibi.cwl` — Performs Individualized Bayesian Inference (IBI) on the selected top-k variant pool.

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

**Access the public project here:**  
[https://sb.biodatacatalyst.nhlbi.nih.gov/projects/f8e1d6a6-6f2e-4fa7-aaa3-d63d663ef9e5/WGS2IBI](https://sb.biodatacatalyst.nhlbi.nih.gov/projects/f8e1d6a6-6f2e-4fa7-aaa3-d63d663ef9e5/WGS2IBI)  
*(Example URL — replace with your actual project link if different.)*

**Create a free Seven Bridges account here:**  
[https://accounts.sb.biodatacatalyst.nhlbi.nih.gov/auth/login](https://accounts.sb.biodatacatalyst.nhlbi.nih.gov/auth/login)


**Advantages:**
- Fully set up and ready to use — no need to install CWL runners, Docker, or manage AWS instances.
- After creating a Seven Bridges account, users can **anonymously copy and test** the public WGS2IBI project without needing to contact the authors.
- Users can **review previously submitted jobs**, including job inputs, parameters, outputs, and runtime settings — making it easy to learn how the workflow operates.
- A **sample dataset** from the **1000 Genomes Project benchmark data** is included in the public project for immediate testing.
- Free trial credits available: Simply email **support@sevenbridges.com** with the subject "Request for Trial Access to WGS2IBI Project."


---

## Sample Data

To facilitate testing and exploration of the WGS2IBI workflow, a **sample dataset** derived from the **1000 Genomes Project benchmark data** is provided.

- The sample dataset includes a subset of variants.
- It is designed for **benchmarking**, **demonstration**, and **exploration** of the tools' full capabilities without needing controlled-access data.
- Available in both:
  - The public project on Seven Bridges (pre-loaded and ready to use)
  - The `/data` directory in this GitHub repository

Researchers can **run the entire WGS2IBI workflow end-to-end** using this sample dataset to understand the preprocessing, population analysis, top-k variant extraction, and individualized Bayesian inference steps.

> **Compatibility Note:**  
> The WGS2IBI workflow is **fully compatible** with large-scale datasets from **TOPMed Freeze 8 and later cohorts**, including **FHS**, **JHS**, **MESA**, **GENOA**, and **WHI**.  
> These cohorts can be accessed directly the Seven Bridges BioData Catalyst platform, eliminating the need to manually download and upload massive genomic files.

> **Access Requirements:**  
> Please note that **TOPMed datasets are controlled access** and require appropriate authorization (e.g., dbGaP approval).  
> Therefore, to enable immediate testing without access restrictions, the 1000 Genomes sample data is provided.


## Getting Started (Manual Setup)

### Prerequisites
- **Docker** (for running containerized environments)
- **CWL Runner** (e.g., `cwltool`)
- **AWS account** (for cloud deployment) — optional if running locally

### Quick Start
```bash
git clone https://github.com/your-repo-name
cd your-repo-name
```

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

