# Cloud-Based Modular Workflow for Individualized Bayesian Genomic Analysis  
**An Efficient Approach for Preprocessing and Association Studies of Whole Genome Sequences**  

## Overview
This repository contains the CWL code and resources for a cloud-based, modular workflow designed to optimize the **Individualized Bayesian Inference (IBI)** approach for **genome-wide association studies (GWAS)** using **Whole Genome Sequencing (WGS)** data. The workflow emphasizes **efficiency**, **scalability**, and **reproducibility**, making IBI more accessible to genomic researchers.  

Designed for datasets such as the **TOPMed Freeze 8 and Freeze 9 cohorts**, this pipeline handles large-scale genomic data, preprocessing variants from **102M to 18M** while ensuring computational and cost efficiency.  

[Read the full study (link pending publication)]  

---

## Features  
- **Scalability**: Handles large datasets with millions of variants across thousands of participants.  
- **Cost Efficiency**: Dynamically optimized AWS instance selection reduces preprocessing costs to **$20.83 per cohort** and individual-level analysis costs to **under $5**.  
- **Memory Optimization**: Utilizes numpy arrays for efficient storage and faster computation, reducing RAM usage significantly.  
- **Reproducibility**: Implements CWL and Docker for consistent workflows across diverse environments.  
- **Compatibility**: Compatible with public datasets, including **MESA**, **JHS**, **FHS**, **GENOA**, and **WHI**.  

---

## Code Files  
### Workflow Components
1. **`cwl/preprocessing.cwl`**: Automates data preprocessing using VCF2SNP, PLINK, and BCFtools.  
2. **`cwl/individual_analysis.cwl`**: Conducts IBI-based individual-level analyses with dynamic batching.  
3. **`cwl/population_analysis.cwl`**: Executes GWAS-level analysis for population-wide insights.

### Scripts
- **`scripts/Confirm_topk_pool.py`**: Manages top-k variant selection for refined genetic analysis.  
- **`scripts/IBI_pool.py`**: Performs pooled IBI analysis across all chromosomes.  

---

## Performance Highlights  
- Preprocessed **102M variants to 18M** with optimized workflows, reducing computational burden.  
- Achieved **4x faster processing times** by leveraging parallel computing on HPC clusters and AWS.  
- Fine-tuned workflows for cost efficiency, achieving preprocessing costs of **$20.83 per cohort** and individual-level analysis costs of **under $5**.  
- Demonstrated scalability with large datasets from cohorts like JHS (2,777 samples) and WHI (11,054 samples).  

---

## Technologies Used  
- **Workflow Language**: CWL  
- **Containerization**: Docker  
- **Cloud Computing**: AWS, Seven Bridges  
- **Programming**: Python (pandas, numpy, custom scripts)  
- **Data Processing Tools**: PLINK, BCFtools, VCF2SNP  
- **Optimization Techniques**: GPU acceleration (CUDA), parallel processing  

---

## Getting Started  
### Prerequisites  
1. **Docker**: Install Docker to run containerized workflows.  
2. **CWL Runner**: Use a compatible runner like `cwltool`.  
3. **AWS Account**: Set up AWS for scalable and cost-efficient analysis.  

### Steps to Run the Workflow  
1. Clone the repository:  
   ```bash
   git clone https://github.com/your-repo-name
   cd your-repo-name
2. Review setup instructions in the `/docs` directory.
3. Run preprocessing workflows in the `/cwl` directory with sample data from `/data`.

---

## Applications
- **Precision Medicine**: Enables personalized insights by analyzing the most statistically significant genetic variants.
- **Population Studies**: Supports large-scale analyses across diverse cohorts, providing race-specific and phenotype-specific insights.
- **Genetic Research**: Facilitates efficient exploration of deep genetic variants using top-k variant selection.

---

## Citation
If you use this workflow in your research, please cite:  
**"Cloud-Based Modular Workflow for Individualized Bayesian Genomic Analysis"**  
by **Yasaman J. Soofi, Jin Ren, Md. Asad Rahman, David Roberson, and Jinling Liu.**

---

## License
This project is licensed under the [MIT License](LICENSE).
