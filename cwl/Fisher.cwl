{
    "class": "CommandLineTool",
    "cwlVersion": "v1.2",
    "baseCommand": [
        "python",
        "myscript.py"
    ],
    "inputs": [
        {
            "id": "Phenotype_name",
            "type": "string",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
            "doc": "The name of the phenotype to be analyzed, guiding the GWAS to focus on the specified trait."
        },
        {
            "id": "txt_file",
            "type": "File[]?",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
            "doc": "Text files listing the subject IDs (`subIDs`) and variant IDs (`varIDs`) that correspond to the SNP data matrix, ensuring accurate data alignment.",
        },
        {
            "id": "npy_file",
            "type": "File[]?",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
            "doc": "This input contains the preprocessed SNP data in a binary format that allows efficient storage and rapid access during analysis.",
        },
        {
            "id": "Pheno",
            "type": "File?",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
            "doc": "A CSV file containing phenotypic information corresponding to the subjects in the SNP matrix.",
        }
    ],
    "outputs": [
        {
            "id": "GWAS_stat",
            "doc": "This output file contains the statistical results of the GWAS, including the p-values calculated for each variant using Fisher's Exact Test.",
            "type": "File?",
            "outputBinding": {
                "glob": "*_fisher_exact_result.csv"
            }
        }
    ],
    "doc": "### Fisher Tool Description\n\nThe Fisher tool is designed to perform genome-wide association studies (GWAS) using Fisher's Exact Test, providing a statistical framework for identifying significant genetic associations. This tool takes as input a compressed SNP matrix in `.npy` format, alongside subject and variant IDs, phenotype data, and the specified phenotype name, and outputs a CSV file with the GWAS statistics for each variant.\n\n&nbsp;\n\n#### Inputs:\n- **Compressed SNP matrix (`.npy`):** This input contains the preprocessed SNP data in a binary format that allows efficient storage and rapid access during analysis.\n- **Subject IDs and Variant IDs (`.txt`):** Text files listing the subject IDs (`subIDs`) and variant IDs (`varIDs`) that correspond to the SNP data matrix, ensuring accurate data alignment.\n- **Phenotype Data (`.csv`):** A CSV file containing phenotypic information corresponding to the subjects in the SNP matrix.\n- **Phenotype Name:** The name of the phenotype to be analyzed, guiding the GWAS to focus on the specified trait.\n\n&nbsp;\n\n#### Outputs:\n- **GWAS stat.csv:** This output file contains the statistical results of the GWAS, including the p-values calculated for each variant using Fisher's Exact Test.\n\n&nbsp;\n\n#### Function:\nThe Fisher tool employs Fisher's Exact Test to assess the association between genetic variants and the specified phenotype across the entire population. This statistical method is particularly useful for analyzing small sample sizes or sparse data, providing a robust and precise evaluation of genetic associations. The resulting `GWAS stat.csv` file lists the p-values for each variant, enabling researchers to identify and prioritize significant genetic markers for further study.\n\n&nbsp;\n\nThis tool is a vital component of the genomic analysis pipeline, complementing the Bayesian methods used in IBI_GDsearch by providing a different statistical perspective on the genetic data. The integration of Fisher's Exact Test within the workflow ensures a comprehensive and multi-faceted approach to identifying significant genetic associations.",
    "label": "Fisher",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "ResourceRequirement",
            "ramMin": "${\r\n  var max_size = 0;\r\n  if (inputs.npy_file) {\r\n    inputs.npy_file.forEach(function(file) {\r\n      if(file.size > max_size) {\r\n        max_size = file.size;\r\n      }\r\n    });\r\n  }\r\n  return Math.max(max_size * 10 * 1e-6, 12000);\r\n}\r\n"
        },
        {
            "class": "DockerRequirement",
            "dockerPull": "jupyter/datascience-notebook:295612d3ade4"
        },
        {
            "class": "InitialWorkDirRequirement",
            "listing": [
                {
                    "entryname": "myscript.py",
                    "entry": "import os\nimport time\nimport numpy as np\nimport pandas as pd\nfrom scipy.stats import fisher_exact\n\nimport tooloptions\n\n\nclass ReadEachVariantsData:\n    \"\"\"\n    Read genotype data\n    \"\"\"\n    def __init__(self,genotype_file_root_path,filename):\n        self.genotype_file_root_path=genotype_file_root_path\n        self.filename=filename\n\n    def read_var_ids(self):\n        var_ids_file = open(os.path.join(self.genotype_file_root_path, \"{}_varIDs.txt\".format(self.filename)), \"r\")\n        var_ids = var_ids_file.read().split(\"\\n\")\n        return var_ids\n\n    def read_sub_ids(self):\n        sub_ids_file = open(os.path.join(self.genotype_file_root_path, \"{}_subIDs.txt\".format(self.filename)), \"r\")\n        sub_ids = sub_ids_file.read().split(\"\\n\")\n        return sub_ids\n\n    def read_genotype_data(self):\n        variants=np.load(os.path.join(self.genotype_file_root_path,\"{}_variants_value.npy\").format(self.filename))\n        return variants\n\n\n\n################# main ############################################\nPhenotype_name = tooloptions.Phenotype_name\nphenotype_file_path = tooloptions.gzpath_pheno\n\n\ndf_phenotype=pd.read_csv(phenotype_file_path)\ntrait=df_phenotype[Phenotype_name].values.reshape(-1,1)\n\n\n\n## Get info for Phenotype and SNP data and etc\ntxt_first_file_path = tooloptions.txt_first_file_path\nnpy_first_file_path = tooloptions.npy_first_file_path\n\n\ngenotype_file_root_path = os.path.dirname(txt_first_file_path) \n\nnpy_files_basename_list = tooloptions.npy_files_basename_list\n# Extract the file list without the '_variants_value.npy' part\nfile_list = [filename.split('_variants_value.npy')[0] for filename in npy_files_basename_list]\nprint(\"File list:\",file_list)\n\n\n\n########### Fisher Exact test for all chromosomes\n\n\nfor filename in file_list:\n    print(filename)\n    result = []\n    start_time = time.time()\n    read_each_variants_data=ReadEachVariantsData(genotype_file_root_path,filename)\n    variants=read_each_variants_data.read_genotype_data()\n    sub_ids=read_each_variants_data.read_sub_ids()\n    var_ids=read_each_variants_data.read_var_ids()\n\n    variant_sum=variants.sum(axis=1).reshape(-1,1)\n    trait_sum=trait.sum(axis=0).reshape(1,-1)\n\n    V1D1=variants@trait\n    V1D0=np.ones_like(V1D1)*variant_sum-V1D1\n    V0D1=np.ones_like(V1D1)*trait_sum-V1D1\n    V0D0=np.ones_like(V1D0)*variants.shape[1]-np.ones_like(V1D0)*trait_sum-V1D0\n    inf_or=[]\n    for variant_index,variant_name in enumerate(var_ids):\n        # print(variant_index,variant_name)\n        v0d0=V0D0[variant_index][0]\n        v0d1=V0D1[variant_index][0]\n        v1d0=V1D0[variant_index][0]\n        v1d1=V1D1[variant_index][0]\n\n        contigency_table = [[v1d1, v0d1], [v1d0, v0d0]]\n        # print(v1d1, v0d1, v1d0, v0d0)\n\n        # perform fisher's exact test\n        odds_ratio, p_value = fisher_exact(contigency_table)\n        if odds_ratio==np.inf:\n            inf_or.append([variant_name,v1d1,v0d1,v1d0,v0d0,odds_ratio])\n\n        # save results\n        result.append([filename,variant_index,variant_name,odds_ratio,p_value])\n\n        # print(contigency_table)\n        # print(\"odds_ratio:{},p_value:{}\".format(odds_ratio,p_value))\n    df_inf_or=pd.DataFrame(inf_or,columns=[\"variant_name\",\"v1d1\",\"v0d1\",\"v1d0\",\"v0d0\",\"odds_ratio\"])\n    elapsed_time = round(time.time() - start_time)\n    df_result =pd.DataFrame(result,columns=[\"Chromosome\",\"Variant Index\",\"Variants ID\",\"Odds Ratio\",\"P Value\"])\n    df_result.to_csv(os.path.join(\"{}_fisher_exact_result.csv\".format(filename)),index=False)\n    print(\"{},elapsed_time:{}\".format(filename,elapsed_time))",
                    "writable": false
                },
                {
                    "entryname": "tooloptions.py",
                    "entry": "txt_first_file_path = \"$(inputs.txt_file[0].path)\"\nnpy_first_file_path = \"$(inputs.npy_file[0].path)\"\n\ngzpath_pheno = \"$(inputs.Pheno.path)\"\nPhenotype_name = \"$(inputs.Phenotype_name)\"\n\nnpy_files_basename_list = ${\n          var basename = [];\n          if(inputs.npy_file){\n            inputs.npy_file.forEach(function(file){\n              basename.push('\"' + file.basename + '\"');\n            });\n          }\n          return '[' + basename.join(\",\") + ']';\n        }\n        \n        \n        ",
                    "writable": false
                }
            ]
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "stdout": "standard.out",
}