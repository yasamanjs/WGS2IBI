{
    "class": "CommandLineTool",
    "cwlVersion": "v1.2",
    "baseCommand": [
        "python",
        "myscript.py"
    ],
    "inputs": [
        {
            "id": "txt_file",
            "type": "File[]",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            }
        },
        {
            "id": "top_k",
            "type": "int",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
            "default": 10
        },
        {
            "id": "npy_file",
            "type": "File[]",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            }
        },
        {
            "id": "stat_file",
            "type": "File[]",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
        }
    ],
    "outputs": [
        {
            "id": "pool_varIDs",
            "type": "File?",
            "outputBinding": {
                "glob": "*topK_pool_varIDs.csv"
            }
        },
        {
            "id": "pool_variants",
            "type": "File?",
            "outputBinding": {
                "glob": "*topK_pool_variants.csv"
            }
        },
        {
            "id": "pool_topGD_index",
            "type": "File?",
            "outputBinding": {
                "glob": "*topK_pool_topGD_index.csv"
            }
        }
    ],
    "doc": "### Confirm_top_k_pool Tool \n\nThe Confirm_top_k_pool tool integrates results from various analyses to consolidate a pool of top genetic variants. This tool is essential for selecting the most significant genetic variants based on specified criteria, facilitating targeted downstream analysis.\n\n&nbsp;\n\n#### Inputs:\n- **SNP matrix (`.npy`):** A binary format file containing the preprocessed SNP data.\n- **Subject IDs (`subIDs.txt`):** A text file listing the subject IDs corresponding to the SNP data matrix.\n- **Variant IDs (`varIDs.txt`):** A text file listing the variant IDs corresponding to the SNP data matrix.\n- **Top_k (integer):** The specified number of top variants to select.\n- **Stat file (`.csv`):** A CSV file containing variant IDs in the first column, IBI global ranks in the second column, and other rankings (e.g., Fisher’s exact test) in subsequent columns. It is crucial that the columns contain rankings, not lgM or p-values, and that the IBI global rank is in the second column.\n\n&nbsp;\n\n#### Outputs:\n- **pool_varIDs.csv:** A CSV file listing the variant IDs in the top K pool.\n- **pool_variants.csv:** A CSV file providing a matrix of the pooled top variants.\n- **pool_topGD_index.csv:** A CSV file detailing the indices of the top genetic determinants in the pool.\n\n&nbsp;\n\n#### Function:\nThe Confirm_top_k_pool tool performs the following key functions:\n\n1. **Variant Selection:**\n   - The tool reads the SNP matrix, subject and variant IDs, and statistical rankings from the provided input files.\n   - It selects the top K variants based on the specified value of K and their rankings, prioritizing the IBI global rank.\n\n2. **Pooling Top Variants:**\n   - It consolidates the selected top K variants into a unified pool, ensuring that the most significant variants are prioritized for further analysis.\n\n3. **Output Generation:**\n   - The tool generates three output files:\n     - `pool_varIDs.csv` lists the variant IDs in the top K pool.\n     - `pool_variants.csv` provides a matrix of the pooled top variants.\n     - `pool_topGD_index.csv` details the indices of the top genetic determinants in the pool.\n\n&nbsp;\n\n#### Benefits:\n- **Efficiency:** By consolidating the most significant variants into a single pool, the Confirm_top_k_pool tool streamlines the selection process for targeted analysis.\n- **Integration:** The tool integrates seamlessly with various analysis results, ensuring a comprehensive approach to identifying top genetic variants.\n- **Flexibility:** The ability to specify the value of K allows researchers to customize the selection process based on their specific analysis needs.\n\n&nbsp;\n\nBy selecting and consolidating the top K variants, the Confirm_top_k_pool tool enhances the precision and focus of downstream genetic analyses. This tool is a critical component of the genomic analysis pipeline, ensuring that only the most relevant variants are prioritized for further study.",
    "label": "confirm_top_k_pool",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "LoadListingRequirement"
        },
        {
            "class": "ResourceRequirement",
            "ramMin": "${\r\n  var max_size = 0;\r\n  if (inputs.npy_file) {\r\n    inputs.npy_file.forEach(function(file) {\r\n      if(file.size > max_size) {\r\n        max_size = file.size;\r\n      }\r\n    });\r\n  }\r\n  return Math.max(max_size * 10 * 1e-6, 12000);\r\n}\r\n"
        },
        {
            "class": "DockerRequirement",
            "dockerPull": "images.sb.biodatacatalyst.nhlbi.nih.gov/jamalyasa/nvidia_cuda_python:cu11.8.0-Py3.8-torch2.1.2"
        },
        {
            "class": "InitialWorkDirRequirement",
            "listing": [
                {
                    "entryname": "myscript.py",
                    "entry": "import os\r\nimport time\r\nimport pandas as pd\r\nimport numpy as np\r\nfrom datetime import datetime\r\nimport scipy.stats as stats\r\nfrom tqdm import tqdm\r\nimport logging\r\nimport torch\r\nimport tooloptions\r\n\r\n\r\nlogging.getLogger().setLevel(logging.INFO)\r\nlogging.info(\"Libraries are loaded\")\r\n\r\n\r\ndef read_variantsF_efficient(variants_path_file, file_name):\r\n    \"\"\"\r\n    read the compressed variants data\r\n    :param variants_path_file:\r\n    :param variants_size: if variants_size=None,select all the variants\r\n    :return: subIDs, varIDs, variants_tensor, df  # list, list, array and dataframe\r\n    \"\"\"\r\n    varIDs_file = open(os.path.join(variants_path_file, \"{}_varIDs.txt\".format(file_name)), \"r\")\r\n    varIDs = varIDs_file.read().split(\"\\n\")\r\n    subIDs_file = open(os.path.join(variants_path_file, \"{}_subIDs.txt\".format(file_name)), \"r\")\r\n    subIDs = subIDs_file.read().split(\"\\n\")\r\n    subIDs = list(str(x) for x in subIDs)\r\n    variants = np.load(os.path.join(variants_path_file, \"{}_variants_value.npy\".format(file_name)))\r\n\r\n    A0 = np.ones(len(subIDs), dtype=np.int8)\r\n    variants = np.row_stack((A0, variants))\r\n    varIDs.insert(0, 'A0')\r\n    # when variants data is large, df will consume a large amount of memory,\r\n    # but df is no being used, set df=pd.Dataframe()\r\n    df = pd.DataFrame()\r\n    # df = pd.DataFrame(variants, index=varIDs, columns=subIDs, dtype=np.int8)\r\n    variants_tensor = torch.as_tensor(variants, dtype=torch.float32) #float32 for gpu and float64 for cpu\r\n\r\n    return subIDs, varIDs, variants_tensor, df\r\n\r\n\r\n\r\ndef get_top_variants_info(file_list, txt_path, stat_all, varIDs_all, device, k=1000,\r\n                              union_flag=True):\r\n    \"\"\"\r\n    select top k variants from all chromosomes, if have multi traits, use union top k variants from all traits\r\n    :param chromosome_file_folder_path:\r\n    :param glgm_all:\r\n    :param varIDs_all:\r\n    :param device:\r\n    :param k:\r\n    :return:\r\n    \"\"\"\r\n    # Loop through each trait and take the top k variants from all chromosomes\r\n    top_k_varIDs_all = []\r\n    topGD_varID_all=[]\r\n    for traits_num in range(stat_all.shape[1]):\r\n        stat_all_single = stat_all[:, traits_num]\r\n        _, sorted_indices = torch.sort(stat_all_single) ### Yasaman: rankings are sorted ascendingly\r\n        sorted_varIDs_all = [varIDs_all[i] for i in sorted_indices]\r\n        top_k_varIDs = sorted_varIDs_all[:k]\r\n        top_k_varIDs_all.append(top_k_varIDs)\r\n        topGD_varID = sorted_varIDs_all[0]\r\n        topGD_varID_all.append(topGD_varID)\r\n\r\n    if union_flag:\r\n        # take union top k variants from all traits\r\n        top_k_varIDs_union = list(set(x for sublist in top_k_varIDs_all for x in sublist))\r\n\r\n    top_k_varID_final = []\r\n    top_k_variants = []\r\n\r\n    # extract top k variants from all files\r\n    for file_name in tqdm(file_list):\r\n        subIDs, varIDs, variants_tensor, df_variants = read_variantsF_efficient(txt_path,\r\n                                                                                file_name)\r\n        indexs = [i for i, value in enumerate(varIDs) if value in top_k_varIDs_union]\r\n        selected_varIDs = [value for i, value in enumerate(varIDs) if value in top_k_varIDs_union]\r\n        top_k_varID_final.extend(selected_varIDs)\r\n        top_k_variants.extend(variants_tensor[indexs, :].tolist())\r\n\r\n    top_variants = np.array(top_k_variants)\r\n    if \"AO\" not in top_k_varID_final:\r\n        A0 = np.ones(top_variants.shape[1], dtype=np.int8)\r\n        top_variants = np.row_stack((A0, top_variants))\r\n        top_k_varID_final.insert(0, 'A0')\r\n    top_variants_tensor = torch.tensor(top_variants, dtype=torch.float32, device=device)\r\n\r\n    # get topGD_index for all traits\r\n    topGD_index = [top_k_varID_final.index(topGD_varID) for topGD_varID in topGD_varID_all]\r\n    return top_k_varID_final, top_variants_tensor, topGD_index\r\n\r\n\r\n\r\n\r\n\r\n#######################################################################################\r\n#######################################################################################\r\n\r\ndevice = torch.device('cuda' if torch.cuda.is_available() else 'cpu')\r\n\r\n\r\n################# main ############################################\r\n## Get info for Phenotype and SNP data and etc\r\ntxt_first_file_path = tooloptions.txt_first_file_path\r\nnpy_first_file_path = tooloptions.npy_first_file_path\r\n\r\ntxt_path = os.path.dirname(txt_first_file_path)\r\n\r\n\r\nnpy_files_basename_list = tooloptions.npy_files_basename_list\r\n# Extract the file list without the '_variants_value.npy' part\r\nfile_list = [filename.split('_variants_value.npy')[0] for filename in npy_files_basename_list]\r\nprint(\"File list:\",file_list)\r\n\r\nk = tooloptions.top_k\r\n\r\n\r\noutput_root_path = \"./\"\r\n\r\n\r\nstat_first_file_path = tooloptions.stat_first_file_path\r\nstat_file_directory = os.path.dirname(stat_first_file_path)\r\n\r\nstat_file_list = tooloptions.stat_file_list\r\n\r\n\r\ntorch.cuda.empty_cache()\r\n\r\n###### read GDsearch results\r\nprint(file_list[0])\r\nstart_time = time.time()\r\n\r\n\r\n# Assuming these lists are initially empty and filled during the loop\r\nvarIDs_all = []\r\nglgm_values = []\r\np_values = []\r\n\r\nfor stat_file_name in tqdm(stat_file_list):\r\n    df = pd.read_csv(os.path.join(stat_file_directory, stat_file_name), index_col=0)\r\n    # Extend lists with values from the DataFrame\r\n    varIDs_all.extend(df.iloc[:, 0].values)  # varIDs\r\n    glgm_values.extend(df.iloc[:, 1].values)  # glgM values are in the second column\r\n    p_values.extend(df.iloc[:, 2].values)  # p-values are in the third column\r\n\r\n# Create tensors from the lists\r\nglgm_tensor = torch.tensor(glgm_values, dtype=torch.float32, device=device).view(-1, 1)\r\np_values_tensor = torch.tensor(p_values, dtype=torch.float32, device=device).view(-1, 1)\r\n# Concatenate these tensors along dimension 1 to form stat_all with shape (N, 2)\r\nstat_all = torch.cat((glgm_tensor, p_values_tensor), dim=1)\r\n\r\n\r\n##get top k variants from all chrom\r\nstart_time = time.time()\r\ntop_k_variantsID, top_k_variants_tensor, topGD_index = get_top_variants_info(file_list,\r\n                                                                                 txt_path,\r\n                                                                                 stat_all, varIDs_all, device, k=k)\r\n                                                                                 \r\nlogging.info(\"get TOP k variants from all chromosome elapsed time:{}\".format(time.time() - start_time))\r\n\r\n\r\n#  Save top k variants to CSV files\r\ndf_top_k_variantsID = pd.DataFrame(top_k_variantsID)\r\ndf_top_k_variants = pd.DataFrame(top_k_variants_tensor.to(\"cpu\"))\r\ndf_topGD_index=pd.DataFrame(topGD_index)\r\ndf_topGD_index = df_topGD_index[:1]\r\n\r\ndf_top_k_variantsID.to_csv(os.path.join(output_root_path, \"topK_pool_varIDs.csv\"), index=False)\r\ndf_top_k_variants.to_csv(os.path.join(output_root_path, \"topK_pool_variants.csv\"), index=False)\r\ndf_topGD_index.to_csv(os.path.join(output_root_path, \"topK_pool_topGD_index.csv\"), index=False)\r\n\r\n\r\n\r\n",
                    "writable": false
                },
                {
                    "entryname": "tooloptions.py",
                    "entry": "txt_first_file_path = \"$(inputs.txt_file[0].path)\"\nnpy_first_file_path = \"$(inputs.npy_file[0].path)\"\n\ntop_k = $(inputs.top_k)\n\n\nnpy_files_basename_list = ${\n          var basename = [];\n          if(inputs.npy_file){\n            inputs.npy_file.forEach(function(file){\n              basename.push('\"' + file.basename + '\"');\n            });\n          }\n          return '[' + basename.join(\",\") + ']';\n        }\n\n\nglobal_marginal_first_file_path = \"$(inputs.global_marginal_file[0].path)\"\nglobal_marginal_file_list = ${\n          var basename = [];\n          if(inputs.global_marginal_file){\n            inputs.global_marginal_file.forEach(function(file){\n              basename.push('\"' + file.basename + '\"');\n            });\n          }\n          return '[' + basename.join(\",\") + ']';\n        }",
                    "writable": false
                }
            ]
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "WorkReuse",
            "enableReuse": false
        }
    ],    
    "stdout": "standard.out",
}