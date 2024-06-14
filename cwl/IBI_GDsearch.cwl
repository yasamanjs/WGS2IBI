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
            },
            "doc": "Text files listing the subject IDs (`subIDs`) and variant IDs (`varIDs`) that correspond to the SNP data matrix, ensuring accurate data alignment."
        },
        {
            "id": "Pheno",
            "type": "File",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
            "doc": "A CSV file containing phenotypic information corresponding to the subjects in the SNP matrix.",
        },
        {
            "id": "Phenotype_name",
            "type": "string",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
            "doc": "The name of the phenotype as a string to be analyzed, guiding the analysis to focus on the specified trait."
        },
        {
            "id": "npy_file",
            "type": "File[]?",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            },
            "doc": "This input contains the preprocessed SNP data in a binary format that allows efficient storage and rapid access during analysis."
        }
    ],
    "outputs": [
        {
            "id": "GDsearch_topGD",
            "type": "File?",
            "outputBinding": {
                "glob": "*_GDsearch_topGD.csv"
            }
        },
        {
            "id": "GDsearch_stat",
            "type": "File?",
            "outputBinding": {
                "glob": "*_GDsearch_stat.csv"
            }
        }
    ],
    "doc": "### IBI_GDsearch Tool Description\n\nThe IBI_GDsearch tool is designed to perform comprehensive global Bayesian analysis at the population level, facilitating the identification of significant genetic determinants. Utilizing a compressed SNP matrix in `.npy` format alongside subject and variant IDs, phenotype data, and the specified phenotype name, this tool calculates the log marginal likelihood for each variant.  \n\n&nbsp;\n\n#### Inputs:\n- **Compressed SNP matrix (`.npy`):** This input contains the preprocessed SNP data in a binary format that allows efficient storage and rapid access during analysis.  \n- **Subject IDs and Variant IDs (`.txt`):** Text files listing the subject IDs (`subIDs`) and variant IDs (`varIDs`) that correspond to the SNP data matrix, ensuring accurate data alignment.  \n- **Phenotype Data (`.csv`):** A CSV file containing phenotypic information corresponding to the subjects in the SNP matrix.  \n- **Phenotype Name:** The name of the phenotype to be analyzed, guiding the Bayesian analysis to focus on the specified trait.  \n\n&nbsp;\n\n#### Outputs:\n- **TopGD.csv:** This file lists the top genetic determinants identified by the Bayesian analysis, providing a ranked output of the most significant variants.  \n- **IBI GDsearch Stat.csv:** This output contains detailed statistics from the Bayesian analysis, including the log marginal likelihood values for each variant.  \n\n&nbsp;\n\n#### Function:\nIBI_GDsearch employs a Bayesian framework to assess genetic variants across the entire population, calculating the log marginal likelihood for each variant. This global analysis method allows for a comprehensive evaluation of genetic determinants, providing a robust alternative or complement to traditional GWAS methods. The resulting outputs enable researchers to identify and prioritize significant variants for further study.  \n\n&nbsp;\n\nThis tool is a crucial component of the genomic analysis pipeline, enhancing the precision and depth of population-level genetic studies through advanced Bayesian methodologies.",
    "label": "IBI_GDsearch",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "LoadListingRequirement"
        },
        {
            "class": "ResourceRequirement",
            "ramMin": "${\r\n  var max_size = 0;\r\n  if (inputs.npy_file) {\r\n    inputs.npy_file.forEach(function(file) {\r\n      if(file.size > max_size) {\r\n        max_size = file.size;\r\n      }\r\n    });\r\n  }\r\n  return Math.max(max_size * 8 * 1e-6, 12000);\r\n}\r\n"
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
                    "entry": "import os\r\nimport time\r\nimport pandas as pd\r\nimport numpy as np\r\nfrom datetime import datetime\r\nimport scipy.stats as stats\r\nfrom tqdm import tqdm\r\nimport logging\r\nimport torch\r\nimport tooloptions\r\n\r\n\r\nlogging.getLogger().setLevel(logging.INFO)\r\nlogging.info(\"Libraries are loaded\")\r\n\r\n\r\n\r\ndef read_traitsF(traits_path_file, trait_string, variant_subIDs):\r\n    \"\"\"\r\n    read the .csv file with traits (subjects_row x traits_column)\r\n    :param traits_path_file:\r\n    :return:subIDs, traitIDs, traits_tensor  # list, array,tensor\r\n    \"\"\"\r\n    # make sure the subIDs become the index column; as default, column names are inferred from the first line of the\r\n    # file\r\n    traits = pd.read_csv(traits_path_file, index_col=0)\r\n    # Convert all columns to float32\r\n    traits = pd.DataFrame(traits[trait_string])\r\n    traitIDs = [trait_string]\r\n    subIDs = list(traits.index)\r\n    # traitIDs = traits.columns\r\n    traits = traits.astype('float32')\r\n    \r\n    # variant_subIDs = set(variant_subIDs)\r\n    # # Filter the traits DataFrame to only include rows with matching subIDs\r\n    # matched_subIDs = variant_subIDs.intersection(traits.index)\r\n    # traits_filtered = traits.loc[matched_subIDs]\r\n    \r\n    # subIDs = list(traits_filtered.index)\r\n    \r\n    traits_tensor = torch.as_tensor(traits.values, dtype=torch.float32)\r\n    # np.int8 has changed the type to int8; when using int8, the subIDs become negative.\r\n    #     print(np.sum(traits)) # sum gives odd results of -94.\r\n    return subIDs, traitIDs, traits_tensor  # list, array\r\n\r\n\r\n\r\ndef read_variantsF_efficient(variants_path_file, file_name):\r\n    \"\"\"\r\n    read the compressed variants data\r\n    :param variants_path_file:\r\n    :param variants_size: if variants_size=None,select all the variants\r\n    :return: subIDs, varIDs, variants_tensor, df  # list, list, array and dataframe\r\n    \"\"\"\r\n    varIDs_file = open(os.path.join(variants_path_file, \"{}_varIDs.txt\".format(file_name)), \"r\")\r\n    varIDs = varIDs_file.read().split(\"\\n\")\r\n    subIDs_file = open(os.path.join(variants_path_file, \"{}_subIDs.txt\".format(file_name)), \"r\")\r\n    subIDs = subIDs_file.read().split(\"\\n\")\r\n    subIDs = list(str(x) for x in subIDs)\r\n    variants = np.load(os.path.join(variants_path_file, \"{}_variants_value.npy\".format(file_name)))\r\n\r\n    A0 = np.ones(len(subIDs), dtype=np.int8)\r\n    variants = np.row_stack((A0, variants))\r\n    varIDs.insert(0, 'A0')\r\n    # when variants data is large, df will consume a large amount of memory,\r\n    # but df is no being used, set df=pd.Dataframe()\r\n    df = pd.DataFrame()\r\n    # df = pd.DataFrame(variants, index=varIDs, columns=subIDs, dtype=np.int8)\r\n    variants_tensor = torch.as_tensor(variants, dtype=torch.float32) #float32 for gpu and float64 for cpu\r\n\r\n    return subIDs, varIDs, variants_tensor, df\r\n\r\n\r\n\r\n\r\ndef GDsearch_all(traits_tensor, variants_tensor):\r\n    \"\"\"\r\n    Get all the stats for all the variants in any given population for multiple traits;\r\n    particulary used for the entire population\r\n\r\n    :param traits_tensor: traits n*k\r\n    :param variants_tensor: variants m*n\r\n    :return:\r\n    \"\"\"\r\n    variants_sum = variants_tensor.sum(axis=1).reshape(-1, 1)\r\n    traits_sum = traits_tensor.sum(axis=0).reshape(1, -1)\r\n    V1D1 = variants_tensor @ traits_tensor\r\n    V1D0 = torch.ones_like(V1D1) * variants_sum - V1D1\r\n    V0D1 = torch.ones_like(V1D1) * traits_sum - V1D1\r\n    V0D0 = torch.ones_like(V1D0) * variants_tensor.shape[1] - torch.ones_like(V1D0) * traits_sum - V1D0\r\n\r\n    # GiVen the Dirichlet Distributions we are using,\r\n    # the expectation of these conditional probabilities is as follows: prior probability\r\n    # P(D=1|V=1) = (alpha11 + V1D1)/(alpha1 + V1D1 + V1D0)*1.0\r\n    cp_D1V1 = (1 + V1D1) / (2 + V1D1 + V1D0) * 1.0\r\n    # P(D=1|V=0) = (alpha01 + V0D1)/(alpha0 + V0D1 + V0D0)*1.0\r\n    cp_D1V0 = (1 + V0D1) / (2 + V0D1 + V0D0) * 1.0\r\n    # RR is risk ratio; OR is oDDs ratio\r\n    RR = cp_D1V1 / cp_D1V0\r\n\r\n    # Calculate the log Marginal Likelihood for this particular SNP based on the collected counts and equation 5 in\r\n    # the worD file\r\n    # when j=0 (V=0)\r\n    lgM = torch.lgamma(torch.tensor(2.0)) - torch.lgamma(2.0 + V0D1 + V0D0)\r\n    lgM += torch.lgamma(1.0 + V0D0) - torch.lgamma(torch.tensor(1.0))\r\n    lgM += torch.lgamma(1.0 + V0D1) - torch.lgamma(torch.tensor(1.0))\r\n\r\n    # when j=1 (V=1)\r\n    lgM += torch.lgamma(torch.tensor(2.0)) - torch.lgamma(2.0 + V1D1 + V1D0)\r\n    lgM += torch.lgamma(1.0 + V1D0) - torch.lgamma(torch.tensor(1.0))\r\n    lgM += torch.lgamma(1.0 + V1D1) - torch.lgamma(torch.tensor(1.0))\r\n\r\n    if variants_tensor.shape[1] == 1:\r\n        # lgM is #traits x 1;otherwise, lgM is, variants x traits.\r\n        lgM = lgM.reshape(1, lgM.shape[0])\r\n\r\n    # get the max and index of TopGD across all the rows of variants for each column of the trait inside the 2-D array\r\n    max_value = torch.max(lgM, dim=0).values\r\n    # thus, max_value or max_index is, one vector with the size of K (# of traits)\r\n    max_index = torch.max(lgM, dim=0).indices\r\n    return RR, lgM, max_value, max_index\r\n\r\n\r\n\r\ndef save_GDsearch_result(output_root_path,traitIDs, rr, glgm, varIDs, topGD_index, glgm_topGD,file_name):\r\n    \"\"\"\r\n    collect the headers for the output file\r\n    :param traitIDs:\r\n    :return:\r\n    \"\"\"\r\n    topGD = []\r\n    for item in topGD_index:\r\n        # currently the wgs SNPs are labeled with numbers, thus varIDs and topGD both are int lists.\r\n        topGD.append(varIDs[item])\r\n    gstat_head = ['RR', 'M']\r\n    if len(traitIDs) == 1:\r\n        gstat_newhead = gstat_head\r\n    else:\r\n        gstat_newhead = []\r\n        for item in gstat_head:\r\n            for trait in traitIDs:\r\n                new = item + '_' + trait\r\n                gstat_newhead.append(new)\r\n    gstat_newhead.extend(['seq', 'varID'])\r\n\r\n    outfile_GDsearch = \"{}_GDsearch_stat.csv\".format(file_name)\r\n    # output the RR and glgm for all the variants\r\n    # output the RR and glgm for all the variants\r\n    with open(os.path.join(output_root_path,\"{}_GDsearch_stat.csv\".format(file_name)),\r\n              \"w\") as outfile:  # more efficient than using dataframe to_csv...\r\n        outfile.write(','.join(gstat_newhead) + '\\n')\r\n        for i in range(0, rr.shape[0]):\r\n            ls = []\r\n            ls.extend(rr[i].tolist())  # row i of rr that is corresponding to the ith variant\r\n            ls.extend(glgm[i].tolist())\r\n            ls.extend([str(i), varIDs[i]])\r\n            outfile.write(','.join(str(item) for item in ls) + '\\n')\r\n\r\n    with open(\r\n            os.path.join(output_root_path,\"{}_GDsearch_topGD.csv\".format(file_name)),\r\n            \"w\") as outfile:\r\n        glgm_topGD = glgm_topGD.cpu().tolist()\r\n        for i in range(0, len(traitIDs)):\r\n            line = [traitIDs[i], str(topGD[i]), str(glgm_topGD[i])]\r\n            #         print(line)\r\n            outfile.write(','.join(str(item) for item in line) + '\\n')\r\n\r\n\r\ndef GD_search_all_chromosome(file_list, txt_path, phenotype_file_path,output_root_path, device=\"cpu\"):\r\n    subIDs_BP, traitIDs, traits_tensor = read_traitsF(phenotype_file_path ,Phenotype_name, variants_subIDs)\r\n    traits_tensor = traits_tensor.to(device=device)\r\n    torch.cuda.empty_cache()\r\n    glgm_all = torch.tensor([], device=device)\r\n    varIDs_all = []\r\n    for file_name in tqdm(file_list):\r\n        start_time = time.time()\r\n        subIDs, varIDs, variants_tensor, df_variants = read_variantsF_efficient(txt_path, file_name)\r\n\r\n        logging.info(\"file_name:{}\".format(file_name))\r\n        logging.info(\"variants:{}\".format(np.shape(variants_tensor)))\r\n        logging.info(\"traits: {}\".format(np.shape(traits_tensor)))\r\n        logging.info(\"read variants data elapsed time:{}\".format(time.time() - start_time))\r\n        # 2 With GDsearch_all, calculate and output the global stats related to all the traits for all the variants\r\n        # move the tensors to GPU\r\n        variants_tensor = variants_tensor.to(device=device)\r\n        torch.cuda.empty_cache()\r\n        start_time = datetime.now()\r\n        rr, glgm, glgm_topGD, topGD_index = GDsearch_all(traits_tensor, variants_tensor)\r\n        logging.info(\"GDsearch all elapsed time: {}s \".format((datetime.now() - start_time).seconds))\r\n\r\n        # save GDsearch_results for each chrom\r\n        save_GDsearch_result(output_root_path,traitIDs, rr, glgm, varIDs, topGD_index, glgm_topGD,file_name)\r\n        glgm_all = torch.concat([glgm_all, glgm])\r\n        varIDs_all.extend(varIDs)\r\n\r\n    return glgm_all, varIDs_all\r\n\r\ndef get_top_variants_info(file_list, txt_path, glgm_all, varIDs_all, device, k=1000,\r\n                              union_flag=True):\r\n    \"\"\"\r\n    select top k variants from all chromosomes, if have multi traits, use union top k variants from all traits\r\n    :param chromosome_file_folder_path:\r\n    :param glgm_all:\r\n    :param varIDs_all:\r\n    :param device:\r\n    :param k:\r\n    :return:\r\n    \"\"\"\r\n    # Loop through each trait and take the top k variants from all chromosomes\r\n    top_k_varIDs_all = []\r\n    topGD_varID_all=[]\r\n    for traits_num in range(glgm_all.shape[1]):\r\n        glgm_all_single = glgm_all[:, traits_num]\r\n        _, sorted_indices = torch.sort(glgm_all_single, descending=True)\r\n        sorted_varIDs_all = [varIDs_all[i] for i in sorted_indices]\r\n        top_k_varIDs = sorted_varIDs_all[:k]\r\n        top_k_varIDs_all.append(top_k_varIDs)\r\n        topGD_varID = sorted_varIDs_all[0]\r\n        topGD_varID_all.append(topGD_varID)\r\n\r\n    if union_flag:\r\n        # take union top k variants from all traits\r\n        top_k_varIDs_union = list(set(x for sublist in top_k_varIDs_all for x in sublist))\r\n\r\n    top_k_varID_final = []\r\n    top_k_variants = []\r\n\r\n    # extract top k variants from all files\r\n    for file_name in tqdm(file_list):\r\n        subIDs, varIDs, variants_tensor, df_variants = read_variantsF_efficient(txt_path,\r\n                                                                                file_name)\r\n        indexs = [i for i, value in enumerate(varIDs) if value in top_k_varIDs_union]\r\n        selected_varIDs = [value for i, value in enumerate(varIDs) if value in top_k_varIDs_union]\r\n        top_k_varID_final.extend(selected_varIDs)\r\n        top_k_variants.extend(variants_tensor[indexs, :].tolist())\r\n\r\n    top_variants = np.array(top_k_variants)\r\n    if \"AO\" not in top_k_varID_final:\r\n        A0 = np.ones(top_variants.shape[1], dtype=np.int8)\r\n        top_variants = np.row_stack((A0, top_variants))\r\n        top_k_varID_final.insert(0, 'A0')\r\n    top_variants_tensor = torch.tensor(top_variants, dtype=torch.float32, device=device)\r\n\r\n    # get topGD_index for all traits\r\n    topGD_index = [top_k_varID_final.index(topGD_varID) for topGD_varID in topGD_varID_all]\r\n    return top_k_varID_final, top_variants_tensor, topGD_index\r\n\r\n\r\n\r\n\r\n\r\n#######################################################################################\r\n#######################################################################################\r\n\r\ndevice = torch.device('cuda' if torch.cuda.is_available() else 'cpu')\r\n\r\n\r\n################# main ############################################\r\nPhenotype_name = tooloptions.Phenotype_name\r\nphenotype_file_path = tooloptions.gzpath_pheno\r\n\r\n\r\n## Get info for Phenotype and SNP data and etc\r\ntxt_first_file_path = tooloptions.txt_first_file_path\r\nnpy_first_file_path = tooloptions.npy_first_file_path\r\n\r\ntxt_path = os.path.dirname(txt_first_file_path)\r\n\r\n\r\nnpy_files_basename_list = tooloptions.npy_files_basename_list\r\n# Extract the file list without the '_variants_value.npy' part\r\nfile_list = [filename.split('_variants_value.npy')[0] for filename in npy_files_basename_list]\r\nprint(\"File list:\",file_list)\r\n\r\n\r\n\r\n\r\nvariants_subIDs_file = open(os.path.join(os.path.dirname(txt_first_file_path), \"{}_subIDs.txt\".format(file_list[0])), \"r\")\r\nvariants_subIDs = variants_subIDs_file.read().split(\"\\n\")\r\nvariants_subIDs = list(str(x) for x in variants_subIDs)\r\n\r\n\r\noutput_root_path = \"./\"\r\n\r\n\r\n\r\n\r\n## 2 GDsearch\r\ntorch.cuda.empty_cache()\r\n\r\nprint(file_list[0])\r\n\r\n\r\n\r\n# Extract subIDs from traits data\r\ntraits_subIDs, traitIDs , traits_tensor = read_traitsF(phenotype_file_path, Phenotype_name, variants_subIDs)\r\n\r\n\r\nstart_time = time.time()\r\n\r\n\r\n\r\nglgm_all, varIDs_all = GD_search_all_chromosome(file_list, txt_path, phenotype_file_path,output_root_path,device=device)\r\n\r\n## save GDsearch all results\r\nprint(\"get GDsearch results from all chromosome elapsed time:{}\".format((time.time() - start_time)))\r\nlogging.info(\"get GDsearch results from all chromosome elapsed time:{}\".format((time.time() - start_time)))\r\n\r\n\r\n# # Convert the tensor to a numpy array\r\n# glgm_all_np = glgm_all.cpu().detach().numpy()\r\n# # Create a DataFrame from the numpy array\r\n# glgm_all_df = pd.DataFrame(glgm_all_np)\r\n# glgm_all_df['varIDs'] = pd.Series(varIDs_all)\r\n# # Save the DataFrame to a CSV file ## extra\r\n# glgm_all_df.to_csv(\"glgm_all.csv\", index=False)\r\n# print(glgm_all)\r\n\r\n\r\n",
                    "writable": false
                },
                {
                    "entryname": "tooloptions.py",
                    "entry": "txt_first_file_path = \"$(inputs.txt_file[0].path)\"\nnpy_first_file_path = \"$(inputs.npy_file[0].path)\"\n\ngzpath_pheno = \"$(inputs.Pheno.path)\"\nPhenotype_name = \"$(inputs.Phenotype_name)\"\n\nnpy_files_basename_list = ${\n          var basename = [];\n          if(inputs.npy_file){\n            inputs.npy_file.forEach(function(file){\n              basename.push('\"' + file.basename + '\"');\n            });\n          }\n          return '[' + basename.join(\",\") + ']';\n        }",
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