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
            "id": "Pheno",
            "type": "File",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            }
        },
        {
            "id": "Phenotype_name",
            "type": "string",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            }
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
            "id": "topK_pool_varIDs",
            "type": "File",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            }
        },
        {
            "id": "topK_pool_variants",
            "type": "File",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            }
        },
        {
            "id": "topK_pool_topGD_index",
            "type": "File",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            }
        }
    ],
    "outputs": [
        {
            "id": "IBI_sGD_stat",
            "type": "File?",
            "outputBinding": {
                "glob": "*_sGD_stat.csv"
            }
        }
    ],
    "doc": "### IBI_sGD \n \n\nThe IBI_sGD tool is designed to perform secondary genetic analysis at the individual chromosome level, leveraging the pool of top variants identified in the earlier stages of the analysis. This tool integrates the top variants with the full genomic data for each chromosome, enabling a detailed and focused investigation of genetic determinants.\n\n&nbsp;\n\n#### Inputs:\n- **SNP matrix (`.npy`):** A binary format file containing the preprocessed SNP data for each chromosome.\n- **Subject IDs (`subIDs.txt`):** A text file listing the subject IDs corresponding to the SNP data matrix.\n- **Variant IDs (`varIDs.txt`):** A text file listing the variant IDs corresponding to the SNP data matrix.\n- **Pool varIDs (`pool_varIDs.csv`):** A CSV file listing the variant IDs in the top K pool, created by the Confirm_topk_Pool tool.\n- **Pool variants matrix (`pool_variants_matrix.csv`):** A CSV file providing a matrix of the pooled top variants.\n- **TopGD index (`pool_topGD_index.csv`):** A CSV file detailing the index of the top genetic variant in the pool.\n\n\n&nbsp;\n\n#### Outputs:\n- **IBI_sGD statistics (`IBI_sGD_stats.csv`):** A CSV file containing the results of the secondary genetic analysis, including the statistical significance and other relevant metrics for each variant.\n\n&nbsp;\n\n#### Function:\nThe IBI_sGD tool performs several key functions to conduct a secondary genetic analysis:\n\n1. **Integration of Top Variants:**\n   - The tool integrates the pool of top variants with the SNP matrix for each chromosome. This integration allows for a comprehensive analysis that includes both the top-ranked variants and the full set of variants from the chromosome.\n\n2. **Secondary Genetic Analysis:**\n   - The tool applies the Individualized Bayesian Inference (IBI) method to perform a detailed analysis of genetic variants, assessing their associations with the specified phenotypes. This method is particularly effective for identifying subtle genetic signals that might be missed by other approaches.\n\n3. **Output Generation:**\n   - The tool generates a `IBI_sGD_results.csv` file that contains the results of the analysis. This file includes metrics such as the log marginal likelihood of each variant, providing researchers with detailed insights into the genetic determinants of the phenotypes under study.\n\n\n&nbsp;\n\nThe IBI_sGD tool is a critical component of the genomic analysis pipeline, enabling researchers to conduct a detailed and focused secondary analysis of genetic variants. By integrating top variants with full chromosome data, the tool enhances the precision and depth of the overall genomic investigation.",
    "label": "IBI_sGD_CPU",
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
                    "entry": "import os\r\nimport time\r\nimport pandas as pd\r\nimport numpy as np\r\nfrom datetime import datetime\r\nimport scipy.stats as stats\r\nfrom tqdm import tqdm\r\nimport logging\r\nimport torch\r\nimport tooloptions\r\n\r\n\r\nlogging.getLogger().setLevel(logging.INFO)\r\nlogging.info(\"Libraries are loaded\")\r\n\r\n\r\n\r\ndef read_traitsF(traits_path_file, trait_string):\r\n    \"\"\"\r\n    read the .csv file with traits (subjects_row x traits_column)\r\n    :param traits_path_file:\r\n    :return:subIDs, traitIDs, traits_tensor  # list, array,tensor\r\n    \"\"\"\r\n    # make sure the subIDs become the index column; as default, column names are inferred from the first line of the\r\n    # file\r\n    traits = pd.read_csv(traits_path_file, index_col=0)\r\n    # Convert all columns to float32\r\n    traits = pd.DataFrame(traits[trait_string])\r\n    traitIDs = [trait_string]\r\n    subIDs = list(traits.index)\r\n    # traitIDs = traits.columns\r\n    traits = traits.astype('float32')\r\n    \r\n    # variant_subIDs = set(variant_subIDs)\r\n    # # Filter the traits DataFrame to only include rows with matching subIDs\r\n    # matched_subIDs = variant_subIDs.intersection(traits.index)\r\n    # traits_filtered = traits.loc[matched_subIDs]\r\n    \r\n    # subIDs = list(traits_filtered.index)\r\n    \r\n    traits_tensor = torch.as_tensor(traits.values, dtype=torch.float32)\r\n    # np.int8 has changed the type to int8; when using int8, the subIDs become negative.\r\n    #     print(np.sum(traits)) # sum gives odd results of -94.\r\n    return subIDs, traitIDs, traits_tensor  # list, array\r\n\r\n\r\n\r\ndef read_variantsF_efficient(variants_path_file, file_name):\r\n    \"\"\"\r\n    read the compressed variants data\r\n    :param variants_path_file:\r\n    :param variants_size: if variants_size=None,select all the variants\r\n    :return: subIDs, varIDs, variants_tensor, df  # list, list, array and dataframe\r\n    \"\"\"\r\n    varIDs_file = open(os.path.join(variants_path_file, \"{}_varIDs.txt\".format(file_name)), \"r\")\r\n    varIDs = varIDs_file.read().split(\"\\n\")\r\n    subIDs_file = open(os.path.join(variants_path_file, \"{}_subIDs.txt\".format(file_name)), \"r\")\r\n    subIDs = subIDs_file.read().split(\"\\n\")\r\n    subIDs = list(str(x) for x in subIDs)\r\n    variants = np.load(os.path.join(variants_path_file, \"{}_variants_value.npy\".format(file_name)))\r\n\r\n    A0 = np.ones(len(subIDs), dtype=np.int8)\r\n    variants = np.row_stack((A0, variants))\r\n    varIDs.insert(0, 'A0')\r\n    # when variants data is large, df will consume a large amount of memory,\r\n    # but df is no being used, set df=pd.Dataframe()\r\n    df = pd.DataFrame()\r\n    # df = pd.DataFrame(variants, index=varIDs, columns=subIDs, dtype=np.int8)\r\n    variants_tensor = torch.as_tensor(variants, dtype=torch.float32) #float32 for gpu and float64 for cpu\r\n\r\n    return subIDs, varIDs, variants_tensor, df\r\n\r\n\r\n\r\n\r\ndef cal_lgM(top_k_variants, variants_batch, weights_batch, traits, traits_xor, variant_type):\r\n    \"\"\"\r\n    Calculate and return the lgM for all the drivers or any driver for any given population for multiple traits\r\n    :param variants:\r\n    :param variants_batch:\r\n    :param weights_batch:\r\n    :param traits:\r\n    :param traits_xor:\r\n    :param variant_type:\r\n    :return:\r\n    \"\"\"\r\n    # BP_V is DW, BP_V_xor is (J-D)W', traits_xor is J-traits, J is all one matrix\r\n    # BP_V represents the weighted traits of all variants, weights_batch is b*n, traits is n*k,BP_V is b*n*k,\r\n    # b is the batch_size\r\n    BP_V = torch.einsum('ij,jk->ijk', weights_batch, traits)\r\n    # BP_V_xor represents the weighted traits_xor of all variants, weights_batch is b*n, traits_xor is n*k,\r\n    # BP_V_xor is b*n*k, b is the batch_size\r\n    BP_V_xor = torch.einsum('ij,jk->ijk', weights_batch, traits_xor)\r\n    # calculate the second dimension summation of BP_V, and add a broadcasting dimension, BP_V_sum is b*1*k\r\n    BP_V_sum = BP_V.sum(dim=1).unsqueeze(dim=1)\r\n    # calculate the second dimension summation of BP_V_xor, and add a broadcasting dimension, BP_V_xor_sum is b*1*k\r\n    BP_V_xor_sum = BP_V_xor.sum(dim=1).unsqueeze(dim=1)\r\n    # variant=1\r\n    if variant_type == 1:\r\n        # Unsqueeze variants_batch to add a broadcasting dimension: from (b, n) to (b, 1, n)\r\n        variants_unsqueezed = variants_batch.unsqueeze(dim=1)\r\n        # three-dimensional matrix multiplication, V1D1, V1D0, V0D1, V0D0 is b*1*k, for all the snps\r\n        # V1D1 = torch.bmm(variants_unsqueezed, BP_V)\r\n        # V1D0 = torch.bmm(variants_unsqueezed, BP_V_xor)\r\n        V1D1 = torch.einsum(\"bqn,bnk->bqk\",variants_unsqueezed, BP_V)\r\n        V1D0 = torch.einsum(\"bqn,bnk->bqk\",variants_unsqueezed, BP_V_xor)\r\n        V0D1 = BP_V_sum - V1D1\r\n        V0D0 = BP_V_xor_sum - V1D0\r\n    # variant=0\r\n    else:\r\n        # top_k_variants t*n\r\n        V1D1 = torch.einsum('tn,bnk->btk', top_k_variants, BP_V)\r\n        V1D0 = torch.einsum('tn,bnk->btk', top_k_variants, BP_V_xor)\r\n        V0D1 = torch.ones_like(V1D1) * BP_V_sum - V1D1\r\n        V0D0 = torch.ones_like(V1D0) * BP_V_xor_sum - V1D0\r\n\r\n    # if variant_type=1, lgM is b*1*k; if variant_type=0, lgM is b*m*k\r\n    lgM = torch.lgamma(torch.tensor(2.0)) - torch.lgamma(2.0 + V0D1 + V0D0)\r\n    lgM += torch.lgamma(1.0 + V0D0) - torch.lgamma(torch.tensor(1.0))\r\n    lgM += torch.lgamma(1.0 + V0D1) - torch.lgamma(torch.tensor(1.0))\r\n\r\n    # when j=1 (V=1)\r\n    lgM += torch.lgamma(torch.tensor(2.0)) - torch.lgamma(2.0 + V1D1 + V1D0)\r\n    lgM += torch.lgamma(1.0 + V1D0) - torch.lgamma(torch.tensor(1.0))\r\n    lgM += torch.lgamma(1.0 + V1D1) - torch.lgamma(torch.tensor(1.0))\r\n\r\n\r\n    if variants_tensor.ndim == 1:\r\n        lgM = lgM.reshape(1, lgM.shape[0])\r\n\r\n    return lgM\r\n\r\n\r\ndef lgMcal(variants, traits, device,varIDs,top_k_variantsID,top_k_variants,use_oneTopGD):\r\n    \"\"\"\r\n    sGDsearch core code\r\n    :param variants: m*n, m:the number of variants, n:the number of individuals\r\n    :param traits:  n*k, n:the number of individuals, k:the number of traits\r\n    :param devices:\r\n    :return:\r\n    \"\"\"\r\n    # set batch size\r\n    batch_size = 250\r\n    variants_num = variants.shape[0]\r\n    split_size = int(np.ceil(variants_num / batch_size))\r\n    logging.info(\"batch_size:{},split_size:{}\".format(batch_size, split_size))\r\n    offset = 0\r\n    # element_run save sGD results\r\n    element_run = []\r\n    # cycle through each batch\r\n    for i in tqdm(range(split_size)):\r\n        # get variants and varIDs for each batch, variants_batch is b*n, b is the batch_size\r\n        variants_batch = variants[offset:offset + batch_size]\r\n        varIDs_batch = np.array(varIDs[offset:offset + batch_size]).reshape(-1, 1)\r\n        # traits_xor is the inverse of traits, the 0 in traits become 1 in traits_xor, and\r\n        # 1 in traits become 0 in traits_xor\r\n        traits_xor = torch.ones(traits.shape, dtype=torch.float32, device=device) - traits\r\n\r\n        # variant=1\r\n        weights_batch_v1 = variants_batch\r\n        lgMv1_SD = cal_lgM(top_k_variants, variants_batch, weights_batch_v1, traits, traits_xor,variant_type=1)\r\n        # remove one dimension,from (b,1,k) to (b,k)\r\n        lgMv1_SD = lgMv1_SD.squeeze(dim=1)\r\n        # released unused variable\r\n        del weights_batch_v1\r\n        torch.cuda.empty_cache()\r\n\r\n        # variants=0\r\n        # select individuals with variants=0 for all the snps\r\n        weights_batch_v0 = torch.ones(variants_batch.shape, dtype=torch.float32, device=device) - variants_batch\r\n        lgMv0 = cal_lgM(top_k_variants, variants_batch, weights_batch_v0, traits, traits_xor,variant_type=0)\r\n        # released unused variable\r\n        del weights_batch_v0\r\n        torch.cuda.empty_cache()\r\n        # collect the lgMv0_topGD for each trait; the lgM value for V0 group when using topGD as the driver\r\n        lgMv0_topGD = []\r\n        # collect the r between SD and topGD for each trait\r\n        r = []\r\n        variants_batch_array = variants_batch.to(\"cpu\").numpy()\r\n        if use_oneTopGD:\r\n            for m in range(0, len(traitIDs)):\r\n                lgMv0_topGD.append(list(lgMv0[:, m, m].cpu().numpy()))\r\n            for j in topGD_index:\r\n                r_each_trait = []\r\n                for row in variants_batch_array:\r\n                    # [0] to get only the coefficient and ignore the p-values\r\n                    r1 = stats.spearmanr(row, variants[j, :].to(\"cpu\").numpy())[0]\r\n                    r_each_trait.append(r1)  # a vector of K\r\n                r.append(r_each_trait)\r\n            lgMv0_sGD = torch.zeros((variants_batch.shape[0], len(traitIDs)), device=device)\r\n            sGD = np.zeros((variants_batch.shape[0], len(traitIDs)))\r\n        else:\r\n            # with sGD, lgMv0 is bath_size*m_variants*k_traits\r\n            lgMv0_sGD = torch.max(lgMv0, dim=1).values\r\n            sGD_index = torch.max(lgMv0, dim=1).indices\r\n\r\n            # collect the variant ID of sGD for each batch and each trait in a 1D array\r\n            sGD = np.array([top_k_variantsID[i] for pair in sGD_index for i in pair]).reshape(sGD_index.shape)\r\n\r\n            k = 0\r\n            # collect the lgMv0_topGD and r for each trait in a 1D array specifically with mxk lgMv0\r\n            # topGD_index is one output from GDsearch_all, a vector of K (#traits ordered in the original trait input file)\r\n            for j in topGD_index:\r\n                # a vector of K\r\n                lgMv0_topGD.append(list(lgMv0[:, j, k].cpu().numpy()))\r\n\r\n                r_each_trait = []\r\n                for row in variants_batch_array:\r\n                    # [0] to get only the coefficient and ignore the p-values\r\n                    r1 = stats.spearmanr(row, top_k_variants[j, :].to(\"cpu\").numpy())[0]\r\n                    r_each_trait.append(r1)  # a vector of K\r\n                r.append(r_each_trait)\r\n                k = k + 1\r\n        lgMv0_topGD = torch.tensor(np.array(lgMv0_topGD).T, device=device)\r\n        r = np.array(r).T\r\n\r\n        if use_oneTopGD:\r\n            lgM_v1v0 = lgMv1_SD + lgMv0_topGD\r\n        else:\r\n            lgM_v1v0 = lgMv1_SD + lgMv0_sGD\r\n        # get sequence number for varIDs_batch\r\n        seq = np.arange(offset, offset + variants_batch.shape[0]).reshape(-1, 1)\r\n\r\n        # save all the batch results to merged_arr\r\n        merged_arr = np.concatenate([lgMv1_SD.cpu().numpy(), lgMv0_sGD.cpu().numpy()], axis=1)\r\n        merged_arr = np.concatenate([merged_arr, lgMv0_topGD.cpu().numpy()], axis=1)\r\n        merged_arr = np.concatenate([merged_arr, lgM_v1v0.cpu().numpy()], axis=1)\r\n        merged_arr = np.concatenate([merged_arr, sGD], axis=1)\r\n        merged_arr = np.concatenate([merged_arr, r], axis=1)\r\n        merged_arr = np.concatenate([merged_arr, seq], axis=1)\r\n        merged_arr = np.concatenate([merged_arr, varIDs_batch], axis=1)\r\n\r\n        element_run.extend(merged_arr)\r\n        offset += batch_size\r\n    return element_run\r\n\r\n\r\n\r\ndef save_sGD_result(element_run,file_name,output_sGD_root_path):\r\n    \"\"\"\r\n    collect the headers for this file\r\n    :return:\r\n    \"\"\"\r\n    # return(lgMv1_SD, lgMv0_sGD, lgMv0_topGD, lgM_v1v0, sGD, r, i, varID)\r\n    outlgM = ['lgMv1_SD', 'lgMv0_sGD', 'lgMv0_topGD', 'lgM_v1v0', 'sGD', 'r']\r\n    if len(traitIDs) == 1:\r\n        outAll = outlgM\r\n    else:\r\n        outAll = []\r\n        for item in outlgM:\r\n            for trait in traitIDs:\r\n                new = item + '_' + trait\r\n                outAll.append(new)\r\n    outAll = outAll + ['seq', 'varID']\r\n    df_res = pd.DataFrame(element_run, columns=outAll)\r\n    output_file_path = os.path.join(output_sGD_root_path, \"{}_sGD_stat.csv\".format(file_name.split(\".\")[0]))\r\n    df_res.to_csv(output_file_path, index=False)\r\n\r\n\r\n\r\n\r\n\r\n#######################################################################################\r\n#######################################################################################\r\n\r\ndevice = torch.device('cuda' if torch.cuda.is_available() else 'cpu')\r\n\r\n\r\n################# main ############################################\r\nPhenotype_name = tooloptions.Phenotype_name\r\nphenotype_file_path = tooloptions.gzpath_pheno\r\n\r\n\r\n## Get info for Phenotype and SNP data and etc\r\ntxt_first_file_path = tooloptions.txt_first_file_path\r\nnpy_first_file_path = tooloptions.npy_first_file_path\r\n\r\ntxt_path = os.path.dirname(txt_first_file_path)\r\n\r\n\r\nnpy_files_basename_list = tooloptions.npy_files_basename_list\r\n# Extract the file list without the '_variants_value.npy' part\r\nfile_list = [filename.split('_variants_value.npy')[0] for filename in npy_files_basename_list]\r\nprint(\"File list:\",file_list)\r\n\r\n\r\noutput_root_path = \"./\"\r\n\r\n\r\nGDsearch_pool_varIDs_path = tooloptions.GDsearch_pool_varIDs_path\r\nGDsearch_pool_variants_path = tooloptions.GDsearch_pool_variants_path\r\nGDsearch_pool_topGD_index_path = tooloptions.GDsearch_pool_topGD_index\r\n\r\n\r\ndf_top_k_variantsID = pd.read_csv(GDsearch_pool_varIDs_path)\r\ndf_top_k_variants = pd.read_csv(GDsearch_pool_variants_path)\r\ndf_topGD_index = pd.read_csv(GDsearch_pool_topGD_index_path)\r\n\r\n\r\n\r\ntop_k_variantsID = df_top_k_variantsID.values.flatten().tolist()\r\ntop_k_variants_tensor = torch.as_tensor(df_top_k_variants.values, dtype=torch.float32, device=device)\r\ntopGD_index=df_topGD_index.values.flatten().tolist()\r\n\r\n\r\n\r\n    # 3 sGD search\r\n    # An important flag to dictate whether using topGD or sGD as the driver for A0 group.\r\n    \r\nuse_oneTopGD = False\r\n\r\n\r\ntraits_subIDs, traitIDs, traits_tensor = read_traitsF(phenotype_file_path, Phenotype_name)\r\ntraits_tensor = traits_tensor.to(device=device)\r\nstart_time=time.time()\r\nfor file_name in file_list:\r\n    start_time_each=time.time()\r\n    subIDs, varIDs, variants_tensor, df_variants = read_variantsF_efficient(txt_path,file_name)\r\n    variants_tensor = variants_tensor.to(device=device)\r\n    logging.info(\"device:{}\".format(device))\r\n\r\n    element_run = lgMcal(variants_tensor, traits_tensor, device,varIDs,top_k_variantsID,top_k_variants_tensor,use_oneTopGD)\r\n    save_sGD_result(element_run,file_name,output_root_path)\r\n    logging.info(\"file_name:{},elapsed time:{}\".format(file_name,time.time()-start_time_each))\r\n\r\n    logging.info(\"use TOP k variants from all chromosome for sGD elapsed time:{}\".format(time.time() - start_time))",
                    "writable": false
                },
                {
                    "entryname": "tooloptions.py",
                    "entry": "txt_first_file_path = \"$(inputs.txt_file[0].path)\"\nnpy_first_file_path = \"$(inputs.npy_file[0].path)\"\n\ngzpath_pheno = \"$(inputs.Pheno.path)\"\nPhenotype_name = \"$(inputs.Phenotype_name)\"\n\n\nnpy_files_basename_list = ${\n          var basename = [];\n          if(inputs.npy_file){\n            inputs.npy_file.forEach(function(file){\n              basename.push('\"' + file.basename + '\"');\n            });\n          }\n          return '[' + basename.join(\",\") + ']';\n        }\n        \n\nGDsearch_pool_varIDs_path = \"$(inputs.topK_pool_varIDs.path)\"\nGDsearch_pool_variants_path = \"$(inputs.topK_pool_variants.path)\"\nGDsearch_pool_topGD_index = \"$(inputs.topK_pool_topGD_index.path)\"",
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