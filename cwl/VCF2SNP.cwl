{
    "class": "Workflow",
    "cwlVersion": "v1.2",
    "doc": "### VCF2SNP Tool Description\n\nThe VCF2SNP tool provides a comprehensive preprocessing pipeline for whole genome sequencing (WGS) data, converting it into a format ready for downstream analysis using IBI or Fisher methods. This tool integrates several key components to ensure that the genomic data is meticulously processed and transformed into a clean, structured SNP matrix.\n\n&nbsp;\n\n#### Inputs:\n- **WGS data (VCF.gz):** The initial raw genomic data in VCF.gz format, containing variant information from whole genome sequencing.\n\n&nbsp;\n\n#### Outputs:\n- **Compressed SNP matrix (`.npy`):** This output contains the preprocessed SNP data in a binary format that allows efficient storage and rapid access during analysis.\n- **Subject IDs and Variant IDs (`subIDs.txt` and `varIDs.txt`):** Text files listing the subject IDs and variant IDs that correspond to the SNP data matrix, ensuring accurate data alignment.\n- **Annotations (`annotations.csv`):** A CSV file containing additional annotations for the SNPs, providing contextual information for downstream analyses.\n\n&nbsp;\n\n#### Function:\nThe VCF2SNP tool performs a series of crucial steps to preprocess the WGS data:\n\n1. **BCFtools Merge and Filter:**\n   - This component merges multiple VCF or BCF files from non-overlapping sample sets into one multi-sample file. It also filters out any monomorphic variants, ensuring that only polymorphic markers are retained. Developed by the Genetic Analysis Center (GAC) at the University of Washington, BCFtools ensures consistency and completeness in the merged dataset.\n   - **Output:** A single merged VCF.gz WGS file.\n\n2. **PLINK2:**\n   - Developed by Shaun Purcell at the Broad Institute and Harvard Medical School, PLINK2 is a widely accepted tool in the genomic research community. In this workflow, a custom app based on PLINK2 simplifies the toolset to include only the necessary functions, ensuring ease of use. PLINK2 refines the merged data by applying several filters, including Minor Allele Frequency (MAF), Hardy-Weinberg Equilibrium (HWE), genotype missingness (Geno), sample missingness (Mind), and quality metrics (Qual). It matches IDs with phenotypes using the `--keep` option, enhancing the quality of the dataset for further processing.\n   - **Output:** A filtered WGS dataset in VCF.gz format.\n\n3. **input_prep_python:**\n   - This component processes the filtered VCF data into a clean SNP*Sample matrix, applying dominant encoding to simplify the representation of genotype data. It separates the SNP matrix from annotations and compresses the data into formats tailored for downstream analyses.\n   - **Output:** A compressed SNP matrix in `.npy` format, text files containing subject and variant IDs (`subIDs.txt` and `varIDs.txt`), and a CSV file with additional annotations.\n\n&nbsp;\n\nBy transforming raw genomic data into standardized, analysis-ready formats, the VCF2SNP tool serves as the foundational step in the genomic analysis pipeline. It ensures that the resulting datasets are primed for high-precision analysis required in identifying significant genetic markers associated with various diseases.",
    "label": "VCF2SNP",    
    "inputs": [
        {
            "id": "in_variants",            
            "type": "File[]",
            "label": "Variants VCF.GZ files",
            "doc": "Input files which will be merged.",
            "secondaryFiles": [
                {
                    "pattern": "${ \n    if(self.basename.split('.').pop() == 'gz'){\n        if(self.nameroot.split('.').pop() == 'bcf'){\n            return self.nameroot + \".gz.csi\"}\n        else{\n            return self.nameroot + \".gz.tbi\"\n    }\n}  else{\n        if(self.basename.split('.').pop() == 'bcf'){\n            return self.basename + \".csi\"\n    }\n        else{\n            return self.basename + \".tbi\"}\n}\n\n}",
                    "required": true
                }
            ],
            
        },
        {
            "id": "Missingness_per_individual",
            "type": "float?",
            "doc": "--mind\n\nexclude individuals with too much missing genotype data.",            
        },
        {
            "id": "Missingness_per_marker",
            "type": "float?",
            "doc": "--geno\n\nexclude SNPs on the basis of missing genotype rate, with the --geno option: the default is to include all SNPS (i.e. --geno 1). To include only SNPs with a 90% genotyping rate (10% missing) use --geno 0.1",            
        },
        {
            "id": "Hardy_Weinberg_equilibrium",
            "type": "float?",
            "doc": "--hwe\n\nTo exclude markers that failure the Hardy-Weinberg test at a specified significance threshold",            
        },
        {
            "id": "Mendel_error_rates",
            "type": "float[]?",
            "doc": "example:\nplink --file mydata --me 0.05 0.1\n\nFor family-based data only, to exclude individuals and/or markers on the basis on Mendel error rate.\n\n\n-  the first parameter determines that families with more than 5% Mendel errors (considering all SNPs) will be discarded.\n-  the second parameter indicates that SNPs with more than 10% Mendel error rate will be excluded (i.e. based on the number of trios);",            
        },
        {
            "id": "Biallelic",
            "type": "boolean?",
            "doc": "--max-alleles 2 filters out the multiallelic variants",            
        },
        {
            "id": "SNPs_only",
            "type": "boolean?",
            "doc": "Set TRUE if you want to include this filter:\n\n--snps-only excludes all variants with one or more multi-character allele codes. With 'just-acgt', variants with single-character allele codes outside of {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', <missing code>} are also excluded.",            
        },
        {
            "id": "MAF_Threshold",
            "type": "float?",
            "doc": "--maf\n\nexclude SNPs on the basis of MAF (minor allele frequency)",            
        },
        {
            "id": "QUAL",
            "type": "int?",
            "doc": "--var-min-qual causes all variants with QUAL value smaller than the given number, or with no QUAL value at all, to be skipped.",            
        }
    ],
    "outputs": [
        {
            "id": "annotation_csv",
            "outputSource": [
                "datafiltering2/annotation_csv"
            ],
            "type": "File[]?",
            "label": "Annotations CSV",            
        },
        {
            "id": "compressed_SNP_info",
            "outputSource": [
                "datafiltering2/compressed_SNP_info"
            ],
            "type": "File[]?",
            "label": "compressed_SNP_info",            
        },
        {
            "id": "out_variants",
            "outputSource": [
                "bcftools_merge_and_filter/out_variants"
            ],            
            "type": "File[]?",
            "label": "Merged VCF.GZ file",
            "doc": "Merged output file",            
        }
    ],
    "steps": [
        {
            "id": "bcftools_merge_and_filter",
            "in": [
                {
                    "id": "in_variants",
                    "source": [
                        "in_variants"
                    ]
                }
            ],
            "out": [
                {
                    "id": "out_variants"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",                
                "baseCommand": [],
                "inputs": [
                    {                        
                        "id": "in_variants",
                        "type": "File[]",
                        "label": "Input variants files",
                        "doc": "Input files which will be merged.",                        
                        "secondaryFiles": [
                            {
                                "pattern": "${ \n    if(self.basename.split('.').pop() == 'gz'){\n        if(self.nameroot.split('.').pop() == 'bcf'){\n            return self.nameroot + \".gz.csi\"}\n        else{\n            return self.nameroot + \".gz.tbi\"\n    }\n}  else{\n        if(self.basename.split('.').pop() == 'bcf'){\n            return self.basename + \".csi\"\n    }\n        else{\n            return self.basename + \".tbi\"}\n}\n\n}",
                                "required": true
                            }
                        ]
                    },
                    {
                        
                        "id": "output_type",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "CompressedBCF",
                                    "UncompressedBCF",
                                    "CompressedVCF",
                                    "UncompressedVCF"
                                ],
                                "name": "output_type"
                            }
                        ],
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 90,
                            "valueFrom": "${\n    \n    //if(inputs.in_variants.length == 1){\n    //    return \"\"\n    //}else{\n        \n        if (inputs.output_type === 'CompressedBCF') return '--output-type b'\n        if (inputs.output_type === 'UncompressedBCF') return '--output-type u'\n        if (inputs.output_type === 'CompressedVCF') return '--output-type z'\n        if (inputs.output_type === 'UncompressedVCF') return '--output-type v'\n    \n    //}\n}\n    \n"
                        },
                        "label": "Output type",
                        "doc": "Output types: b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v].",
                        "default": "CompressedVCF"
                    },
                    {
                        
                        "id": "regions",
                        "type": "string[]?",
                        "inputBinding": {
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 18,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        return \"--regions \"+ inputs.regions\n    }\n}"
                        },
                        "label": "Regions for processing",
                        "doc": "Restrict to comma-separated list of regions."
                    },
                    {
                        
                        "id": "regions_file",
                        "type": "File?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 19,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        return \"--regions-file \" + inputs.regions_file\n    }\n}"
                        },
                        "label": "Regions file",
                        "doc": "Restrict to regions listed in a file.",
                        
                    },
                    {
                        
                        "id": "threads",
                        "type": "int?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 95,
                            "valueFrom": "${\n    //if(inputs.in_variants.length == 1){\n    //    return \"\"\n    //}else{\n        return \"--threads \" + inputs.threads\n    //}\n}"
                        },
                        "label": "Threads",
                        "doc": "Number of extra output compression threads [0]."
                    },
                    {
                        
                        "id": "tag",
                        "type": "string?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 28,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        return \"--tag \" + inputs.tag\n    }\n}"
                        },
                        "label": "Tagged GEN file",
                        "doc": "Tag to take values for .gen file: GT,PL,GL,GP."
                    },
                    {                        
                        "id": "mem_per_job",
                        "type": "int?",
                        "label": "Memory per job",
                        "doc": "Memory per job in MB. Appropriate instance will be chosen based on this parameter."
                    },
                    {
                        
                        "id": "cpu_per_job",
                        "type": "int?",
                        "label": "CPU per job",
                        "doc": "Number of CPUs per job. Appropriate instance will be chosen based on this parameter."
                    },
                    {                        
                        "id": "force_samples",
                        "type": "boolean?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 2,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        \n        if(inputs.force_samples == true){\n        return \"--force-samples\"\n        }else{\n            return \"\"\n        }\n    }\n}"
                        },
                        "label": "Force samples",
                        "doc": "If the merged files contain duplicate samples names, proceed anyway. Duplicate sample names will be resolved by prepending index of the file as it appeared on the command line to the conflicting sample name."
                    },
                    {                        
                        "id": "use_header",
                        "type": "File?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        return \"--use-header \" + inputs.use_header\n    }\n}"
                        },
                        "label": "Use provided header",
                        "doc": "Use the provided header.",                        
                    },
                    {                        
                        "id": "missing_genotypes",
                        "type": "boolean?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 5,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        if(inputs.missing_genotypes == true){\n        return \"--missing-to-ref\"\n        }else{\n            return \"\"\n        }\n        \n    }\n}"
                        },
                        "label": "Assume missing genotypes",
                        "doc": "Assume genotypes at missing sites are 0/0."
                    },
                    {                        
                        "id": "apply_filters",
                        "type": "string[]?",
                        "inputBinding": {
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 6,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        return \"--apply-filters \" + inputs.apply_filters\n    }\n}"
                        },
                        "label": "Apply filters",
                        "doc": "Require at least one of the listed FILTER strings (e.g. \"PASS,.\")."
                    },
                    {                        
                        "id": "filter_logic",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "x",
                                    "+"
                                ],
                                "name": "filter_logic"
                            }
                        ],
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 7,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        return \"--filter-logic \" + inputs.filter_logic\n    }\n}"
                        },
                        "label": "Filters logic",
                        "doc": "Remove filters if some input is PASS (\"x\"), or apply all filters (\"+\")."
                    },
                    {                       
                        "id": "info_rules",
                        "type": "string[]?",
                        "inputBinding": {
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 9,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        return \"--info-rules \" + inputs.info_rules\n    }\n}"
                        },
                        "label": "Info rules",
                        "doc": "Rules for merging INFO fields (scalars or vectors) or - to disable the default rules. METHOD is one of sum, avg, min, max, join. Default is DP:sum,DP4:sum if these fields exist in the input files. Fields with no specified rule will take the value from the first input file. The merged QUAL value is currently set to the maximum. This behaviour is not user controllable at the moment."
                    },
                    {                        
                        "id": "merge",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "snps",
                                    "indels",
                                    "both",
                                    "all",
                                    "none",
                                    "id"
                                ],
                                "name": "merge"
                            }
                        ],
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 10,
                            "valueFrom": "${\n    if(inputs.in_variants.length == 1){\n        return \"\"\n    }else{\n        return \"--merge \" + inputs.merge\n    }\n}"
                        },
                        "label": "Merge option",
                        "doc": "The option controls what types of multiallelic records can be created: none - no new multiallelics, output multiple records instead; snps - allow multiallelic SNP records; indels - allow multiallelic INDEL records; both - both SNP and INDEL records can be multiallelic; all - SNP records can be merged with INDEL records; id - merge by id."
                    },
                    {                        
                        "id": "output_name",
                        "type": "string?",
                        "label": "Output name",
                        "doc": "Output file name."
                    }
                ],
                "outputs": [
                    {
                        "id": "out_variants",
                        "doc": "Merged output file",
                        "label": "Merged output file",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "${  var files_array = [].concat(inputs.in_variants)\n    var in_file = files_array[0]\n    var fname = in_file.basename\n    var fext = in_file.nameext\n    var froot = in_file.nameroot\n    var array_length = files_array.length\n    //if (array_length != 1){\n    \n        if (fext == '.gz') {\n            if (froot.split('.').pop() == 'vcf'){\n                froot = froot.split('.vcf')[0]}\n            else if (froot.split('.').pop() == 'bcf'){\n                froot = froot.split('.bcf')[0]}\n        }\n    \n        if(in_file.metadata.sample_id){\n            froot = in_file.metadata.sample_id}\n    \n        if (inputs.output_name) {\n            var out = inputs.output_name\n            if (inputs.output_type == 'UncompressedVCF') {\n                out += \".merged.vcf\"} \n            else if (inputs.output_type == 'CompressedVCF') {\n                out += \".merged.vcf.gz\"} \n            else if (inputs.output_type == 'UncompressedBCF') {\n                out += \".merged.bcf\"} \n            else if (inputs.output_type == 'CompressedBCF') {\n                out += \".merged.bcf.gz\"} \n            else {\n                out += \".merged.vcf\"}\n        }\n        \n        else if (inputs.output_type == 'UncompressedVCF') {\n            var out = froot + '.merged' + '.vcf'} \n        else if (inputs.output_type == 'CompressedVCF') {\n            var out = froot + '.merged' + '.vcf.gz'} \n        else if (inputs.output_type == 'UncompressedBCF') {\n            var out = froot + '.merged' + '.bcf'} \n        else if (inputs.output_type == 'CompressedBCF') {\n            var out = froot + '.merged' + '.bcf.gz'} \n        else var out = froot + '.merged.vcf'\n    \n        return out\n    //}\n}"
                        },                        
                    }
                ],
                "doc": "**BCFtools Merge and Filter**: Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file, and filter out any monomorphic variants.\n\n**BCFtools** is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, they must be indexed and therefore also compressed. [1]\n\nA list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.\n\n### Common Use Cases\n\n**BCFtools Merge and Filter** should be used when the input files (and expected output files) contain monomorphic variants, e.g. TOPMed datasets.\n\n\n### Changes Introduced by Seven Bridges\n\n* BCFtools works in all cases with gzipped and indexed VCF/BCF files. To be sure BCFtools works in all cases, we added subsequent `bgzip` and `index` commands if a VCF file is provided on input. If VCF.GZ is given on input only indexing will be done. Output type can still be chosen with the `output type` command.\n* If only file file is provided on input, there will be no output.\n\n### Common Issues and Important Notes\n\n* This tool is a third party app developed by the University of Washington Genetic Analysis Center (UW GAC).\n\n* Note that this app is equivalent to running two apps from the BCFtools toolkit sequentially (both of them available in the Public Apps Gallery):  merging of multiple VCF/BCF with **BCFtools Merge** and filtering with **BCFtools view** apps. This app uses piping to perform these tasks more efficiently without writing intermediate data to disk. The command line has the form:  *bcftools merge args | bcftools view --min-ac 1 args*.\n\n* Original BCFtools Merge is faster than BCFtools Merge and Filter since the latter has an additional filtering step included.\n\n* Note that it is the responsibility of the user to ensure that the sample names are unique across all files. If they are not, the program will exit with an error unless the **force samples** (`--force-samples`) option is used. The sample names can also be given explicitly using the print header (`--print-header`) option.\n\n* Note that only records from different files can be merged, but not from the same file. For \"vertical\" merge take a look at BCFtools concat. At least two files must be given to merge.\n\n* BCFtools Merge doesn't work with regular BCFs in some cases.\n\n### Performance Benchmarking\n\nIt took 4 minutes to execute this tool on AWS c4.2xlarge instance using inputs of 12.4 MB and 56 KB. The price is negligible ($0.02).\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1 - BCFtools page](https://samtools.github.io/bcftools/bcftools.html)",
                "label": "Bcftools Merge and Filter_modified",
                "arguments": [
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "${\n    if (inputs.in_variants.length == 1) {\n        return \"python copy_files.py && echo 'only one file specified, nothing to merge' && bcftools view --min-ac 1 \";\n    } else {\n        return \"python copy_files.py && bcftools merge \";\n    }\n\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 101,
                        "valueFrom": "${ \n    //if(inputs.in_variants.length >1){\n    var files_array = [].concat(inputs.in_variants)\n    var in_file = files_array[0]\n    var fname = in_file.basename\n    var fext = in_file.nameext\n    var froot = in_file.nameroot\n    if (fext == '.gz') {\n        if (froot.split('.').pop() == 'vcf'){\n        froot = froot.split('.vcf')[0]}\n        else if (froot.split('.').pop() == 'bcf'){\n        froot = froot.split('.bcf')[0]}\n    }\n\n    if(in_file.metadata.sample_id){\n        var froot = in_file.metadata.sample_id}\n\n    if (inputs.output_name) {\n        var out = inputs.output_name\n        if (inputs.output_type == 'UncompressedVCF') {\n            out += \".merged.vcf\"} \n        else if (inputs.output_type == 'CompressedVCF') {\n            out += \".merged.vcf.gz\"} \n        else if (inputs.output_type == 'UncompressedBCF') {\n            out += \".merged.bcf\"} \n        else if (inputs.output_type == 'CompressedBCF') {\n            out += \".merged.bcf.gz\"} \n        else {\n            out += \".merged.vcf\"}\n    }\n    \n    else if (inputs.output_type == 'UncompressedVCF') {\n        var out = froot + '.merged' + '.vcf'} \n    else if (inputs.output_type == 'CompressedVCF') {\n        var out = froot + '.merged' + '.vcf.gz'} \n    else if (inputs.output_type == 'UncompressedBCF') {\n        var out = froot + '.merged' + '.bcf'} \n    else if (inputs.output_type == 'CompressedBCF') {\n        var out = froot + '.merged' + '.bcf.gz'} \n    else var out = froot + '.merged.vcf'\n\n    return \"--output-file \"+ out\n        \n    //} else{ return \"\" }\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 41,
                        "valueFrom": "${  \n    if(inputs.in_variants.length > 1){\n    \n    var files_array = [].concat(inputs.in_variants);\n    var in_file = files_array[0];\n    var fname = in_file.basename;\n    var fext = in_file.nameext;\n    var froot = in_file.nameroot;\n    \n    if(fext == '.vcf'){\n        return \"./input_files/*.vcf.gz\"\n    }\n    if(fext == '.bcf'){\n        return \"./input_files/*.bcf\"\n    }\n    else if(fext == '.gz'){\n        if (froot.split('.').pop() == 'bcf'){\n        return \"./input_files/*.bcf.gz\"}\n        else if (froot.split('.').pop() == 'vcf'){\n        return \"./input_files/*.vcf.gz\"}\n    }}else{\n        return inputs.in_variants[0].basename;\n    }\n\n\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 51,
                        "valueFrom": "${\n    if(inputs.in_variants.length==1){\n        return \"\"\n    }else{\n    \n     return \"-O u | bcftools view --min-ac 1\"}\n    \n}    "
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": "${\n    if (inputs.mem_per_job) {\n        return inputs.mem_per_job;} \n    else {\n        return 1000;}\n\n}",
                        "coresMin": "${\n    if (inputs.cpu_per_job) {\n        return inputs.cpu_per_job;} \n    else {\n        return 1;}\n}"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerImageId": "21caaa02f72e",
                        "dockerPull": "staphb/bcftools:1.10.2"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "file_paths.txt",
                                "entry": "${  var i;\n    var text = ''\n    var file_array = [].concat(inputs.in_variants);\n    for (i = 0; i < file_array.length; i++) {\n       text += file_array[i].path + '\\n';}\n    \n    return text;\n}",
                                "writable": false
                            },
                            {
                                "entryname": "copy_files.py",
                                "entry": "import os\n\nos.system('mkdir input_files')\n\nwith open('./file_paths.txt') as file:\n    for counter, line in enumerate(file):\n        mv_line = line.strip('\\n')\n        basename = line.split('/')[-1].strip('\\n')\nfile.close()\n\n\nif (counter > 0):\n    with open('./file_paths.txt') as txt:\n        for counter, line in enumerate(txt):\n            mv_line = line.strip('\\n')\n            basename = line.split('/')[-1].strip('\\n')\n            if basename.endswith('.vcf.gz'):\n                nameroot = basename.split('.vcf.gz')\n                new_name = nameroot[0] + \".\" + str(counter) + '.vcf.gz'\n            elif basename.endswith('.vcf'):\n                nameroot = basename.split('.vcf')\n                new_name = nameroot[0] + \".\" + str(counter) + '.vcf'\n            elif basename.endswith('.bcf.gz'):\n                nameroot = basename.split('.bcf.gz')\n                new_name = nameroot[0] + \".\" + str(counter) + '.bcf.gz'\n            elif basename.endswith('.bcf'):\n                nameroot = basename.split('.bcf')\n                new_name = nameroot[0] + \".\" + str(counter) + '.bcf'\n\n            cmd = \"cp \" + mv_line + \" ./input_files/\" + new_name\n            status = os.system(cmd)\n\n    bash_cmd = \"bash /opt/gzip_vcf.sh\"\n    os.system(bash_cmd)\n\nelse:\n    cmd = \"cp \" + mv_line + \" ./\" + basename\n    print(cmd)\n    os.system(cmd)\n            ",
                                "writable": false
                            },
                            {
                                "entryname": "Staging",
                                "entry": "${\n    return inputs.in_variants\n}",
                                "writable": true
                            }
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "var updateMetadata = function(file, key, value) {\n    file['metadata'][key] = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict = {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n        else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict) {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar orderBy = function(files, key, order) {\n    var compareFunction = function(a, b) {\n        if (a['metadata'][key].constructor === Number) {\n            return a['metadata'][key] - b['metadata'][key];\n        } else {\n            var nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n            if (nameA < nameB) {\n                return -1;\n            }\n            if (nameA > nameB) {\n                return 1;\n            }\n            return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n    if (order == undefined || order == \"asc\")\n        return files;\n    else\n        return files.reverse();\n};",
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                
                "stdout": "job.out.log",                
        {
            "id": "plink_filter_makevcffile",
            "in": [
                {
                    "id": "vcf_file",
                    "source": "bcftools_merge_and_filter/out_variants"
                },
                {
                    "id": "Missingness_per_individual",
                    "source": "Missingness_per_individual"
                },
                {
                    "id": "Missingness_per_marker",
                    "source": "Missingness_per_marker"
                },
                {
                    "id": "Hardy_Weinberg_equilibrium",
                    "source": "Hardy_Weinberg_equilibrium"
                },
                {
                    "id": "Mendel_error_rates",
                    "source": [
                        "Mendel_error_rates"
                    ]
                },
                {
                    "id": "Biallelic",
                    "source": "Biallelic"
                },
                {
                    "id": "SNPs_only",
                    "source": "SNPs_only"
                },
                {
                    "id": "MAF_Threshold",
                    "source": "MAF_Threshold"
                },
                {
                    "id": "QUAL",
                    "source": "QUAL"
                }
            ],
            "out": [
                {
                    "id": "plink_output_vcf"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",                
                "baseCommand": [
                    "plink2"
                ],
                "inputs": [
                    {
                        "format": [
                            "VCF",
                            "GZIP"
                        ],
                        "id": "vcf_file",
                        "type": "File",
                        "inputBinding": {
                            "prefix": "--vcf",
                            "shellQuote": true,
                            "position": 0
                        },
                        "doc": "Input VCF.gz file."
                    },
                    {
                        "id": "Missingness_per_individual",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--mind",
                            "shellQuote": false,
                            "position": 3
                        },
                        "doc": "--mind\n\nexclude individuals with too much missing genotype data."
                    },
                    {
                        "id": "Missingness_per_marker",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--geno",
                            "shellQuote": false,
                            "position": 3
                        },
                        "doc": "--geno\n\nexclude SNPs on the basis of missing genotype rate, with the --geno option: the default is to include all SNPS (i.e. --geno 1). To include only SNPs with a 90% genotyping rate (10% missing) use --geno 0.1"
                    },
                    {
                        "id": "Hardy_Weinberg_equilibrium",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--hwe",
                            "shellQuote": false,
                            "position": 3
                        },
                        "doc": "--hwe\n\nTo exclude markers that failure the Hardy-Weinberg test at a specified significance threshold"
                    },
                    {
                        "id": "Mendel_error_rates",
                        "type": "float[]?",
                        "inputBinding": {
                            "prefix": "--me",
                            "shellQuote": false,
                            "position": 3
                        },
                        "doc": "example:\nplink --file mydata --me 0.05 0.1\n\nFor family-based data only, to exclude individuals and/or markers on the basis on Mendel error rate.\n\n\n-  the first parameter determines that families with more than 5% Mendel errors (considering all SNPs) will be discarded.\n-  the second parameter indicates that SNPs with more than 10% Mendel error rate will be excluded (i.e. based on the number of trios);"
                    },
                    {                        
                        "id": "Biallelic",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--max-alleles",
                            "shellQuote": false,
                            "position": 2,
                            "valueFrom": "${\r\n    if (inputs.Biallelic) {\r\n      return 2;\r\n    } else {\r\n      return null;\r\n    }\r\n  }"
                        },
                        "doc": "--max-alleles 2 filters out the multiallelic variants"
                    },
                    {
                        "id": "SNPs_only",
                        "type": "boolean?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 2,
                            "valueFrom": "${\r\n    if (inputs.SNPs_only) {\r\n      return \"--snps-only\";\r\n    } else {\r\n      return null;\r\n    }\r\n  }"
                        },
                        "doc": "Set TRUE if you want to include this filter:\n\n--snps-only excludes all variants with one or more multi-character allele codes. With 'just-acgt', variants with single-character allele codes outside of {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', <missing code>} are also excluded."
                    },
                    {
                        "id": "MAF_Threshold",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--maf",
                            "shellQuote": false,
                            "position": 3
                        },
                        "doc": "--maf\n\nexclude SNPs on the basis of MAF (minor allele frequency)"
                    },
                    {
                        "id": "Pass_Filter",
                        "type": "boolean?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 2,
                            "valueFrom": "${\r\n    if (inputs.PASS_FILTER) {\r\n      return \"--var-filter\";\r\n    } else {\r\n      return null;\r\n    }\r\n  }\r\n  \r\n"
                        },
                        "doc": "To skip variants which failed one or more filters tracked by the FILTER field, use --var-filter. This can be combined with one or more (space-delimited) filter names to ignore."
                    },
                    {
                        "id": "QUAL",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--var-min-qual",
                            "shellQuote": false,
                            "position": 2
                        },
                        "doc": "--var-min-qual causes all variants with QUAL value smaller than the given number, or with no QUAL value at all, to be skipped."
                    },
                    {
                        "id": "Sample_ID",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--keep",
                            "shellQuote": false,
                            "position": 1
                        }
                    },
                    {
                        "id": "Pheno",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--pheno",
                            "shellQuote": false,
                            "position": 1,
                            "valueFrom": "$(inputs.Pheno.path) --1"
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "plink_output_vcf",
                        "doc": "Output PLINK binary BED file.",
                        "type": "File[]",
                        "outputBinding": {
                            "glob": "*_plink_output.vcf.gz"
                        }
                    }
                ],
                "doc": "### PLINK Tool Description\n\nThe PLINK tool is a widely recognized and used software in the genomic research community, developed by Shaun Purcell at the Broad Institute and Harvard Medical School. It provides a comprehensive range of functionalities for analyzing and manipulating genomic data. Within this workflow, a custom version of PLINK2 is used to simplify the toolset, ensuring that only the necessary functions for this specific application are included. This approach makes the pipeline both user-friendly and efficient.\n\n&nbsp;\n\n#### Inputs:\n- **Merged WGS data (VCF.gz):** The filtered and merged genomic data in VCF.gz format, prepared in the previous step using BCFtools.\n\n&nbsp;\n\n#### Outputs:\n- **Filtered WGS data (VCF.gz):** The refined dataset after applying various filters and quality controls, ready for further processing and analysis.\n\n&nbsp;\n\n#### Function:\nThe custom PLINK2 app performs several key functions to refine the merged WGS data:\n\n1. **ID Matching:**\n   - The tool matches IDs with phenotypes using the `--keep` option, ensuring that the data corresponds accurately to specific study participants.\n\n2. **Filtering:**\n   - **Minor Allele Frequency (MAF):** Filters variants based on their minor allele frequency to retain variants that are common enough to be statistically significant.\n   - **Hardy-Weinberg Equilibrium (HWE):** Applies filters to ensure that the genotype frequencies are consistent with the expected Hardy-Weinberg proportions.\n   - **Genotype Missingness (Geno):** Removes variants with a high proportion of missing genotype data.\n   - **Sample Missingness (Mind):** Excludes samples with excessive missing data.\n   - **Quality Metrics (Qual):** Applies quality control measures to ensure that the retained variants meet the necessary standards for reliable analysis.\n   - **Biallelic Variants:** If set to true, adds `--max-alleles 2` to include only biallelic variants.\n   - **SNPs Only:** If set to true, adds `--snps-only` to include only SNPs.\n   - **Mendelian Error Rates:** Applies `--me` to filter based on Mendelian error rates.\n   - **Pass Filter:** If set to true, adds `--var-filter` to include only variants that pass all filters.\n   - **QUAL:** An integer variable that applies `--var-min-qual` to include variants with a minimum quality score.\n\n3. **Optional Phenotype Application:**\n   - If a phenotype file is provided, it applies `--pheno pheno_file.txt --1`. However, this step is not necessary for the current analysis.\n\n&nbsp;\n\n#### Benefits:\n- **Efficiency:** The custom PLINK2 app reduces the size of the dataset and improves its quality, ensuring that only high-quality data is passed on for further processing.\n- **Ease of Use:** By simplifying the toolset, the custom PLINK2 app makes it easier for users to perform essential genomic data refinement without the need for extensive command-line knowledge.\n\n&nbsp;\n\nBy refining and filtering the merged WGS data, the PLINK2 tool prepares a high-quality dataset that is essential for downstream analyses in the genomic analysis pipeline. This step ensures the integrity and reliability of the data, which is crucial for accurate and meaningful genetic association studies.",
                "label": "PLINK",
                "arguments": [
                    {
                        "prefix": "--out",
                        "shellQuote": false,
                        "position": 5,
                        "valueFrom": "$(inputs.vcf_file.nameroot)_plink_output"
                    },
                    {
                        "prefix": "--export",
                        "shellQuote": false,
                        "position": 4,
                        "valueFrom": "vcf bgz"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 4,
                        "valueFrom": "--missing"
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "biocontainer/plink2"
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    }
                ],
                "stdout": "plink.log",
                "stderr": "plink.err",                
            },
            "label": "PLINK",
            "scatter": [
                "vcf_file"
            ],            
        },
        {
            "id": "datafiltering2",
            "in": [
                {
                    "id": "FilteredVCF",
                    "source": "plink_filter_makevcffile/plink_output_vcf",
                    "pickValue": "first_non_null"
                }
            ],
            "out": [
                {
                    "id": "compressed_SNP_info"
                },
                {
                    "id": "annotation_csv"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",                
                "baseCommand": [
                    "python",
                    "myscript.py"
                ],
                "inputs": [
                    {
                        "id": "FilteredVCF",
                        "type": "File",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 0
                        },                        
                    }
                ],
                "outputs": [
                    {
                        "id": "compressed_SNP_info",
                        "label": "compressed_SNP_info",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": [
                                "*_subIDs.txt",
                                "*_varIDs.txt",
                                "*_variants_value.npy"
                            ]
                        }
                    },
                    {
                        "id": "annotation_csv",
                        "label": "Annotations CSV",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "*_annotations.csv"
                        }
                    }
                ],
                "doc": "### SNP_prep Tool Description\n\nThe SNP_prep tool processes the filtered whole genome sequencing (WGS) data, transforming it into a clean, compressed format suitable for downstream analysis using IBI or Fisher methods. This step is crucial for ensuring that the genomic data meets the stringent requirements of advanced analytical tools, enhancing its utility and accuracy for genetic studies.\n\n&nbsp;\n\n#### Inputs:\n- **Filtered WGS data (VCF.gz):** The refined dataset in VCF.gz format, prepared in the previous steps using BCFtools and PLINK2.\n\n&nbsp;\n\n#### Outputs:\n- **Compressed SNP matrix (`variants_values.npy`):** A binary format file containing the preprocessed SNP data, allowing efficient storage and rapid access during analysis.\n- **Subject IDs (`subIDs.txt`):** A text file listing the subject IDs corresponding to the SNP data matrix.\n- **Variant IDs (`varIDs.txt`):** A text file listing the variant IDs corresponding to the SNP data matrix.\n- **Annotations (`annotations.csv`):** A CSV file containing additional annotations for the SNPs, providing contextual information for downstream analyses.\n\n&nbsp;\n\n#### Function:\nThe SNP_prep tool performs several key functions to process the filtered WGS data:\n\n1. **Separation of SNP Matrix and Annotations:**\n   - The tool separates the SNP matrix from the annotations, ensuring that the dataset is structured and organized for analysis.\n\n2. **Dominant Encoding:**\n   - It applies dominant encoding to the genotype data, where each allele is encoded as 0 or 1. In this encoding scheme, 1 represents the presence of the dominant allele, and 0 represents its absence. This simplification reduces data size and complexity.\n\n3. **Data Compression:**\n   - The tool compresses the SNP data into a `.npy` format. This binary format, used by the Python programming language, allows for efficient storage and rapid read/write operations, making it ideal for handling large genomic datasets.\n\n&nbsp;\n\n#### Benefits:\n- **Efficiency:** By transforming the data into a binary format and applying dominant encoding, the SNP_prep tool ensures that the dataset is both compact and easily accessible for analysis.\n- **Integration:** The tool produces outputs that are directly compatible with IBI and Fisher methods, facilitating seamless integration within the genomic analysis pipeline.\n- **Accuracy:** The meticulous processing and organization of the data enhance its accuracy and utility for genetic studies, ensuring that downstream analyses are based on high-quality input.\n\n&nbsp;\n\nBy processing the filtered WGS data into standardized, analysis-ready formats, the SNP_prep tool serves as a vital link in the genomic analysis pipeline. It ensures that the resulting datasets are primed for high-precision analysis, enabling researchers to identify significant genetic markers associated with various diseases.",
                "label": "SNP_prep",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": "$(inputs.FilteredVCF.size *250 * 1e-6)"
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
                                "entry": "import gzip\nimport os\nimport time\nimport pandas as pd\nimport numpy as np\nfrom datetime import datetime\nimport scipy.stats as stats\nfrom tqdm import tqdm\nimport logging\nimport tooloptions\n\nimport gc\n\n## Find the number of rows to skip\n\n# Define the header keywords to search for\nheader_keywords = [\"#CHROM\", \"POS\", \"ID\", \"REF\", \"ALT\"]\n\n# Initialize variables\nheader_row_number = None\nline_number = 0\n\n# Open the file with 'latin-1' encoding\nwith gzip.open(tooloptions.gzpath, 'rb') as file:\n    # Loop through the lines in the file\n    for line_bytes in file:\n        line_number += 1\n        # Decode the line from bytes to a string using 'latin-1' encoding\n        line = line_bytes.decode('latin-1')\n        # Check if the line contains all the header keywords\n        if all(keyword in line for keyword in header_keywords):\n            header_row_number = line_number\n            break  # Stop searching once the header is found\n        \n        \n        \n\n\nwith gzip.open(tooloptions.gzpath,'r') as gz_file:        \n      df = pd.read_csv(gz_file, header=0, delimiter=\"\\t\", index_col=0, skiprows=header_row_number - 1)\n        \n        \n\n##################################################################################\nprint(\"raw df:\", df)\n\n\n#################### double check on PASS filter\n# non_PASS = df[~(df[\"FILTER\"]==\"PASS\")]\n# if non_PASS.shape[0] == 0:\n#     print(\"all variants have PASS quality\")\n# else:\n#     df = df[(df[\"FILTER\"]==\"PASS\")]\n    \n\n################## assigning unique identifier chr:index\ndf = df.reset_index()\ndf['index'] = list(range(df.shape[0]))\n# Adding a new column by combining 'chr' and 'index' columns\ndf['index'] = df.apply(lambda row: f\"{row['#CHROM']}:{row['index']}\", axis=1)\ndf = df.set_index('index')\n\n\nvcf_col = df[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']]\nvcf_SNPs = df[df.columns[~df.columns.isin(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])]]\n\n\ndel df\ngc.collect()\n\nprint(\"vcf_SNPs shape\", vcf_SNPs.shape)\n\n###########################################################\n################################# Output file name\n\nfile_name = tooloptions.input_file_basename\nfile_name = file_name.split('.')[0]\n\n\nChr_num = vcf_col[\"#CHROM\"].mode()[0]\nchromosome = f\"chr{Chr_num}\"\n\n\nfile_name = f\"{file_name}_{chromosome}\"\n\n############## annotation columns\n###############################################################\nannot_file_name = '{}_annotations.csv'.format(file_name)\nvcf_col.to_csv(annot_file_name)\n\n\n##################### handling missing values\n##################################################################\n##### TopMEd freeze 9 cohorts do not have any missing values ######\n# # Record the start time\n# start_time = time.time()\n\n# def fill_mode(row):\n#     # calculate the mode\n#     mode_value=row.mode()[0] if not row.mode().empty else 0\n#     return row.fillna(mode_value)\n\n\n# # unique_value = np.unique(vcf_SNPs.values)\n# null_counts = vcf_SNPs.isna().sum()\n\n\n# print(\"The percent of each SNP with less than 50 missing values: \", \n#       100*np.round(null_counts[null_counts<50].shape[0]/null_counts.shape[0],2))\n\n\n# if null_counts[null_counts>0].shape[0] > 0:\n#     print(null_counts[null_counts>0].shape[0], \" SNPs have missing values.\")\n#     print(\"Filling missing values withthe mode\")\n#     # Fill missing values with the mode\n#     df_filled = df.apply(lambda row:fill_mode(row),axis=1)\n    \n# # Record the end time\n# end_time = time.time()\n# # Calculate the difference\n# elapsed_time = end_time - start_time\n\n# print(f\"The code for handling missing values ran for {elapsed_time} seconds.\")\n\n\n\n\n################### dominant encoding \n#############################################################\ndf_ready = vcf_SNPs.copy()\n\ndf_ready_replaced = df_ready.replace(to_replace ='0/0', value =0)\ndf_ready_replaced = df_ready_replaced.replace(to_replace ='0/1', value =1)\ndf_ready_replaced = df_ready_replaced.replace(to_replace ='1/0', value =1)\ndf_ready_replaced = df_ready_replaced.replace(to_replace ='1/1', value =1)\n\n\n# vcf_dominant012.to_csv(\"SNPs_matrix_dominant_encoded.csv\") ##final output data\n\n\n################### compressing files \n#############################################################\ndef compress_file(df):\n    st = time.time()\n    subIDs = list(df.columns)\n    data = np.array(df, dtype=np.int8)\n    varIDs = df.index.astype(str)\n\n    # Use the output file names from tooloptions\n    subIDs_file_name = '{}_subIDs.txt'.format(file_name.split('.')[0])\n    varIDs_file_name = '{}_varIDs.txt'.format(file_name.split('.')[0])\n    npy_file_name = '{}_variants_value.npy'.format(file_name.split('.')[0])\n\n    with open(subIDs_file_name , 'w') as file:\n        file.write('\\n'.join(subIDs))\n    with open(varIDs_file_name, 'w') as file:\n        file.write('\\n'.join(varIDs))\n    np.save(npy_file_name, data)\n    print(\"compressing data elapsed time:{}\".format(time.time() - st))\n\n\n\n\ncompress_file(df_ready_replaced)\n",
                                "writable": false
                            },
                            {
                                "entryname": "tooloptions.py",
                                "entry": "gzpath = \"$(inputs.FilteredVCF.path)\"\n\ninput_file_basename = \"$(inputs.FilteredVCF.basename)\"\n",
                                "writable": false
                            }
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    }
                ],
                
                "stdout": "standard.out",
                
            },
            "label": "SNP_prep",
            "scatter": [
                "FilteredVCF"
            ],
            
        }
    ],
    "requirements": [
        {
            "class": "ScatterFeatureRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "StepInputExpressionRequirement"
        }
    ],
    
}