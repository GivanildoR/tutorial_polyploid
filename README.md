# Pipeline to working with polyploid genomic Flex-Seq data

A personal pipeline to auxiliary on the manipulation of polyploid genomic data 🥔🫐 🍠  
It was created to auxiliar some processes in the lab, I would highlight that a bunch of scripts in R, and Python languages and command lines found here can be found in different repositories/journals.   
Any comments or suggestions are welcome.  

## Organizing the genomic data

Used data  
```bash
  Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf
```
Obtaining and verifying the COD genotype names in the VCF file:  
```bash
awk '/^#CHROM/ {print}' Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf
```  
Let's presume that you have a phenotypic table and you want to consider just the samples with phenotypic information. Then, the genotypes with no information should be removed from the VCF file, to do this, let's consider:    
We should create a .txt file containing the new codes, in this case, "Customer_Code". The first column is the "old name" and the second column is the "new name",  **without header**:  
The first 10 code names of the `genotype_names_mapping.txt` file: 

<img width="156" alt="Imagem2" src="https://github.com/GivanildoR/tutorial_polyploid/assets/167666189/cbc07953-fa9d-47c3-bc1f-303bfc0e7fa5">

Based on the first column of `genotype_names_mapping.txt` let's filter the VCF file:

Step #1

```bash
cut -f1 genotype_names_mapping.txt > lista_ids.txt
```

Step #2 using the `bcftools` to filter:

```bash
bcftools view -S lista_ids.txt CopyOfFinal_DP10_Corrected_UFL_137106_RAW_SNPs.vcf -o
 1_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf
```
Put the renomear_amostras.py in the same directory as your genomic data and run
import sys

```python
# Argumentos do script: [1] Arquivo de mapeamento, [2] Arquivo VCF original
map_file_path = sys.argv[1]
vcf_file_path = sys.argv[2]

# Ler o arquivo de mapeamento para construir o dicionário de renomeação
rename_map = {}
with open(map_file_path, 'r') as map_file:
    for line in map_file:
        old_name, new_name = line.strip().split('\t')
        rename_map[old_name] = new_name

# Gerar novo cabeçalho VCF
with open(vcf_file_path, 'r') as vcf:
    for line in vcf:
        if line.startswith('#CHROM'):
            columns = line.strip().split('\t')
            new_columns = [rename_map.get(col, col) for col in columns]
            print('\t'.join(new_columns))
        else:
            print(line.strip())
```

Running the python script:

```python
#python renomear_amostras.py genotype_names_mapping.txt 1_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf > 2_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf
```

Just to check the genotype names in the final VCF

```bash
bcftools query -l 2_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf
```

## Filtering by chromosomes  
### Sometimes the VCF file has specific parts from another reference genome just to verify any particular aspect in your panel.   
### Sometimes you want to remove some scaffolds, anyways, you can verify the chromosomes names using:  
### if it isn't your case, skip this step.  

If you don't have the `.gz` format  
```bash
bgzip -c 2_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf > 2_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf.gz
```

Checking:  
```bash
awk '!/^#/ {print $1}' 2_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf | sort -u
```

Filtering by chromosomes` names, in my case I'm interested just on the chr01 to chr12  
```bash
#the vcf file should be indexed
bcftools index 2_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf.gz

#Now we should keep just the interested chromosomes
bcftools view -r chr01,chr02,chr03,chr04,chr05,chr06,chr07,chr08,chr09,chr10,chr11,chr12 2_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf.gz -o 3_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf

#verfying the VCF file obtained
awk '!/^#/ {print $1}' 3_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf | sort -u
```

## STRUCTURE analysis
### Filtering SNPs by LD
To realize the Structure analysis we should filter using the LD (the Structure software doesn't work with marker in high LD.  
To prune the SNPs by the LD, we can use the following command:  
Here was considered LD=0.2 in windows with 200 bp, because there are clusters with ~200bp, remember it is not a whole sequence GBS.  
```bash
module load bcftools
bcftools +prune -m 0.2 -w 200 Vcf_Sampled.vcf -Ob -o myvariants.vcf
```
Using the 




