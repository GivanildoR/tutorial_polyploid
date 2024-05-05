# Working with polyploid genomic Flex-Seq data - an introductory pipeline  

It is a personal pipeline to auxiliary on the manipulation of polyploid genomic data 🥔  
It was created to auxiliar and to evaluate the processes. A bunch of scripts in R, and Python languages and command lines found here can be found in different repositories/journals. Here, I could save some considerable important information for me 💻  

Firstly, [here](https://eriqande.github.io/eca-bioinf-handbook/handle-vcf.html), there is a tutorial for noobs on the command line like me.   
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
Put the [rename_samples.py](https://github.com/GivanildoR/tutorial_polyploid/blob/main/rename_samples.py) in the same directory as your genomic data and run

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
To realize the Structure analysis we should filter using the LD (the Structure software doesn't work with markers in high LD).  
To prune the SNPs by the LD, we can use the following command:  
Here was considered LD=0.2 in windows with 200 bp, because there are clusters with ~200bp, remember it is not a whole genome sequence GBS.  
```bash
module load bcftools
bcftools +prune -m 0.2 -w 200 3_Final_DP10_Corrected_UFL_137106_RAW_SNPs.vcf -Ob -o myvariants.vcf
```
We can verify the content of `myvariants.vcf` to do this, firstly we should compress it to `.gz` format     
```bash
# 1 step
bgzip -c myvariants.vcf > myvariants.vcf.gz

# 2 step
bcftools index myvariants.vcf.gz

# 3 how many variants are in this file?  
bcftools stats myvariants.vcf | less   
```

Now, we can select using a random way to reduce the SNP number   
Go to R software:  

```R
setwd("/blue/mresende/share/Givanildo/Tutorial")
library(GWASpoly) #Endelman package

VCF2dosage(
  VCF.file = "myvariants.vcf",
  dosage.file = "myvariants_allele.csv",
  geno.code = "GT",
  ploidy = 4,
  samples = NULL,
  min.DP = 1,
  max.missing = 0.6,
  min.minor = 5
)

```
Now use the select_random_snp_per_interval to reduce the quantity of SNP  
[select_random_snp_per_interval.R](https://github.com/GivanildoR/tutorial_polyploid/blob/main/select_random_snp_per_interval.R)  

### Converting to Structure format

The matrix should be in dimension genotypes in the row names and markers in the column names  
[numeric2structure.R](https://github.com/GivanildoR/tutorial_polyploid/blob/main/numeric2structure.R) 

## Running the STRUCTURE analysis  
### mainparams  

[mainparams](https://github.com/GivanildoR/tutorial_polyploid/blob/main/mainparams.R)

### extraparams  
[extraparams](https://github.com/GivanildoR/tutorial_polyploid/blob/main/extraparams.R)

### submitting the Structure job on the HPG  



