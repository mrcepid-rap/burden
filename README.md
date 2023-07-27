# RunAssociationTesting – Burden

This is the source code for a module that is implemented in mrcepid-runassociationtesting.
For more information about how to run or modify it, see https://github.com/mrcepid-rap/mrcepid-runassociationtesting.

### Table of Contents

- [Introduction](#introduction)
  * [Changelog](#changelog)
  * [Background](#background)
    + [Burden Tools Implemented](#burden-tools-implemented)
      - [1. [BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)](#1-bolt-lmmhttpsalkesgroupbroadinstituteorgbolt-lmmbolt-lmm_manualhtml)
      - [2. [SAIGE-GENE+](https://github.com/saigegit/SAIGE)](#2-saige-genehttpsgithubcomsaigegitsaige)
      - [3. [STAAR](https://github.com/xihaoli/STAAR)](#3-staarhttpsgithubcomxihaolistaar)
      - [4. [REGENIE](https://rgcgithub.github.io/regenie/)](#4-regeniehttpsrgcgithubgithubioregenie)
      - [5. Generalised Linear Models (GLMs)](#5-generalised-linear-models-glms)
- [Methodology](#methodology)
    + [BOLT](#bolt)
    + [SAIGE-GENE+](#saige-gene)
    + [STAAR](#staar)
    + [REGENIE](#regenie)
    + [GLMs](#glms)
- [Running on DNA Nexus](#running-on-dna-nexus)
  * [Inputs](#inputs)
    + [Association Tarballs](#association-tarballs)
  * [Outputs](#outputs)
    + [Per-gene output](#per-gene-output)
    + [Per-marker output](#per-marker-output)

## Introduction

This module performs rare variant burden testing implemented with one of five different methods:

* [BOLT](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)
* [SAIGE-GENE+](https://github.com/saigegit/SAIGE)
* [STAAR](https://github.com/xihaoli/STAAR)
* [REGENIE](https://rgcgithub.github.io/regenie/)
* GLMs – vanilla linear/logistic models implemented with python's [statsmodels module](https://www.statsmodels.org/stable/index.html)

This README makes use of DNANexus file and project naming conventions. Where applicable, an object available on the DNANexus
platform has a hash ID like:

* file – `file-1234567890ABCDEFGHIJKLMN`
* project – `project-1234567890ABCDEFGHIJKLMN`

Information about files and projects can be queried using the `dx describe` tool native to the DNANexus SDK:

```commandline
dx describe file-1234567890ABCDEFGHIJKLMN
```

### Changelog

* v1.1.4
  * Implemented the CommandExecutor class

* v1.1.3
  * Implemented _check_opts method for interface compatibility
  * Minor README changes to reflect changes in linear model code refactoring 

* v1.1.2
  * Internally changed return types from str to Path to enable compatibility with updated utilities
  * Minor method parameters change to covariate processing when running burden tests for compatibility
    * No changes were made to methodology

* v1.1.1
  * Minor changes to BOLT commands to facilitate BOLT upgrade to v2.4.1

* v1.1.0
  * Start of major effort to annotate code better and provide unit tests for functionality

* v1.0.0
  * Initial release. The functionality of this module was formally embedded in RunAssociationTesting

### Background

This module performs rare-variant burden tests to 

#### Burden Tools Implemented

Aggregation can be done in multiple different ways. In this applet, we have implemented five different approaches. 
Each of the five implemented approaches has advantages and disadvantages. I have summarised some primary examples in 
the sections below for each tool.

##### 1. [BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)

Original Publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4342297/ 

BOLT was originally designed to perform association testing for common variants. BOLT does not, by default, allow for 
rare variant burden testing. Nonetheless, here we have implemented an aggregation approach whereby we create "dummy" variables
to trick BOLT into performing rare variant burden tests. For more information, please see the [section in the methodology](#bolt) 
that covers BOLT in more detail.

Advantages:

* Relatively Fast
* Handles cryptic relatedness

Disadvantages:

* Cannot handle quantification of number of qualifying variants – only has a binary "yes/no" that an individual has a qualifying variant

##### 2. [SAIGE-GENE+](https://github.com/saigegit/SAIGE)

Original Publication: https://www.medrxiv.org/content/10.1101/2021.07.12.21260400v2

Like BOLT, SAIGE was originally designed to perform association testing for common variants. Unlike BOLT, the authors of
SAIGE have adapted it for running rare variant association tests in the form of SAIGE-GENE+. Crucially, SAIGE-GENE+ implements 
the [SKAT-O](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3415556/) methodology to handle known issues with p. value
approximation with extreme case-control imbalance of binary traits which reduces false-positive rates when performing
association testing. SAIGE-GENE+ performs two tests: 1) it first performs a test where the total number of qualifying
variants (e.g. PTV, missense, etc.) per individual per gene are summed and 2) a per-variant test to look for an
association at all individual qualifying variants. For more information, please see the [section in the methodology](#saige-gene) 
that covers SAIGE-GENE+ in more detail.

Advantages:

* Handles cryptic relatedness
* Adjusts for variance of rare variants
* Natively handles rare variant aggregation

Disadvantages

* Relatively slow
* High memory usage

##### 3. [STAAR](https://github.com/xihaoli/STAAR)

Original Publication: https://www.nature.com/articles/s41588-020-0676-4

STAAR is a rare variant burden testing method that we have implemented to provide an alternate approach to the two methods
listed above. In brief, it uses a similar approach to SAIGE-GENE+, whereby qualifying variants are aggregated per 
individual per gene. The primary difference in this approach is that it uses a pre-built set of familial relationships
to handle cryptic relatedness between individuals being tested rather than calculating this relationship during execution.
STAAR also ostensibly contains additional functionality to assign different weights to variants based on a variety of 
different annotations (i.e. methylation, open chromatin, etc.). We have not implemented this as we are 1) primarily dealing
with coding sequence and 2) aggregating variants ourselves. We may implement this model as time or necessity allows. For
more information, please see the [section in the methodology](#staar) that covers STAAR in more detail.

Advantages:

* Very fast

Disadvantages:

* Implemented in R as a package rather than as a command-line application (requires creation of wrapper script)
* Does not natively handle cryptic relatedness

##### 4. [REGENIE](https://rgcgithub.github.io/regenie/)

Original Publication: https://www.nature.com/articles/s41588-021-00870-7

REGENIE is very similar to, and has a similar history to, BOLT and SAIGE. Originally developed for GWAS, it now has 
methodology for performing rare variant burden tests.

Advantages:

* Controls for relatedness
* Runs gene-burden tests by default
* Implements several different statistical models (e.g. SKAT, SKAT-O, ACAT, ACAT-O)

Disadvantages:

* Memory, time, and CPU intensive

##### 5. Generalised Linear Models (GLMs)

Original Publication: N/A

GLMs are among the most basic approach that can be implemented for rare variant burden testing. In the case of GLMs, we 
have implemented a simple aggregation approach whereby number of qualifying variants per individual per gene are counted
and tested for association with a given phenotype using either a linear or logistic model. The current approach has been
implemented using the [statsmodels module](https://www.statsmodels.org/stable/index.html) for python3. For more information, 
please see the [section in the methodology](#glms) that covers GLMs in more detail.

Advantages:

* Relatively fast
* The least complex model (less of a black box)

Disadvantages:

* Does not control for case-control imbalance when calculating p. value
* Does not control for cryptic relatedness

## Methodology

In theory, each tool could have been implemented using individual applets. This was decided against in order to simplify
the inputs provided to each tool. By including covariate and sample filtering within the applet as the first step, we can ensure
that covariates, samples, and phenotypes are run **identically** for all four tools/methods.

Here we have documented in more detail covariate processing and then the implementation of each of the four tools documented
in the [background](#background) section of this README. Note that when we use the term "<file_prefix>", that refers to the
prefix that was provided to mrcepid-collapsevariants and mrcepid-mergecollapsevariants.

#### BOLT

BOLT roughly proceeds in two steps. First, BOLT computes a genetic relatedness matrix (GRM) from genetic data curated during
the previous step of this applet. Next, it tests for the association a variant provided in a bgen file for an association
with the phenotype of interest. Normally, this bgen file is a collection of phased SNP data from genotying chips. Here,
we have created a "dummy" bgen file that instead encodes presence/absence of a qualifying variant per-individual. This dummy
variable is then used by BOLT to test for an association with a given phenotype. We test for both per-gene and per-marker 
(i.e. SNV/InDels) association with our phenotype of interest.

##### Inputs

A bgen file of genes to test for association with our phenotype of interest. This file is initially formated as a plink 
.ped file like the following:

```text
1000000 1000000 0 0 0 -9 A A A A
1000001 1000001 0 0 0 -9 A C A A
1000002 1000002 0 0 0 -9 A C A A
1000003 1000003 0 0 0 -9 A A A C
1000004 1000004 0 0 0 -9 A A A A
```

The first 6 columns are standard [.ped format](https://www.cog-genomics.org/plink/1.9/formats#ped). The next four columns
are the representation of the dummy variable. Each pair of columns (i.e. 7 & 8, 9 & 10) represent one gene that we are testing. 
Individuals WITH a qualifying variant (e.g. a PTV, missense, etc) are labelled as heterozygous (A C); individuals without a 
qualifying variant are labelled as homozygous (A A). To be clear, this means that individuals 1000001 and 1000002 have a 
qualifying variant in Gene1, and 1000003 has a qualifying variant in Gene2. An accompanying [.map file](https://www.cog-genomics.org/plink/1.9/formats#map)
is provided with this of the format:

```text
chr1 ENST000000001 0 1000000
chr1 ENST000000002 0 2000000
```

This file is then converted to .bed using plink2 and then .bgen 1.2 8-bit, ref-last format using plink2:

```commandline
plink --make-bed --file <file_prefix>.BOLT --out <file_prefix>.BOLT
plink2 --export bgen-1.2 'bits='8 --bfile <file_prefix>.BOLT --out <file_prefix>.BOLT"
```

The above steps are done per-chromosome and provided via the `--bgenSamplesFileList` argument as described below.

To perform per-marker tests, we simply convert the VCF file used for [SAIGE-GENE+](#saige-gene) into bgen format using plink2, 
and run exactly the same command as described below. Inputs and outputs are essentially identical. 

##### Command Line Example

```commandline
bolt --bfile=UKBB_200K_Autosomes_QCd_WBA /                                       # Filtered genetic data
        --exclude=UKBB_200K_Autosomes_QCd.low_MAC.snplist /                      # List of SNP IDs 
        --phenoFile=phenotypes_covariates.formatted.txt /                        # formated phenotype + covariate file generated during step 1
        --phenoCol=<pheno_name> /                                                # phenotype name extracted from the provided phenotype file
        --covarFile=phenotypes_covariates.formatted.txt /                        # BOLT requires you to provide this twice...
        --covarCol=batch /                                                       # batch is field 22000, we provide here as a categorical variable (covarCol)
        --covarCol=sex /                                                         # sex as a categorical variable
        --covarCol=wes_batch /                                                   # Batch a given individual's WES is derived from (50k / 200k / 450k) 
        --qCovarCol=age /                                                        # age as a continuous variable
        --qCovarCol=PC{1:10} /                                                   # First 10 principal components as a continuous variable
        --covarMaxLevels=110 /                                                   # Since there are 110 batches of arrays (batch), we have to tell BOLT to allow for more categorical bins
        --LDscoresFile=BOLT-LMM_v2.3.5/tables/LDSCORE.1000G_EUR.tab.gz /         # LD scores pre-computed by BOLT
        --geneticMapFile=BOLT-LMM_v2.3.5/tables/genetic_map_hg19_withX.txt.gz /  # Genetic map pre-computed by BOLT
        --lmmInfOnly /                                                           # Tell bolt to compute only infinitesimal mixed model association statistics
        --numThreads=32 /                                                        # Number of threads
        --statsFile=<file_prefix>.tarball_prefix.stats.gz /                      # This sets the output file of basic stats
        --verboseStats /                                                         # Give verbose stats
        --bgenSampleFileList=bolt_input.bgen /                                   # This is a list (one file for each chromosome) of dummy bgen files created above
        --statsFileBgenSnps=<file_prefix>.bgen.stats.gz                          # output statistics for testing dummy variables
```

**Note:** BOLT does not differentiate between binary and continuous traits.

#### SAIGE-GENE+

##### Inputs

SAIGE-GENE+ requires 2 input files that are created during mrcepid-mergecollapsevariants:

1. A VCF format file of rare variants we want to test.

This file was generated during mrcepid-collapsevariants for each VCF file and merged across VCFs during mrcepid-collapsevariants
into a plink .ped format file. This .ped file was then converted into a VCF that contained ONLY variants we want to test using
a command like:

```commandline
plink2 --pfile <file_prefix> --export vcf --out <file_prefix>.SAIGE
```

2. A "groupFile" that tells SAIGE which variants to combine for a given gene.

Each gene we are testing is represented by one tab-delimited row, where the first column is the gene name and subsequent columns
are variants we want to test for that gene. Variant ID is represented as the VCF fields CHROM:POS_REF/ALT and **NOT** by the ID
column.

```text
ENST000000001   1:1000000_A/T   1:1000010_T/G   1:1000020_G/A
ENST000000001   foo foo foo
ENST000000002   1:2000000_G/C   1:2000030_A/C   1:2000050_ATC/A 1:2000000_G/GATC
ENST000000002   foo foo foo foo 
```

Note the interleaved row for each gene. This is used by SAIGE to define a mask for each variant. We create a dummy value
so we can define our own masks.

Files are created per-chromosome to enable faster parallelization on DNA Nexus.

##### Command Line Example

SAIGE-GENE+ proceedes in two steps:

1. Fitting the null GLMM (done for the whole genome):

```commandline
step1_fitNULLGLMM.R
          --phenoFile=phenotypes_covariates.formatted.txt /                             # formated phenotype + covariate file generated during step 1
          --phenoCol=<pheno_name> /                                                     # phenotype name extracted from the provided phenotype file
          --isCovariateTransform=FALSE /                                                # Should SAIGE inverse normalise covariates (NO!)
          --sampleIDColinphenoFile=IID /                                                # Sample ID name in our covariates file
          --outputPrefix=<phenoname>.SAIGE_OUT /                                        # Output name for step 1
          --sparseGRMFile=sparseGRM_200K.sparseGRM.mtx /                                # Sparse GRM pre-computed during mrcepid-buildgrms
          --sparseGRMSampleIDFile=sparseGRM_200K.sampleIDs.txt /                        # Associated sample file for the GRM
          --nThreads=64 /                                                               # Number of threads (default 64)
          --LOCO=FALSE /                                                                # Should we do leave-one-chrom-out (NO!)
          --skipModelFitting=FALSE /                                                    # Should we skip model fitting and go straight to step2 (NO!)
          --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age,sex,wes_batch /   # Covariates to control for. Note – Sex is excluded when testing one sex (see source code)
          --qCovarColList=wes_batch                                                     # Categorical covariates to control for from the overal covarColList
          --traitType=binary/quantitative \                                             # set according to type of trait as provided by user
          --useSparseGRMtoFitNULL=TRUE                                                  # Use the sparse GRM to fit the null model rather than estimating the variance ratio from genetic data (massively improves runtime)
```     

**Note:** I have shortened the name of the sparseGRMFile for readability (see source code).

2. Performing rare variant burden tests (done per chromosome):

```commandline
step2_SPAtests.R
          --vcfFile=saige_input.vcf.gz \                                          # Input vcf file I document above
          --vcfField=GT \                                                         # Hardcoded INFO field to check for presence absence in the vcf file
          --GMMATmodelFile=SAIGE_OUT.rda \                                        # File generated by step1 above
          --sparseGRMFile=sparseGRM_200K.sparseGRM.mtx /                          # Sparse GRM pre-computed during mrcepid-buildgrms
          --sparseGRMSampleIDFile=sparseGRM_200K.sampleIDs.txt /                  # Associated sample file for the GRM
          --LOCO=FALSE \                                                          # Should we do leave-one-chrom-out (NO!)
          --SAIGEOutputFile=<file_prefix>.<chromosome>.SAIGE_OUT.SAIGE.gene.txt \ # Output file from this step
          --groupFile=<file_prefix>.<chromosome>.SAIGE.groupFile.txt \            # Input groupFile I document above
          --is_output_moreDetails=TRUE                                            # Output additional information in the output file (het/hom counts, carrier status, etc.)
          --maxMAF_in_groupTest=0.5 \                                             # Minimum allele count to include a variant. 0.5 means we include all variants (even singletons)
          --maxMissing=1                                                          # We define our own missingness, so set a value to ignore
          --chrom=<chromosome>                                                    # Dummy chromosome value
          --annotation_in_groupTest=foo                                           # Dummy value from the groupFile
```

**Note:** I have shortened the name of the sparseSigmaFile for readability (see source code).

#### STAAR

STAAR is unique here in that it does not have a command-line tool associated with it and is instead provided as an R package.
To run STAAR, we created two R wrapper scripts that 

1. Generates a null model (`runSTAAR_Null.R`).
2. Ingests the data required to run burden tests, and formats the phenotypic and genetic data into tables / matrices that are compatible with
the STAAR package (`runSTAAR_Genes.R`).
   
These commented scripts are provided in this repository at:

`resources/usr/bin/`

##### Inputs

STAAR requires the following inputs:

1. A sparse matrix created by the R package 'Matrix'. This file is created in the mrcepid-collapsevariants step of 
   this workflow and saved in the binary .rds format. Specifically, it a lgCMatrix object that stores only information for
   individuals that have given variant as a value of "1" or "2" depending on het/hom status.
   
```text
        chr1:1000000:A:T chr1:1000010:T:G chr1:1000020:G:A chr1:2000000:G:C chr1:2000030:A:C
1000000                .                .                .                .                1
1000001                .                2                .                .                .
1000002                .                .                1                .                .
1000003                .                .                .                .                1
1000004                1                .                .                2                . 
```

Column identifiers are _identical_ to that found in the bgen file created for BOLT. Note the differences from the input
[provided to SAIGE-GENE+](#saige-gene).

2. A tab-delimited file of variant information with a header. varID is identical to that stored for SAIGE:

```text
varID   chrom   pos ENST    column
1:1000000:A:T   1   1000000 ENST00000000001 1
1:1000010:T:G   1   1000010 ENST00000000001 2
1:1000020:G:A   1   1000020 ENST00000000001 3
1:2000000:G:C   1   2000000 ENST00000000002 4
1:2000030:A:C   1   2000030 ENST00000000002 5
```

3. A sparse GRM. 

Currently, we use the same GRM created for SAIGE as part of the [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms). 
This input is hard-coded into the `runSTAAR.R` script.

##### Command line example

This example command-line is to run both scripts that we have created for this applet for one chromosome. All commands are
space delimited:

```commandline
Rscript runSTAAR_Null.R \
          phenotypes_covariates.formatted.txt \               # Formated phenotypes/covariates to derive the NULL STAAR model
          <pheno_name> \                                      # phenotype name extracted from the provided phenotype file
          <is_binary> \                                       # Is the trait binary (false for quantitative)?
          <quant_covars> \                                    # Names of additional quantitative covariates to include in model (NULL for none)
          <catagorical_covars>                               # Names of additional catagorical covariates to include in model (NULL for none)

Rscript runSTAAR_Genes.R \                                     
          <file_prefix>.<chr>.STAAR.matrix.rds \              # File (1) from above
          <file_prefix>.<chr>.variants_table.STAAR.tsv \      # File (2) from above
          <pheno_name>.STAAR_null.rds                         # Output null model from runSTAAR_Null.R
          <pheno_name> \                                      # phenotype name extracted from the provided phenotype file
          <file_prefix> \                                     # prefix from mrcepid-collapsevariants
          <chr>                                               # Current chromosome we are running to name output files
```

Please see the source code for the R scripts cited above, but in brief, we first fit a null model (in `runSTAAR_Null.R)`:

```r
obj_nullmodel <- fit_null_glmmkin(formated.formula, data=data_for_STAAR, id="FID", family=binomial(link="logit"), kins = sparse_kinship)
```

and then run a loop that tests the effect of having a rare variant on the given phenotype using STAAR (in `runSTAAR_Genes.R)`:

```r
for (gene in genes) {
  # We subset the genotype matrix (genotypes; file (1) above) using the rownumbers in our variant file (file (2) above)
  current_GENE <- Matrix(genotypes[,variants[GENEID == gene, rownum]])
  # And then run STAAR using this subset matrix
  staar_result <- STAAR(genotype = current_GENE, obj_nullmodel = obj_nullmodel, rare_maf_cutoff = 1)
}
```

#### REGENIE

##### Inputs

Like BOLT, REGENIE requires the genetic data to first estimate a null model and control for relatedness between study 
participants. It then uses the derived null model to perform rare variant burden tests individually on masks/chromosome
combinations. Specific inputs for REGENIE are derived from the filtered WES data listed in `project-G6BJF50JJv8p4PjGB9yy7YQ2:file-G86GJ3jJJv8fbXVB9PQ2pjz6`.
Specifically, for burden tests, REGENIE requires:

1. A bgen format file of variants. Unlike for the other tools, we do not create this file during mrcepid-collapsevariants
and instead directly use the filtered bgen files listed in  `project-G6BJF50JJv8p4PjGB9yy7YQ2:file-G86GJ3jJJv8fbXVB9PQ2pjz6`,
filtered to individuals in accordance with inclusion/exclusion lists. This is due to how REGENIE handles variant masks.

2. An annotation file. This file lists the variant ID, identical to that provided to SAIGE-GENE, a gene, and the mask name
provided by the collapsevariants tarball. No header is required:

```text
1:1000000:A:T   ENST00000000001 HC_PTV-MAF_01
1:1000010:T:G   ENST00000000001 HC_PTV-MAF_01
1:1000020:G:A   ENST00000000001 HC_PTV-MAF_01
1:2000000:G:C   ENST00000000002 HC_PTV-MAF_01
1:2000030:A:C   ENST00000000002 HC_PTV-MAF_01
```

3. A gene-to-variant file. This file is almost identical to that provided to SAIGE-GENE, except with additional columns for
chromosome and gene start coordinate, and the list of variant IDs in comma-separated formated:

```text
ENST000000001 1 12345   1:1000000:A:T,1:1000010:T:G,1:1000020:G:A
ENST000000002 2 67890   1:2000000:G:C,1:2000030:A:C,1:2000050:ATC:A,1:2000000:G:GATC 
```

4. A file with one row that provides the name of the mask as listed in file (2) from above. 

```text
HC_PTV-MAF_01   HC_PTV-MAF_01
```

##### Command line example

Step One of REGENIE uses a command like:

```commandline
regenie  \                                                                  
  --step 1  \                                                               # Indicates to REGENIE to run step 1 of 2
  --bed /test/genetics/UKBB_450K_Autosomes_QCd_WBA  \                       # UKBB genetic data filtered to individuals analysed in this specific run
  --covarFile /test/phenotypes_covariates.formatted.txt  \                  # The processed and formated covariate/pheno file from above
  --phenoFile /test/phenotypes_covariates.formatted.txt  \                  # The processed and formated covariate/pheno file from above (required to be listed twice by REGENIE)
  --extract /test/REGENIE_extract.snplist  \                                # A set of variants with MAC > 100 & < ((N_Indv*2) - 100). The latter is because plink does not appear to calculate AC based on the minor allele, but rather on the alternate allele
  --bsize 100  \                                                            # REGENIE computation parameter
  --out /test/fit_out  \                                                    # Name of the file that contains the null model for REGENIE
  --threads 64 \                                                            # Number of threads. Default is to use 64, will be modified by the instance type selected
  --phenoCol <phenoname> \                                                  # Name of the phenotype in covarFile/phenoFile
  --covarColList PC{1:10},age,age_squared,sex \                             # Quantitative covariates to include. Only standard default covariates are listed here.
  --catCovarList wes_batch \                                                # Categorical covariates to include. Only standard default covariates are listed here.
  --bt                                                                      # Indicates running a binary trait. Only included for binary phenotypes
```

Step Two of REGENIE uses a command like:

```commandline
regenie  \                                                                                  
  --step 2  \                                                                               # Indicates to REGENIE to run step 2 of 2    
  --bgen /test/<chromosome>.markers.bgen  \                                                 # QCd bgen file                                    
  --sample /test/<chromosome>.markers.bolt.sample  \                                        # Matching .sample file for file provided by .bgen                                                 
  --covarFile /test/phenotypes_covariates.formatted.txt  \                                  # File identical to that provided to step 1 above                                               
  --phenoFile /test/phenotypes_covariates.formatted.txt  \                                  # File identical to that provided to step 1 above                                                                                                            
  --firth --approx  \                                                                       # Use Firth regression to calculate p. values with the approximation speedup (--approx)            
  --pred /test/fit_out_pred.list  \                                                         # Output file from step 1                        
  --anno-file /test/<tarball_prefix>.<chromosome>.REGENIE.annotationFile.tsv  \             # Annotations file (see inputs section above)                                                                            
  --set-list /test/<tarball_prefix>.<chromosome>.REGENIE.setListFile.tsv  \                 # Set list file (see inputs section above)                                                                       
  --mask-def /test/<tarball_prefix>.<chromosome>.REGENIE.maskfile.tsv  \                    # Mask definition file (see inputs section above)                                                                        
  --aaf-bins 1  \                                                                           # Tells REGENIE to include ALL variants (MAF < 100%) as we define the MAF bin ourselves       
  --vc-tests skato-acat,acato-full  \                                                       # Provide p. values for skat-o and acato-full (see REGENIE full documentation for more information)                       
  --bsize 200  \                                                                            # REGENIE computation parameter
  --threads 1  \                                                                            # Run 1 thread per step 2 job      
  --covarColList PC{1:10},age,age_squared,sex \                                             # Quantitative covariates to include. Only standard default covariates are listed here. We have to include again or REGENIE tries to include the phenotype as a covariate and fails.
  --catCovarList wes_batch \                                                                # Categorical covariates to include. Only standard default covariates are listed here.  
  --out /test/<output_prefix>.<tarball_prefix>.<chromosome>                                 # Name the outfile                                                                
```

REGENIE also (when requested with the `run_marker_tests` flag) performs per-marker tests. These are run identically to step 2
as shown above, without mask and annotation definitions.

#### GLMs

The method to perform Generalised Linear Models (GLMs) as part of this applet is implemented using the [statsmodels](https://www.statsmodels.org/stable/index.html)
python package.

##### Inputs

Variant data is identical to that for [STAAR](#staar). To ingest this data into Python, we convert the R object for 
staar into a sparse matrix using the script `sparseMatrixProcessor.R` in the `general_utilities.linear_model.R_resources` 
package. This script creates a tab-delimited file of all non-reference genotypes with ENST information like:

```text
FID varID   gt  ENST
1000000 1:12345:A:T 1   ENST000000001
1000001 1:12345:A:T 2   ENST000000001
1000003 1:12367:G:C 1   ENST000000002
```

Where 'gt' is 1 if an individual is heterozygous and 2 if an individual is homozygous.

##### Command line example

This method has no command line as it is run within the applet's source code. For specifics on the method used, please see 
the applet source code and the [statsmodels glm README](https://www.statsmodels.org/stable/glm.html). Briefly, this runs in 
two basic steps:

```python
import pandas as pd
import statsmodels.api as sm

# 1. read file (1) from above as a pandas dataframe:
geno_table = pd.read_csv("lm.tsv",
                         sep = "\t",
                         names = ['FID', 'varID', 'gt', 'ENST'])

genes = ["ENST0000000001", "ENST0000000002"]
is_binary = True
pheno_name = 'T2D'
pheno_covars = pd.DataFrame() # This is actually filled with covariate and phenotype information per-participant

# 2. Loop through all genes from file (lm.tsv):
for gene in genes:
    if is_binary == True:
       family = sm.familes.Binomial()
    else:
       family = sm.families.Gaussian()
    
    sm.GLM.from_formula(pheno_name + ' ~ has_var + sex + age + batch + wes_batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10', 
                        data=pheno_covars,
                        family=family).fit()
```

The 'has_var' variable in the generalised linear model is an additive variable where individuals are coded as 0, 1, 2, 3, ...
depending on the number of variants they have in a given gene/gene set.

## Running on DNA Nexus

### Inputs

These inputs are specific to this module. If the option is not required, defaults for 
each option are provided in **[bold brackets]**. Boolean options are flags, and change to 'true' when provided. 

| input                | Boolean? | Required? | description                                                                                                                                                                                                                 |
|----------------------|----------|-----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| association_tarballs | False    | **True**  | Hash ID(s) of the output from [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants) that you wish to use for rare variant burden testing. See below for more information.                     |
| tool                 | False    | **True**  | Tool to use for the burden testing module. **MUST** be one of 'bolt', 'saige', 'staar', 'glm', or 'regenie'. Case must match.                                                                                               |
| run_marker_tests     | **True** | False     | run SAIGE/BOLT/REGENIE per-marker tests? Note that tests with SAIGE currently take a VERY long time. **[False]**                                                                                                            |
| bgen_index           | False    | **True**  | index file with information on filtered and annotated UKBB variants                                                                                                                                                         |
| array_bed_file       | False    | **True**  | plink .bed format file from UKBB genetic data, filtered according to [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms)                                                                                  |
| array_fam_file       | False    | **True**  | corresponding .fam file for 'bed_file'                                                                                                                                                                                      |
| array_bim_file       | False    | **True**  | corresponding .bim file for 'bed_file'                                                                                                                                                                                      |
| low_MAC_list         | False    | **True**  | list of low MAC (<100) variants in 'bed_file'                                                                                                                                                                               |
| sparse_grm           | False    | **True**  | a sparse GRM for all individuals in 'bed_file' provided by [Bycroft et al.](https://www.nature.com/articles/s41586-018-0579-z)                                                                                              |
| sparse_grm_sample    | False    | **True**  | corresponding samples in 'sparse_grm'                                                                                                                                                                                       |
| bolt_non_infinite    | **True** | False     | Should BOLT be run with the flag `--lmmForceNonInf`? Only affects BOLT runs and may substantially increase runtime. **[False]**                                                                                             |
| regenie_smaller_snps | False    | False     | Run step1 of REGENIE with the smaller set of relatedness SNPs? This file is typically located at: `/Bulk/Genotype Results/Genotype calls/ukb_snp_qc.txt`. Only affects REGENIE runs and may substantially decrease runtime. |

#### Association Tarballs

`association_tarballs` takes either a single DNANexus file-ID that points to a single output from the 'mrcepid-collapsevariants' tool **OR** 
a file with multiple such DNANexus file-IDs. This file is one file-ID per line, like:

```text
file-G7z31B0J6F3ZbpGK6J0y2xxF
file-G7z2zP8JQy2PjKjY97Zvb3z9
file-G7z2pXjJ91Jx1ZJ4335bX52Z
file-G7z2gBjJ9ZpVqP3098X3xK6Y
file-G7z2VZQJ845v91751z2k9v2B
file-G7yv39QJPY6bYJKY9776JJBG
file-G7yqY10JkX6Y4fzvJf80z7P6
file-G7yq4ZjJV8qPjKjY97ZvZJ85
file-G7ypg98Jq133zgzY1yVy8gFz
```

Where each line is a different mask. An example file that includes 16 variant masks (8 variant types at two MAF cutoffs) 
is available at `collapsed_variants_new/variant_mask_list.txt (file-G7zPvZ0JJv8v06j8Gv2ppxpJ)` in project 
`project-G6BJF50JJv8p4PjGB9yy7YQ2`. Individual masks are available in `collapsed_variants_new/`. Please see the 
[high-level documentation](https://github.com/mrcepid-rap#collapsed-variants) for all apps for more information on
pre-collapsed variant files.

### Outputs

1. `<file_prefix>.genes.<TOOL>.stats.tsv.gz` (per-gene output)
2. `<file_prefix>.genes.<TOOL>.stats.tsv.gz.tbi` (per-gene output index)
3. `<file_prefix>.marker.<TOOL>.stats.tsv.gz` (per-marker output [when requested for BOLT / SAIGE / REGENIE])
4. `<file_prefix>.marker.<TOOL>.stats.tsv.gz.tbi` (per-marker output index [when requested for BOLT / SAIGE / REGENIE])

Note that some tools provide additional log/stat files that are not documented here, but are discussed in 
tool-specific documentation.

#### Per-gene output

A tab-delimited, gzipped file named like `<output_prefix>.genes.<tool>.stats.tsv.gz` (where `<output_prefix>` is identical
to that provided to the `output_prefix` input and `<tool>` is the name of the tool requesed by the `tool` input parameter) 
containing per-gene burden tests. An index for easy querying with tabix is also provided (`<output_prefix>.genes.<tool>.stats.tsv.gz.tbi`). 
Columns include those in the standard tool output. Additional columns contain per-gene information derived from in the
file `transcripts.tsv.gz (file-G7xyzF8JJv8kyV7q5z8VV3Vb)` in project `project-G6BJF50JJv8p4PjGB9yy7YQ2`. These columns include:
   
| column name       | description                                                                                                                                                                                                                                                                   |
|-------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ENST              | ENSMBL ENST ID. Will normally be the ENST ID that corresponds to MANE transcript or the ENSEMBL canonical transcript except in rare cases. When running in SNP or GENE list mode, this column will contain the dummy value of ENST0000000000 or ENST99999999999, respectively |
| chrom             | chromosome of this gene *without* the 'chr' prefix                                                                                                                                                                                                                            |
| start             | transcription start coordinate in hg38                                                                                                                                                                                                                                        |
| end               | transcription stop coordinate in hg38                                                                                                                                                                                                                                         |
| ENSG              | ENSEMBL ENSG corresponding to ENST                                                                                                                                                                                                                                            |
| MANE              | MANE v0.93 transcript                                                                                                                                                                                                                                                         |
| transcript length | end - start                                                                                                                                                                                                                                                                   |
| SYMBOL            | HGNC gene name                                                                                                                                                                                                                                                                |
| CANONICAL         | Is ENST the ENSEMBL canonical transcript?                                                                                                                                                                                                                                     |
| BIOTYPE           | Should *always* be protein_coding                                                                                                                                                                                                                                             |
| cds_length        | translation stop - translation start accounting for intron length                                                                                                                                                                                                             |
| coord             | formatted 'chrom:start-end' with chr prefix for aid in lookup on databases                                                                                                                                                                                                    |
| manh.pos          | relative position in the genome on a scale of 0-1 (chr1:1 = 0, chrY:14522573 = 1) for easy manhattan plot creation                                                                                                                                                            |

**Note:** If a transcript ID is NOT found in `transcripts.tsv.gz`, this information **WILL NOT** be included in final
output.

Additional columns derived from the prefix of input tarfiles from `collapsevariants` will also be included. If the 
tar is named like "HC_PTV-MAF_01", these additional columns will be headed as 'MASK' and 'MAF' respectively. This 
functionality is trigged by the second value delimited by '-' having the prefix of 'MAF' or 'AC' (e.g. MAF_01, MAF_1,
MAF_005, etc). Otherwise, will be labelled as `var1`, `var2`, etc. as DELIMITED by '-'. For example, if the tarfile 
name is "Foo-Bar.tar.gz", column `var1` will include "Foo" and column `var2` will include "Bar". The software 
currently does not have a method for naming these columns otherwise except for the special case as mentioned above.

#### Per-marker output

A tab-delimited, gzipped file named like `<output_prefix>.markers.<tool>.stats.tsv.gz` (where `<output_prefix>` is identical
to that provided to the `output_prefix` input and `<tool>` is the name of the tool requesed by the `tool` input parameter) 
containing per-marker burden tests. An index for easy querying with tabix is also provided 
(`<output_prefix>.genes.BOLT.stats.tsv.gz.tbi`). Additional columns contain per-marker information contained in the
file `450k_vep.sorted.tsv.gz (file-G857Z4QJJv8x7GXfJ3y5v1qV)` in project `project-G6BJF50JJv8p4PjGB9yy7YQ2`. These columns 
are identical to those provided by [mrcepid-annotatecadd](https://github.com/mrcepid-rap/mrcepid-annotatecadd#outputs).
Note that this output is only produced for BOLT, SAIGE, or REGENIE, where requested.

## Example Command

This is a module for the mrcepid-runassociationtesting app. Example command to run a BOLT burden test:

```commandline
dx run app-mrcepid-runassociationtesting --priority low --destination results/ \ 
        -imode=burden
        -ioutput_prefix="T2D.bolt" \
        -iinput_args=' \
                --tool bolt \
                --run_marker_tests \
                --association_tarballs file-1234567890ABCDEFGHIJKLMN \
                --phenofile file-1234567890ABCDEFGHIJKLMN \
                --is_binary \
                --inclusion_list file-1234567890ABCDEFGHIJKLMN \
                --sex 2 \
                --transcript_index file-1234567890ABCDEFGHIJKLMN \
                --base_covariates file-1234567890ABCDEFGHIJKLMN \
                --bgen_index file-1234567890ABCDEFGHIJKLMN \
                --bed_file file-1234567890ABCDEFGHIJKLMN \
                --fam_file file-1234567890ABCDEFGHIJKLMN \
                --bim_file file-1234567890ABCDEFGHIJKLMN \
                --low_MAC_list file-1234567890ABCDEFGHIJKLMN \
                --sparse_grm file-1234567890ABCDEFGHIJKLMN \
                --sparse_grm_sample file-1234567890ABCDEFGHIJKLMN`
```