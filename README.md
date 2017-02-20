## Readme

FarmCPUpp provides an efficient reimplementation of the FarmCPU R scripts found [here](http://zzlab.net/FarmCPU/index.html). Through the use of the `bigmemory`, `Rcpp`, `RcppEigen`, and `RcppParallel` packages, FarmCPUpp decreases the memory usage of the original FarmCPU scripts and improves the runtime through parallelization. FarmCPUpp will provide the greatest benefits on large datasets and on machines with many CPU cores.

### Installation

FarmCPUpp relies on multiple HPC packages for efficient memory usage and parallel processing. To install FarmCPUpp itself, you will need the `devtools` package. Run this code to ensure that you have the appropriate packages installed and updated:

```
packages <- c("Rcpp", "RcppEigen", "RcppParallel", "parallel", "doParallel", 
              "foreach", "bigmemory", "devtools")
install.packages(packages)
devtools::install_github(repo = "amkusmec/FarmCPUpp")
```

Following installation the package can be loaded for use with

```
library(bigmemory)
library(FarmCPUpp)
```

### Data Formats

FarmCPUpp requires three data files with two optional files. Please note that the order of data in these files is important; results obtained will be incorrect if files are not ordered properly. Please see section [Data Dependencies](#data-dependencies) for more information.

#### Phenotype Data

Phenotype data should be provided as a two column dataframe where the first column contains sample names and the second column contains numeric phenotypic values. Missing values are allowed and should be specified as `NA`. The column name of the phenotypic data will be used for output file names. A dataframe with more than two columns may be provided to FarmCPUpp, but only the phenotype in the second column will be analyzed. The code below loads a sample phenotype file.

```
myY <- read.table(system.file("extdata", "mdp_traits_validation.txt",
                              package = "FarmCPUpp"),
                  header = TRUE, stringsAsFactors = FALSE)
```

#### Genotype Data

Genotype data should be provided in numerical format where 0 indicates no copies of the minor allele, 1 indicates a single copy of the minor allele, and 2 indicates two copies of the minor allele. Any value in the range [0,2] is accepted. The file contains marker scores in the columns and samples in the rows. The file should also include a header line, and the first column should contain taxa names. Note that the marker IDs from the genotype information file are used to redefine the header line. Because genotype files will often be very wide and large, the `bigmemory` package is used to import the data:

```
myGD <- read.big.matrix(system.file("extdata", "mdp_numeric.txt",
                                    package = "FarmCPUpp"),
                        type = "double", sep = "\t", header = TRUE,
                        col.names = myGM$SNP, ignore.row.names = FALSE,
                        has.row.names = TRUE)
```

If the same genotype file will be used in multiple GWASs, it may be advantageous to create a file-backed big matrix that can be reattached, since creating a big matrix from a large genotype file can be time consuming. The following code demonstrates how to create the file-backed big matrix and access it later.

```
# Load the data into a big.matrix.
# Note the use of the backingfile and descriptorfile arguments.
myGD <- read.big.matrix(system.file("extdata", "mdp_numeric.txt",
                                    package = "FarmCPUpp"),
                        type = "double", sep = "\t", header = TRUE,
                        col.names = myGM$SNP, ignore.row.names = FALSE,
                        has.row.names = TRUE, backingfile = "mdp_numeric.bin",
                        descriptorfile = "mdp_numeric.desc")

# Save the pointer for access later
dput(describe(myGD), "mdp_numeric_pointer.desc")

# The big.matrix can be reattached in a different R session using
desc <- dget("mdp_numeric_pointer.desc")
myGD <- attach.big.matrix(desc)
```

Please note that if you are using a file-backed big matrix, you must remain in the same working directory as the backing file until the GWAS is completed or the worker processes for bin selection will not be able to access the genotype data.

#### Genotype Information

Genotype information should be provided as a three column dataframe with colum names. The first column, SNP, should contain unique IDs for each marker. The second column, Chromosome, should contain integer or numeric IDs for each chromosome. The third column, Position, shoudl contain the integer base-pair positions of each marker. The code below loads a sample genotype information file.

```
myGM <- read.table(system.file("extdata", "mdp_SNP_information.txt",
                               package = "FarmCPUpp"),
                   header = TRUE, stringsAsFactors = FALSE)
```

#### Covariates

User-specified covariates can also be included in the GWAS. These may include principal components of the genotype data, results from programs such as STRUCTURE, or other experiment-related covariates that may be important for the GWAS. These data should be provided as a matrix with column names specifying the names of the covariates and row names specifying the sample names.

#### Genotype Prior

FarmCPUpp also accepts a dataframe of prior probabilities that a marker may be selected as a pseudo-QTN. These should be provided in a dataframe with the same format and information as the [genotype information](#genotype-information) with a fourth column named Probability that contains numeric probabilities.

#### Data Dependencies

FarmCPUpp assumes that the numbers of samples and markers are the same across all input data. Therefore, the number of rows in the genotype information should be equal to the number of rows of the prior and the number of columns of the genotype data, and all three should be sorted in the same order. Additionally, the number of rows of the phenotype data, covariates, and genotype data should be equal and sorted in the same order. This means that there may be samples in the genotype table that are not present in the phenotype table, for example. These samples should be added and missing values entered as appropriate. Missing values will be removed as necessary by FarmCPUpp.

### GWAS

GWAS is performed through one main function called `farmcpu`. This function requires phenotype data, genotype information, and genotype data. The simplest GWAS is run using

```
myResults <- farmcpu(Y = myY, GD = myGD, GM = myGM)
```

Please see the documentation (`?farmcpu`) for more information on the various arguments and options. By default FarmCPUpp will run on one core. The user may specify different numbers of cores to use for single-marker regressions and bin selection.

### Working with the Results

FarmCPUpp returns a results list. If `iteration.output = TRUE`, the list will contain one element for each iteration plus a final element with the same name as the phenotype that contains the final results. Each element of the results list is itself a list and may contain up to two elements. The first element, named GWAS, contains the GWAS results, including the information for each marker along with its effect estimate, standard error, t-statistic, and p-value. The second element, named CV, will contain a matrix with the covariate effect estimates if covariates were provided. The results object may be conveniently written to disk by

```
write_results(myResults)
```

FarmCPUpp also contains a function for creating manhattan plots from a single GWAS results table. The plotting colors and y-axis label may be customized. The user may also input a desired significance threshold. For example, to create a manhattan plot for the EarDia example data at a 0.01 significance threshold, use the following code.

```
manhattan_plot(myResults$EarDia$GWAS, cutoff = 0.01)
```

### Fine-tuning Performance

The speed-ups provided by this package are derived from two places: the use of C++ for single-marker regressions and parallel processing. The user has the option to control the degree of parallel processing used at two steps in the program in order to optimize performance through the use of the `ncores.glm` and `ncores.reml` parameters. There are several considerations for setting these parameters.

1. The maximum value for either parameter should be the number of cores present on your machine. Using more cores than you have available will not decrease runtime beyond the runtime achieved using the maximum number of cores available.
2. `ncores.reml` and `ncores.glm` are used at different points in the program. Therefore, the sum of these two parameters may be greater than the number of cores on your machine without affecting performance.
3. Overhead due to interprocess communication decreases the effectiveness of parallelization when the problem size is small. Therefore, single-marker regressions will not be run in parallel if the total number of markers is less than 10000.
4. The upper bound on number of pseudo-QTNs added in each iteration is $\frac{\sqrt{n}}{log_{10}\sqrt{n}}$, where $n$ is the number of non-missing phenotypic observations. Therefore, you can calculate the number of bin combinations you have prior to execution and set `ncores.reml` accordingly.
