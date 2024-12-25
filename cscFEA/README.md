# Cell-type-specific scFEA (cscFEA)

Our custom-made metabolic flux simulation was utilized to reveal altered flux between samples obtained at pre-flight and return time points. It is constructed based on scFEA by applying Recon3D metabolic model. The single cell Flux Estimation Analysis (scFEA) provided a cell-wise metabolic flux model employing graph neural networks by optimizing constraint-based loss functions. Based on the scFEA backbone, we applied additional inputs for reflecting more information of metabolisms and reactions: tissue specificity, gene reaction rule, and direction of reaction. Details are available in manuscript.

# Requirements

cscFEA has its base on Python3 and we suggest installing it in separated conda environment.
Installation example:
```
conda create -n cscFEA python=3.8
```

And install following requirements:
```
matplotlib
numpy
pandas
torch
tqdm
```

Magic is an imputation tool and optional.
If you want to use it:
```
pip install magic-impute
```
or if you do not want to use it, you can remove a following line from src/cscFEA.py:
```
import magic
```


# Usage

You can see the input arguments for cscFEA by help option:

```
python src/cscFEA.py --help
usage: cscFea.py [-h] [--data_dir <data_directory>]
                 [--input_dir <input_directory>] [--res_dir <data_directory>]
                 [--test_file TEST_FILE]
                 [--stoichiometry_matrix STOICHIOMETRY_MATRIX]
                 [--cName_file CNAME_FILE] [--sc_imputation {True,False}]
                 [--output_flux_file OUTPUT_FLUX_FILE]
                 [--output_balance_file OUTPUT_BALANCE_FILE]
                 [--train_epoch [TRAIN_EPOCH]]
                 [--gene_reaction_rule_path <data_directory>]
                 [--reaction_direction_path <data_directory>]
                 [--scoring {True,False}]

cscFEA: A graph neural network model to estimate cell type-specific cell-wise
metabolic flux using single cell RNA-seq data

options:
  -h, --help            show this help message and exit
  --data_dir <data_directory>
                        The data directory for cscFEA model files.
  --input_dir <input_directory>
                        The data directory for single cell input data.
  --res_dir <data_directory>
                        The data directory for result [output]. The output of
                        cscFEA includes two matrices, predicted metabolic flux
                        and metabolites stress at single cell or single
                        cluster resolution.
  --test_file TEST_FILE
                        The test SC file [input]. The input of cscFEA is a
                        single cell profile matrix, where row is gene and
                        column is cell. An example dataset is provided in
                        /data/ folder. The input can be raw counts or
                        normalised counts. The logarithm would be performed if
                        value larger than 50.
  --stoichiometry_matrix STOICHIOMETRY_MATRIX
                        The table describes relationship between compounds and
                        reactions. Each row is an intermediate metabolite and
                        each column is metabolic reaction. Stoichiometry
                        matrix of OXPHOS from Recon3D is provided in /data/
                        folder.
  --cName_file CNAME_FILE
                        The name of compounds. The table contains two rows.
                        First row is compounds name and second row is
                        corresponding id.
  --sc_imputation {True,False}
                        Whether perform imputation for SC dataset.
  --output_flux_file OUTPUT_FLUX_FILE
                        User defined predicted flux file name.
  --output_balance_file OUTPUT_BALANCE_FILE
                        User defined predicted balance file name.
  --train_epoch [TRAIN_EPOCH]
                        User defined EPOCH (training iteration).
  --gene_reaction_rule_path <data_directory>
                        The table of gene reaction rules. It contains two
                        columns. The first column is name of reaction, and the
                        second column is corresponding reaction formula. Genes
                        linked with OR are separated into different rows with
                        same name. Genes linke with AND are written, e.g. A
                        and B and C. For the details, please check enclosed
                        sample.
  --reaction_direction_path <data_directory>
                        The table of directionalities of reactions. It
                        contains two columns. The first column is name of
                        reaction, and the second column is corresponding
                        directionality based on Recon3D. 1 for forward only, 0
                        for bidirectional, -1 for backward reaction only.
  --scoring {True,False}
                        To ease visualizing relative differences between
                        groups, grouping, calculating mean, normalizing, and
                        scoring are introduced.
```

Run code with default parameters:

```
python src/cscFEA.py
```


### References
##### scFEA
N. Alghamdi, W. Chang, P. Dang, X. Lu, C. Wan, Z. Huang, J. Wang, M. Fishel, S. Cao, C. Zhang. scFEA: A graph neural network model to estimate cell-wise metabolic using single cell RNA-seq data, under review at Genome Research, 2020.

##### Recon3D
King ZA, Lu JS, Dräger A, Miller PC, Federowicz S, Lerman JA, Ebrahim A, Palsson BO, and Lewis NE. BiGG Models: A platform for integrating, standardizing, and sharing genome-scale models (2016) Nucleic Acids Research 44(D1):D515-D522. doi:10.1093/nar/gkv1049
