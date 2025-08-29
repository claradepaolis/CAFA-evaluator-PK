# CAFA-evaluator (CAFA5 Edition)

This repo was branched from the [CAFA-evaluator](https://github.com/BioComputingUP/CAFA-evaluator) and additional functionality was added.
Visit the original [CAFA-evaluator Wiki](https://github.com/BioComputingUP/CAFA-evaluator/wiki) for more information about the algorithm.


The two main new functionalities are:
* **Flexible terms-of-interest**
  A file can be passed in with a list of ontology terms that will be evaluated for all proteins. All terms no included in the file
  will not be evaluated. This can be used to exclude terms that were added to the ontology since predictions were collected or to
  exclude terms that have been deemed uninformative. In the figure below, the terms of interest are shown with a red outline.
  These terms will be evaluated for _all proteins_ if they appear in the ground truth file.
* **Protein-specific known annotations**
  To evaluate under the "Partial Knowledge" evaluation setting, any annotations known previous to prediction should be excluded from
  evaluation. This exclusion is done similarly to the terms-of-interst, but for each individual protein.
  In the figure below, the annotations are shown for a single protein. Terms with annotations known before the prediction phase are
  shown in yellow in the Cellular Component and Molecular Function aspects of the Gene Ontology. These terms should be excluded from
  evaluation. New annotations are shown to the right in blue. Evaluation should happen _only_ for these terms. If the newly annotated
  terms do not appear in the terms-of-interest file, they will not be evaluated.

 ![Evaluting CAFA Partial Knowledge setting and using terms of interest](EvalutingCAFAPartial.jpg)



## Citation
Please cite the original CAFA Evaluator and the forthcoming CAFA 5 papers if you use this code in published research  
[CAFA-evaluator: A Python Tool for Benchmarking Ontological Classification Methods](https://doi.org/10.1093/bioadv/vbae043)  
*D Piovesan, D Zago, P Joshi, MC De Paolis Kaluza, M Mehdiabadi, R Ramola, AM Monzon, W Reade, I Friedberg, P Radivojac, SCE Tosatto*  
**Bioinformatics Advances (2024)** - DOI: [10.1093/bioadv/vbae043](https://doi.org/10.1093/bioadv/vbae043)

Crowdsourcing the fifth critical assessment of protein function annotation algorithms (CAFA 5) yields improvement in protein function prediction  
*TBD*  
**TBF** - DOI: TBD


## Usage

The program can be executing the command line interface or as a library.
Both the command line and the library accept the following positional arguments:

* **Ontology file** in OBO format

* _**NEW**_: **Terms of Interst file** contains term names that appear in the OBO file to be included in evaluation
  
* _**NEW**_: **Known annotations file** contains annotations known before the evaluation phase

* **Prediction folder** contain prediction files. Files can be organized into sub-folders, 
sub-folders are processed recursively and the sub-folder name is used as prefix for the method name

* **Ground truth file** containing targets and associated ontology terms

Example input files are provided inside the `data/example` folder. 

### Command line

When executed from the command line the script logs information about the calculation in the console (standard error) and
will create a folder named `results` containing the evaluation results. 
A different folder can be specified using the `-out_dir` option. 
 
The original CAFA-evaluator functionality works as before without additional imput arguments

```bashcon
python3 /path/to/CAFA-evaluator/src/cafaeval/__main__.py ontology_file prediction_folder ground_truth_file 
```

_**NEW**_: By default, all terms in the ontology will be considered in the evaluation. To include only specific terms, provide a terms-of-interest file
```bashcon
python3 /path/to/CAFA-evaluator/src/cafaeval/__main__.py ontology_file prediction_folder ground_truth_file -toi terms_of_interest_file
```

_**NEW**_: To evaluate Partial Knowledge annotations, the known annotations file must be passed in with the `-known` option.
```bashcon
python3 /path/to/CAFA-evaluator/src/cafaeval/__main__.py ontology_file prediction_folder ground_truth_file -known known_annotations_file
```

_**NEW**_: You can pass in both terms-of-interest and known annotations:
```bashcon
python3 /path/to/CAFA-evaluator/src/cafaeval/__main__.py ontology_file prediction_folder ground_truth_file -toi terms_of_interest_file -known known_annotations_file
```

### Library

The `cafa_eval` function is the main entry point of the package. 
Below is reported an example using the example files provided in the `data/example` folder.


```pycon
>>> import cafaeval
>>> from cafaeval.evaluation import cafa_eval
>>> cafa_eval("IDPO_disorder_function.obo", "predictions", "ground_truth.tsv")
(                                       n        tp        fp        fn        pr  ...         f         s  pr_micro  rc_micro   f_micro
filename   ns                tau                                                  ...                                                  
pred_5.tsv disorder_function 0.01  168.0  2.928571  7.083333  0.130952  0.292532  ...  0.448380  7.084544  0.292509  0.957198  0.448087
                             0.02  168.0  2.928571  7.083333  0.130952  0.292532  ...  0.448380  7.084544  0.292509  0.957198  0.448087
                             0.03  168.0  2.928571  7.077381  0.130952  0.292695  ...  0.448571  7.078592  0.292683  0.957198  0.448292
                             0.04  168.0  2.928571  7.077381  0.130952  0.292695  ...  0.448571  7.078592  0.292683  0.957198  0.448292
                             0.05  168.0  2.928571  7.077381  0.130952  0.292695  ...  0.448571  7.078592  0.292683  0.957198  0.448292
...                                  ...       ...       ...       ...       ...  ...       ...       ...       ...       ...       ...
pred_1.tsv disorder_function 0.41    1.0  0.005952  0.017857  3.053571  0.250000  ...  0.003937  3.053624  0.250000  0.001946  0.003861
                             0.42    1.0  0.005952  0.017857  3.053571  0.250000  ...  0.003937  3.053624  0.250000  0.001946  0.003861
                             0.43    1.0  0.005952  0.017857  3.053571  0.250000  ...  0.003937  3.053624  0.250000  0.001946  0.003861
                             0.44    1.0  0.005952  0.017857  3.053571  0.250000  ...  0.003937  3.053624  0.250000  0.001946  0.003861
                             0.45    1.0  0.005952  0.011905  3.053571  0.333333  ...  0.003945  3.053595  0.333333  0.001946  0.003868

[352 rows x 14 columns], {'f':                                        n        tp        fp        fn        pr  ...         s  pr_micro  rc_micro   f_micro  cov_max
filename   ns                tau                                                  ...                                                 
pred_1.tsv disorder_function 0.04  166.0  1.744048  2.095238  1.315476  0.466566  ...  2.473964  0.454264  0.570039  0.505608      1.0
pred_2.tsv disorder_function 0.84  163.0  1.755952  1.845238  1.303571  0.504499  ...  2.259248  0.487603  0.573930  0.527256      1.0
pred_3.tsv disorder_function 0.89  168.0  2.113095  1.601190  0.946429  0.638889  ...  1.859983  0.568910  0.690661  0.623902      1.0
pred_4.tsv disorder_function 0.06  168.0  2.333333  0.666667  0.726190  0.777778  ...  0.985798  0.777778  0.762646  0.770138      1.0
pred_5.tsv disorder_function 0.38  167.0  2.345238  1.952381  0.714286  0.596671  ...  2.078941  0.545706  0.766537  0.637540      1.0

[5 rows x 15 columns], 's':                                        n        tp        fp        fn        pr  ...         s  pr_micro  rc_micro   f_micro  cov_max
filename   ns                tau                                                  ...                                                 
pred_1.tsv disorder_function 0.06  102.0  1.047619  0.297619  2.011905  0.816667  ...  2.033799  0.778761  0.342412  0.475676      1.0
pred_2.tsv disorder_function 0.91  124.0  1.196429  0.839286  1.863095  0.638978  ...  2.043410  0.587719  0.391051  0.469626      1.0
pred_3.tsv disorder_function 0.89  168.0  2.113095  1.601190  0.946429  0.638889  ...  1.859983  0.568910  0.690661  0.623902      1.0
pred_4.tsv disorder_function 0.06  168.0  2.333333  0.666667  0.726190  0.777778  ...  0.985798  0.777778  0.762646  0.770138      1.0
pred_5.tsv disorder_function 0.42  156.0  1.642857  0.714286  1.416667  0.776221  ...  1.586552  0.696970  0.536965  0.606593      1.0

[5 rows x 15 columns], 'f_micro':                                        n        tp        fp        fn        pr  ...         s  pr_micro  rc_micro   f_micro  cov_max
filename   ns                tau                                                  ...                                                 
pred_1.tsv disorder_function 0.04  166.0  1.744048  2.095238  1.315476  0.466566  ...  2.473964  0.454264  0.570039  0.505608      1.0
pred_2.tsv disorder_function 0.84  163.0  1.755952  1.845238  1.303571  0.504499  ...  2.259248  0.487603  0.573930  0.527256      1.0
pred_3.tsv disorder_function 0.89  168.0  2.113095  1.601190  0.946429  0.638889  ...  1.859983  0.568910  0.690661  0.623902      1.0
pred_4.tsv disorder_function 0.06  168.0  2.333333  0.666667  0.726190  0.777778  ...  0.985798  0.777778  0.762646  0.770138      1.0
pred_5.tsv disorder_function 0.38  167.0  2.345238  1.952381  0.714286  0.596671  ...  2.078941  0.545706  0.766537  0.637540      1.0

[5 rows x 15 columns]})
```

The output of `cafa_eval` is a tuple containing:
* A pandas DataFrame with the evaluation results, one row per prediction file, namespace and threshold.
* A dictionary with the best scores (max F-measure, max Weighted F-measure, min Semantic similarity). For each 
score the dictionary contain a pandas DataFrame with one row per prediction file, namespace and threshold. 

The `write_results` function generates the output files.

```pycon
>>> import cafaeval
>>> from cafaeval.evaluation import cafa_eval, write_results
>>> res = cafa_eval("IDPO_disorder_function.obo", "predictions", "ground_truth.tsv")
>>> write_results(*res)
```


## Input files
**Prediction file** - Tab separated file with the target ID, term ID and score columns.

~~~txt
T_1	IDPO:00501	0.06
T_1	IDPO:00506	0.05
T_1	IDPO:00507	0.03
T_2	IDPO:00501	0.04
T_2	IDPO:00506	0.02
...
~~~

**Ground truth file** - Tab separated file with the target ID and term ID. 
Additional columns are discarded.
~~~
T_1	IDPO:00024
T_2	IDPO:00506
T_3	IDPO:00502
T_4	IDPO:00025
...
~~~

**Information accretion file (optional)** - If not provided, the weighted and S statistics are not generated.
Information accretion (IA) can be calculated as described in
[Wyatt and Radivojac, Bioinformatics, 2013](https://pubmed.ncbi.nlm.nih.gov/23813009/) 
and implemented in [https://github.com/claradepaolis/InformationAccretion](https://github.com/claradepaolis/InformationAccretion)

```
IDPO:00024  6.32
IDPO:00506  12.04
IDPO:00502  1.34
IDPO:00025  0.56
...
```

_**NEW**_: **Known annotations (optional)** - Tab separated file with target ID and term ID. 
File containing known annotations to exclude from partial-knowledge evaluation. 
If not provided, all terms will be used in evaluation

```
A0A009IHW8	GO:0072523	BPO
A0A009IHW8	GO:0046700	BPO
A0A021WW32	GO:0048869	BPO
A0A021WW32	GO:0006996	BPO
...
```

_**NEW**_: **Terms of Interest (optional)** - File with term ID to include in evaluation for all proteins, one ID per line.  
If not provided, all terms will be used in evaluation. 
This file is used to specify terms that will be evaluated, usually used to exclude terms in the ontology that have since
been obsoleted or are not of interest for the evaluation.

```
GO:0055039
GO:0072523
GO:0003882
GO:0010139
...
```

## Output files

Output files are generated in the `results` folder. The same files are gerated by both
the command line and the `write_results` function.

* `evaluation_all.tsv` corresponds to the first object returned by the `cafa_eval` function.
* `evaluation_best_< metric >.tsv` corresponds to the second object returned by the `cafa_eval` function. 
A different file for each metric is created.

**Note**: Weighted scores are generated only if the *Information Accretion* file is provided.


## Optional parameters

|  Argument   | Default value | Description                                                                                                                                                                                                                                                               |
|:-----------:|------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|  -out_dir   |   'results' | Output directory (tsv files + log). Either relative to current path or absolute                                                                                                                                                                                           |
|     -ia     |            | Information accretion file                                                                                                                                                                                                                                                |
|   -known    |            | Known annotations for each protein that should be excluded in partial-knowledge evaluation                                                                                                                                                                               |
|    -toi     |            | Terms of interest file   (terms considered for all proteins)                                                                                                                                                                                                              |
| -no_orphans |  False (flag) | Exclude orphans nodes (e.g. roots) from the calculation                                                                                                                                                                                                                   |
|    -norm    |     'cafa'  | Normalization strategy. `cafa` normalize precision by the number of predicted targets and recall by the number of targets in the ground truth. `pred` normalize by the number of  predicted targets. `gt` normalize by the number of ground truth proteins                |
|    -prop    |     'max'  | Ancestor propagation strategy. `max` propagate the max score of the traversed subgraph iteratively. `fill` propagate with max until a different score is found                                                                                                            |
|  -th_step   |      0.01  | Step size of prediction score thresholds to consider in the range [0, 1). A smaller step, means more calculation                                                                                                                                                          |
| -max_terms  |            | Number of terms for protein and namespace to consider in the evaluation. Parsing stops when the target limit for every namespace is reached. The score is not checked, meaning that terms are not sorted before the check, and the check is performed before propagation. |
|  -threads   |       4    | Parallel threads. `0` means use all available CPU threads. Do not use multi thread if you are short in memory                                                                                                                                                             |

