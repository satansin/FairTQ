# Fair Top-k Query on Alpha-Fairness

This repository contains the code and the [technical report](FairTopK-TechnicalReport.pdf) for the submission "Fair Top-k Query on Alpha-Fairness".

## Requirements
- Linux Environment with GCC

## Usage

1. Download this repository

2. Change the working direction into FairTQ-Exact/.

3. Compile with `make` command.

    - This will automatically compile the dependency package qhull under qhull/ and compile our main code under UTK/.

4. Run our tests with `./run_all_tests.sh`.

    - This will run the two tests under folder test1/ and test2/, using dataset XING and COMPAS, respectively.
    - For both tests, the exact algorithm FairTQ-Exact will first be executed, and then follows the approximation algorithm beta-FairTQ-Appro.
    - After running both tests with `./run_all_tests.sh`, the result files (to be introduced in detail later) will be generated under each folder.
    - To run each test individually, as shown in the scripts of run_all_tests.sh, follow the steps below:
        - First, copy all files in the corresponding test folder to UTK/bin/.
        - Then, to run the exact algorithm, use `./ftq`, and to run the approximation algorithm, use `./ftq appr1`.
        - After running each algorithm, the result file (to be introduced in detail later) will be generated under folder UTK/bin/.

## Input Format

The detailed format of the input files are introduced as follows. The full examples can be found in any of the test folders.

- The config file (i.e., "config.txt") containing the input parameters (each in a row):
    - dataset file: path to the dataset file (e.g., "XING.txt")
    - n: size of the dataset (e.g., 2280)
    - k: parameter for top-k (e.g., 10)
    - d: number of dimensions (e.g., 2)
    - |P|: number of groups (e.g., 2)
    - alpha: parameter for alpha-fairness (e.g., 0.4)
- The dataset file (e.g., "XING.txt") containing a number of rows, each corresponding to a d-dimensional tuple. For each row:
    - The first d values show the values of this tuple in each dimension plus 0.0001 (e.g., 0.1161 0.0719)
    - The next d values show the values of this tuple in each dimension minus 0.0001 (e.g., 0.1159 0.0717)
    - The last value show the group membership of this tuple (e.g., 1 or 0)
- The input (initial) utility function (i.e., "originalF.txt") containing (d - 1) numbers, which are the values of the initial utility function in the first (d - 1) dimensions (e.g., 0.383139 0.191711 for d = 3). Note that all values of a utility function in the d dimensions sum up to 1, and thus the value of the last dimension can be easily computed.

## Output Format

The output of each of our algorithms contains the total time cost, the best utility function w*, the top-k set w.r.t. w*, the alpha-fairness of this top-k set and the modification penalty of w* to the initial utility function. An example of the output is shown as follows.
```
Total time cost: 0.040000 SEC
Best Utility Function: 0.657716 0.342284
Best Top-k: 108 1612 1652 585 1617 1816 1622 1831 665 2061
Best alpha-Fairness: 0.004082
Best Modification Penalty: 0.250512
```

Moreover, the output is also saved in the output files. Specifically, the time costs of the exact algorithm and the approximation algorithm are saved in "time_ftq_exact.txt" and "time_ftq_appr1.txt", respectively. The other results of the exact algorithm and the approximation algorithm are saved in statistic files "stat_ftq_exact.txt" and "stat_ftq_appr1.txt", respectively. In each statistic file, the four rows correspond to the best utility function w*, the top-k set w.r.t. w*, the alpha-fairness of this top-k set and the modification penalty of w* to the initial utility function, respectively. An example of the statistic file is shown as follows.
```
0.657716 0.342284 
108 1612 1652 585 1617 1816 1622 1831 665 2061 
0.004082
0.250512
```

## Additional Description

Additionally, we describe each folder in our code.

- qhull/: the qhull library (http://www.qhull.org/) for computing the convex hull which is a dependency of this project
- test1/ & test2/: the two tests
- UTK/: the folder of the main code of this project
<!-- - UTK/bin/: the folder for the compiled binary
- UTK/src/: the folder for the source code of this project
- UTK/src/headers/: the folder the C++ headers of this project -->
