# Baselines

This folder contains the source code, data and test samples of the baselines.

## Requirements
- Linux Environment with GCC for baselines Greedy, Prop and SR-Adapt written in C++
- For baseline FA*IR written in Python, please refer to [`FA-IR`](FA-IR/)

## Usage

1. Compile with `make` command.
   - This will automatically compile the dependency package qhull under qhull/ and compile all the C++ based baselines.
   - There is no need to compile baseline FA*IR, please refer to [`FA-IR`](FA-IR/) for how to run it.
2. Run our tests with `./run_all_tests.sh`.
   - This will run the two tests under folder test1/ and test2/, using dataset XING and COMPAS, respectively.
   - For both tests, baselines Greedy, Prop and SR-Adapt will be executed one-by-one.
   - After running both tests with `./run_all_tests.sh`, the result files (introduced [here](https://github.com/satansin/FairTQ/tree/main/ftq))
     will be generated under each folder.
   - To run each test individually, as shown in the scripts of run_all_tests.sh, follow the steps below (where "xx" should be replaced by the name of each C++ based baseline):
        - First, copy all files in the corresponding test folder to `xx/bin/` and go into folder `xx/bin/`.
        - Then, to run the each baseline algorithm, use `./ftq xx` (e.g., `./ftq Prop`).
        - After running each algorithm, the result file (to be introduced in detail later) will be generated under folder `FairTQ-Exact/bin/`.

The input format, output format and the steps to run each baseline on other datasets are the same as those introduced in [`ftq/`](../ftq).
