# Baseline FA\*IR

This folder includes the source code and datasets for baseline FA\*IR.

## Requirements

This program was developed and tested in [Python 3.5](https://www.python.org/downloads/release/python-350/). It depends on the following modules:

* cycler 0.10.0
* DateTime 4.1.1
* guacamole 0.9.2
* Jinja2 2.9.6
* MarkupSafe 1.0
* matplotlib 2.0.0
* numpy 1.12.0
* padme 1.1.1
* pandas 0.19.2
* pip 9.0.1
* plainbox 0.34.0
* pyparsing 2.1.10
* PyPDF2 1.26.0
* python-dateutil 2.6.0
* pytz 2016.10
* requests 2.13.0
* scipy 0.18.1
* setuptools 28.8.0
* six 1.10.0
* utils 0.9.0
* zope.interface 4.3.3

To install all the dependencies, use the following command:
```bash
pip install -r requirements.txt
```

## Usage

1. Unzip the dataset file `data_fairtq_processed.zip` into the current folder.
2. Go into the source code folder `cd src/`.
3. Run the code with specific dataset. For example, the below command will execute all the experiments on dataset XING:
```bash
python main.py -c xing
````
4. The results of an execution will be saved in folder `results/`.
