# TRIAGE
## **T**hroughput __R__anking by __I__terative __A__nalysis of __G__enomic __E__nrichment)

TRIAGE is designed to analyze the user input for enrichment and predicted interactions using locally saved repositories from publicly available databases. The expected user input is a genome-scale list of gene candidates (~20,000 observations) with an integer based variable to break down the gene candidates into groups of high, medium, and low confidence according to a cutoff set by the user. The expected output is a subset of gene candidates suggested for further analysis (<1,000), enriched pathways of biological processes represented by the new data set, and HTML file of interactive graph for further analysis of novel and canonical candidates. Unless otherwise noted, all input and output files are in CSV format.


## Getting Started

TRIAG is available either as the [development version](https://triage.niaidawsqa.net) or/and [production version](https://triage.niaid.nih.gov/). TRIAGE can also be downloaded and run as a standalone application. The instructions below will get you a copy of TRIAGE up and running on your local machine. 

### Prerequisites

To run TRIAGE on your local machine, R and Rstudio should be installed and running on your machine. Instructions are available for [R installation](https://cran.r-project.org/bin/) and for [RStudio installation](https://www.rstudio.com/products/rstudio/download/). In addition, several R packages are also required for TRIAGE to work on your local machine. 

```
# To install required R packages from R commandline
install.packages('dplyr')
install.packages('leaflet')
install.packages('DT')   
install.packages('shinyjs')
install.packages('shinyBS')
install.packages('readr')
install.packages('stringi')
install.packages('reshape2')
install.packages('data.table')
install.packages('edgebundleR')
install.packages('igraph')
install.packages('shinyAce')
install.packages('rJava')
install.packages('mailR')
```

### Installing

Before installing TRIAGE, you need to install *git* by following [these instructions]: https://gist.github.com/derhuerst/1b15ff4652a867391f03. After installing *git*, you can install TRIAGE on your local machine from command line. You can install it in your home directory or anywhere under your home directory:

```
# For developmental version
git clone https://github.niaid.nih.gov/Signaling-Systems-Unit/TRIAGE.git

# Or for production version
TO BE ADDED HERE!!!
```
After this, you should see a TRIAGE directory which contains all required files and data to run TRIAGE on your local machine.

## Running TRIAGE

To run TRIAGE on your local machine, start RStudio, open 'app.R' file in your TRIAGE directory, and then click 'Run App' button. You should see TRIAGE running in your default web browser.


### Running Tests

A sample input file (`TRIAGEinput_HuTNF_CSAfdr_5percCO.csv`) can be used to test TRIAGE. The sample is a part of TRIAGE package and can be found at:

```
# Sample input file
TRIAGE/app/inputOutputs/TRIAGEinputFiles/TRIAGEinput_HuTNF_CSAfdr_5percCO.csv
```

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the current version available, see the [latest version of TRIAGE] https://github.niaid.nih.gov/Signaling-Systems-Unit/TRIAGE.git. 

## Authors

**Sam Katz** 
**Jian Song**
**Iain Fraser**
[Contact us] mailto: triage-team@nih.gov

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* We appreciate all the help we got from of the OEB Platform Team, NIH\NIAID\OCICB.
* We also want to thank LISB colleagues at NIH for their supports.


