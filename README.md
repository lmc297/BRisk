# BRisk (beta)
Exposure/risk assessment tool for *B. cereus* group isolates

## Overview

Example code for a *B. cereus* group risk assessment app

**This app is just an example of what a *B. cereus* group exposure/risk assesment app might look like. It is not an actual exposure/risk assessment, the parameters for all growth conditions are randomly generated, and the output of this demo app is completely meaningless. Please do not treat it as a real result.**

### Citation
If you find this code useful, please cite:
Carroll, Laura M. 2017. BRisk: *Bacillus cereus* group exposure assessment (beta).

## Launching BRisk Demo App (Latest Version)

1. Download R, if necessary: https://www.r-project.org/

2. Dowload R Studio, if necessary: https://www.rstudio.com/products/rstudio/download/

3. Open R Studio, and install the following packages, if necessary, by typing the following commands into R Studio's console:

```
install.packages("shiny")
install.packages("ggplot2")
install.packages("readr")
install.packages("stringr")
install.packages("fitdistrplus")
install.packages("mc2d")
```

4. Load shiny package

```
library(shiny)
```

5. Launch the app by typing the following command into R Studio's console:
```
runGitHub("BRisk","lmc297")
```

You're ready to go!

------------------------------------------------------------------------

# Uploading Files

1. Upload a BTyper final results file to BRisk by clicking the "Browse" button under "Choose a BTyper final results file to analyze" in the left panel (final results files are text files with the extension "_final_results.txt"). 
Note: several sample BTyper final results files are provided in BRisk's sample_data directory: anthrax_final_results.txt, diarrheal_final_results.txt, emetic_final_results.txt, and manliponensis_final_results.txt. Feel free to use any of those, or upload your own BTyper final results file!

2. Select a food matrix in the drop-down menu titled "Select a food matrix"

3. Uplaod a bacterial count data file (a text file that has *B. cereus* group isolate counts in log(CFU/g) or log(CFU/mL), with one count per line).
Note: several sample bacterial count data files are provided in BRisk's sample_data directory: counts1.txt, counts2.txt, counts3.txt, and counts4.txt. Feel free to use any of those, or upload your own count data files!

------------------------------------------------------------------------
