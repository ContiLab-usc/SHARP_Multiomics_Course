# R script to install all packages needed for workshop

# R packages
if (!requireNamespace("knitr", quietly = TRUE)) install.packages("knitr")
if (!requireNamespace("gesso", quietly = TRUE)) install.packages("gesso")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("CCA", quietly = TRUE)) install.packages("CCA")
if (!requireNamespace("corrplot", quietly = TRUE)) install.packages("corrplot")
if (!requireNamespace("factoextra", quietly = TRUE)) install.packages("factoextra")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
if (!requireNamespace("bama", quietly = TRUE)) install.packages("bama")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("networkD3", quietly = TRUE)) install.packages("networkD3")
if (!requireNamespace("mediation", quietly = TRUE)) install.packages("mediation")
if (!requireNamespace("sm", quietly = TRUE)) install.packages("sm")
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
if (!requireNamespace("PMA", quietly = TRUE)) install.packages("PMA")
if (!requireNamespace("gap", quietly = TRUE)) install.packages("gap")
if (!requireNamespace("kableExtra", quietly = TRUE)) install.packages("kableExtra")
if (!requireNamespace("table1", quietly = TRUE)) install.packages("table1")
if (!requireNamespace("lmtest", quietly = TRUE)) install.packages("lmtest")
if (!requireNamespace("stargazer", quietly = TRUE)) install.packages("stargazer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
if (!requireNamespace("gtools", quietly = TRUE)) install.packages("gtools")
if (!requireNamespace("interactionR", quietly = TRUE)) install.packages("interactionR")
if (!requireNamespace("epiR", quietly = TRUE)) install.packages("epiR")
if (!requireNamespace("jtools", quietly = TRUE)) install.packages("jtools")
if (!requireNamespace("interactions", quietly = TRUE)) install.packages("interactions")
if (!requireNamespace("msm", quietly = TRUE)) install.packages("msm")
if (!requireNamespace("glasso", quietly = TRUE)) install.packages("glasso")
if (!requireNamespace("lbfgs", quietly = TRUE)) install.packages("lbfgs")
if (!requireNamespace("lavaan", quietly = TRUE)) install.packages("lavaan")
if (!requireNamespace("RGCCA", quietly = TRUE)) install.packages("RGCCA")
if (!requireNamespace("LUCIDus", quietly = TRUE)) install.packages("LUCIDus")
if (!requireNamespace("rsvd", quietly = TRUE)) install.packages("rsvd")
if (!requireNamespace("pls", quietly = TRUE)) install.packages("pls")
if (!requireNamespace("summarytools", quietly = TRUE)) install.packages("summarytools")
if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc")
if (!requireNamespace("janitor", quietly = TRUE)) install.packages("janitor")
if (!requireNamespace("BinaryDosage", quietly = TRUE)) install.packages("BinaryDosage")
if (!requireNamespace("GxEScanR", quietly = TRUE)) install.packages("GxEScanR")


# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# This warning is expected and can be safely ignored:
# "Warning message: package(s) not installed when version(s) same as or 
# greater than current; use `force = TRUE` to re-install: 'MultiAssayExperiment'"
BiocManager::install("MultiAssayExperiment")
BiocManager::install("Biobase")
BiocManager::install("ComplexHeatmap")
BiocManager::install("qvalue")
BiocManager::install("mogsa")

install.packages("devtools")
# devtools install 
devtools::install_github("kupietz/kableExtra")
if (!requireNamespace("IntNMF", quietly = TRUE)) install.packages("IntNMF")


