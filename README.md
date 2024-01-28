# Tropical Logistic Regression

1. Coalescent model
First generate data with python's dendropy library (this may take a while)
```
cd coalescent_model
python coalescent_data_generation.py
python gen_err_dataset.py
```
Then, run the following in R 
```
source("./coalescent_model/coalescent_model_data.R")
```
To explore the distance distributions in Euclidean, tropical and BHV geometry
```
source("./coalescent_model/covariate_distribution/distance_distribution.R") 
```
To compare classical, tropical and BHV logistic regression, run the following
```
source("./coalescent_model/covariate_distribution/gen_error.R")
```
2. Mr Bayes convergence criterion
Run the following in the terminal (warning; it takes a long time to run)
```
./mr_bayes/mb.sh mr_bayes/primates.nex mb_data  100000 1 200 100
./mr_bayes/tlr.sh mb_data/ mr_bayes/primates_iterations.txt 
python3 mr_bayes/plotting.py mb_data mr_bayes/primates_iterations.txt 200
```
