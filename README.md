### GP-BART

## Installation

You can install the development version:

``` r
devtools::install_github("MateusMaiaDS/gpbart")
```


### Setting a example

This a vignette to explain how to run a simple example of the model, setting its own prior and its hyperparameters. To start we going to use the `friedman` example as the dataset to be used.

```{r setup, eval = FALSE}
library(gpbart)
# Setting simualation parameters
n <- 100
seed <- 42
sd <- 0.1
# Loading the data
fried_data <- sim_friedman(n = n,seed = seed,sd = sd)
# Setting a cross-validation object
cv_obj <- k_fold(data = fried_data,dependent_variable = "y",
                 k_partitions = 5,seed = seed)
# Selecting one fold
fold <- 1
x_train <- cv_obj[[fold]]$x_train
y_train <- cv_obj[[fold]]$y_train
x_test <- cv_obj[[fold]]$x_test
y_test <- cv_obj[[fold]]$y_test
```

### Running the model

To run the model we would have:

```{r, eval = FALSE}
gp_bart_mod <- gp_bart(x_train = x_train,
                       y_train = c(y_train),
                       x_test = x_test,
                       n_tree = 10,
                       gp_variables_ = colnames(x_train), # Selecting all var.
                       rotation_variables_ = colnames(x_train)) 
```