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
library(mlbench)

# Setting simulation parameters
n <- 100
seed <- 42
sd <- 0.1

# Loading the data
fried_data <- as.data.frame(mlbench::mlbench.friedman1(n = n,sd = 0.1))

# Setting a cross-validation object
cv_obj <- k_fold(data = fried_data,dependent_variable = "y",
                 k_partitions = 5,seed = seed,as_data_frame = TRUE)

# Selecting one fold
fold <- 2
x_train <- cv_obj[[fold]]$x_train
y_train <- c(unlist(cv_obj[[fold]]$y_train))
x_test <- cv_obj[[fold]]$x_test
y_test <- cv_obj[[fold]]$y_test


```

### Running the model

To run the model we would have:

```{r, eval = FALSE}
gp_bart_mod <- gpbart(x_train = x_train,
                       y = c(unlist(y_train)),
                       x_test = x_test,
                       n_tree = 20)
```
