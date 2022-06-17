pvcm
================
Joonho Gong

``` r
library(pvcm)
data(yeast)
```

## Functions

### pvcm

`pvcm` is a function to fit a principal varying coefficient model with
adaptive weights. Parallel computting for multi-core CPU is built-in to
reduce computational time.

#### Usage

``` r
set.seed(12345)


## Preparation for design matrix X and the long y
# 1. y_vec
y_ready <- c( t( yeast$y ) )
n <- length(y_ready)

# 2. X_mat
p_org <- dim(yeast$x)[2]
X_temp <- apply( scale(yeast$x), 1, FUN = function(x){
  matrix( rep(x, 18), 18, p_org, byrow=TRUE)
}, simplify = FALSE)

X_mat <- matrix(NA, n ,p_org)
for(i in 1:length(X_temp)){
  X_mat[ ( 1 + 18 * (i- 1) ):(18 * i) , ] <- X_temp[[i]]
}

X_ready <- cbind( rep(1, n), X_mat )
colnames(X_ready) <- c( "intercept", colnames(yeast$x) )

## Fit PVCM
pvcm_fit <- pvcm(y_sample = y_ready, X_sample = X_ready, 
                 U_sample = rep(seq(0,1,length.out = 18), dim(yeast$x)[1]), 
                 knots_val = 7, 
                 control = list(num_tau_1=20, num_tau_2=20, max.iter=100, tol=1e-5)
                 )

true_TFs <- c("intercept",
              "ACE2_YPD", "SWI4_YPD", "SWI5_YPD", "SWI6_YPD", "MBP1_YPD", "STB1_YPD", "FKH1_YPD",
              "FKH2_YPD", "NDD1_YPD", "MCM1_YPD", "ABF1_YPD", "BAS1_YPD", "CBF1_YPD", "GCN4_YPD",
              "GCR1_YPD", "GCR2_YPD", "LEU3_YPD", "MET31_YPD", "REB1_YPD", "SKN7_YPD", "STE12_YPD")
# Experimentally confirmed covariates
sort(true_TFs)
#>  [1] "ABF1_YPD"  "ACE2_YPD"  "BAS1_YPD"  "CBF1_YPD"  "FKH1_YPD"  "FKH2_YPD" 
#>  [7] "GCN4_YPD"  "GCR1_YPD"  "GCR2_YPD"  "intercept" "LEU3_YPD"  "MBP1_YPD" 
#> [13] "MCM1_YPD"  "MET31_YPD" "NDD1_YPD"  "REB1_YPD"  "SKN7_YPD"  "STB1_YPD" 
#> [19] "STE12_YPD" "SWI4_YPD"  "SWI5_YPD"  "SWI6_YPD"

true_positive <-  sum(pvcm_fit$selected_covariates %in% true_TFs) - 1
# Not consider the intercept

true_negative <- p_org - (length(true_TFs) - 1) - sum(!(pvcm_fit$selected_covariates %in% true_TFs))
# except intercept

accuracy <- (true_positive + true_negative) / p_org
# 0.754717
```

#### Argument

**y\_sample:** A numeric vector of response variable.

**X\_sample:** A design matrix of covariates including an intercept if
it is included in your model.

**U\_sample:** A numeric vector of index variable.

**num\_cores:** The number of CPU cores to implement parallel computing.
Default NULL will use the total number of cores - 1).

**knots\_val:** A numeric value for the knots argument of select.knots
function in the face R package. Default is 7 for cubic B-spline.

**control:** A list that contains two numbers of tuning parameters
except 0 and the stopping criteria for DYS algorithm. There are 4
elements: **num\_tau\_1**, **num\_tau\_2**, **max.iter**, and **tol**.
First two are the number of candidates for each penalty except 0.
Default is 20 and zero will be automatically added. The total number of
iterations, **max.iter**, and numeric tolerance, **tol**, are stopping
criteria of DYS algorithm with default values 100 and 1e-5,
respectively.

#### Value

**Theta\_sparse:** The estimated Theta matrix with zero rows.

**Theta\_low\_rank:** The estimated Theta matrix of low rank.

**Theta:** The estimated Theta matrix from DYS algorithm close to
**Theta\_sparse:** and **Theta\_low\_rank**.

**rank:** The estimated rank of Theta matrix.

**U:** The matrix of left singular vectors of **Theta\_low\_rank**.

**V:** The matrix of right singular vectors of **Theta\_low\_rank**.

**D:** The diagonal matrix of singular values of **Theta\_low\_rank**.

**B\_U\_mat:** A matrix of orthogonal basis functions evaluated at
**U\_sample**.

**fitted\_values:** A vector of fitted values corresponding to
**y\_sample**.

**selected\_covariates:** The names of selected covariates by pvcm.

### dys\_pvcm

`dys_pvcm` is a function specifically developed for the principal
varying coefficient models.

#### Argument

**y\_vec:** A numeric vector of response variable.

**M\_mat:** A matrix whose row is the kroneker product between the
covariates X\_i and the evaluated basis functions at U\_i.

**tau:** A vector of the two tuning parameters.

**Theta\_initial:** An initial value for the Theta matrix.

**weights\_ann:** A vector of weights for the adaptive nuclear norm
penalty.

**weights\_agl:** A vector of weights for the adaptive group lasso
penalty.

**gamma\_val:** A numeric value for one over Lipschitz constant.

**max.iter:** A maximum iteration number.

**tol:** A numeric value for tolerance of relative deviation in order to
stop DYS iteration.

#### Value

**Theta\_sparse:** The estimated Theta matrix with zero rows.

**Theta\_low\_rank:** The estimated Theta matrix of low rank.

**Theta:** The estimated Theta matrix from DYS algorithm close to
**Theta\_sparse:** and **Theta\_low\_rank**.

**rank:** The estimated rank of Theta matrix.

**U:** The matrix of left singular vectors of **Theta\_low\_rank**.

**V:** The matrix of right singular vectors of **Theta\_low\_rank**.

**D:** The diagonal matrix of singular values of **Theta\_low\_rank**.

**converge:** A logical value to indicate convergence.

**rank\_trace:** A vector of estimated rank of Theta matrix over
iteration.

**iteration:** The total number of iterations

**rel\_diff:** A vector of relative differences over iteration.

### plot\_pvcm

`plot_pvcm` is a function to display estimated varying coefficient
functions and pricipal functions.

#### Usage

``` r

par(mar = c(3.8, 2, 1.8, 1.8) )
## Varying coefficient functions
plot_pvcm(pvcm_fit, pallet_size = c(4,4), X_sample = X_ready,
          U_sample = rep(seq(0,1,length.out = 18), dim(yeast$x)[1]),
          beta_plot = TRUE, true_covariates = true_TFs
          )

## Principal functions
plot_pvcm(pvcm_fit, pallet_size = c(2,2), 
          U_sample = rep(seq(0,1,length.out = 18), dim(yeast$x)[1]),
          beta_plot = FALSE
          )
```

![](markdown_example_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->![](markdown_example_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

#### Argument

**pvcm\_fit:** A list output of pvcm function.

**pallet\_size:** A vector of two positive integer elements for the
number of rows and columns of pallet.

**X\_sample:** A design matrix of covariates with an intercept if it is
included in your model.

**U\_sample:** A numeric vector of index variable.

**beta\_plot:** A logical value to indicate whether the coefficient
functions (TRUE) or the principal functions (FALSE) are displayed.

**true\_covariates:** A vector of names for the known true covariates.
Default NULL will display the plots for all selected covariates when
**beta\_plot** is TRUE.
