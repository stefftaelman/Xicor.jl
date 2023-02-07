# Xicor.jl

An implementation of the Xi Correlation coefficient and hypothesis test as originally described by [Chatterjee (2021)](https://doi.org/10.1080/01621459.2020.1758115). The implementation is based on the XICOR CRAN package by [Chatterjee and Holmes](https://cran.r-project.org/web/packages/XICOR/index.html).

## Installation

Once the package gets added to the Julia registry, you can install it locally by

```julia
] add Xicor
```

## Running

Xicor.jl allows to eaily calculate the Xi correlation coefficient between two vectors $X$ and $Y$. 

### Calculating Xi correlation

```julia
xvec = rand(100)
yvec = rand(100)
xicor(xvec, yvec)
```

### Dependence testing
Aside retrieving a correlation coefficient, Xicor.jl also allows to test for dependence between two vectors $X$ and $Y$. This is possible by using the `xicor` fucntion with the keyword argument 'pvalue' set to `true`. The function will then return a tuple of the Xi correlation coefficient, the standard deviation, and the associated p-value.

```julia
xicor(xvec, yvec; pvalue=true, method="asymptotic")
```

The `method` keyword argument allows to choose between an asymptotic approximation test and a Monte Carlo permutation 
test. The asymptotic approximation test is based on the asymptotic normal distribution of the Xi correlation coefficient. The Monte Carlo permutation test is based on the permutation distribution of the Xi correlation coefficient. The Monte Carlo permutation test also uses a keyword argument `nperm` to specify the number of permutations to be used. The default value is 1000.

```julia
xicor(xvec, yvec; pvalue=true, method="permutation", nperm=1000)
```

## Application