# Xicor.jl

An implementation of the Xi Correlation coefficient and hypothesis test as originally described by [Chatterjee (2021)](https://doi.org/10.1080/01621459.2020.1758115). The implementation is based on the XICOR CRAN package by [Chatterjee and Holmes](https://cran.r-project.org/web/packages/XICOR/index.html).

## Installation

Xicor.jl is available from the Julia registry, you can install it locally from your Julia REPL by typing:

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
Due to the asymmetry of the Xi correlation coefficient, it finds an interesting application in the inference of directionality of cause-effect relationships. To test its strengths in this endeavor, we apply the coefficient on the [Tuebingen cause-effect pairs database](https://webdav.tuebingen.mpg.de/cause-effect/) described in [Mooij _et al._ (2017)](http://jmlr.org/papers/v17/14-518.html). The data can be downloaded from 
https://webdav.tuebingen.mpg.de/cause-effect/pairs.zip, and consists of a
number of txt files, each containing a list of cause-effect pairs and 
description files containing information about the data sources and the 
ground truth direction of the causal relationships.

We test it here on the first pair file, which contains data from the German Weather Service on the relationship between temperature and altitude. The data is loaded in and the Xi correlation coefficient is calculated for both directions of the causal relationship. The direction with the higher Xi correlation coefficient is then chosen as the predicted direction of the causal relationship.

```julia
    # load in the data
    pair_values = hcat(
        [parse.(
            Float64, 
            [q for q in split(i, r"\s+") if q != ""]
            ) for i in split(chomp(readchomp("./data/pair0001.txt")), "\n")
        ]...
        )
    altitudes = pair_values[1, :]
    temperatures = pair_values[2, :]

    # calculate the Xi correlation coefficient for both directions
    using Xicor
    xy_xi, xy_sd, xy_p = xicor(altitudes, temperatures; pvalue=true)
    yx_xi, yx_sd, yx_p = xicor(temperatures, altitudes; pvalue=true)
    predicted_direction = xy_xi > yx_xi ? "Altitude causes temperature" : "Temperature causes altitude"
    println("Predicted direction: $predicted_direction")
```

And it correctly predicts the direction of the causal relationship in this case.

For a deeper exploration of directional inference using Xi correlation, check out the [causal direction prediction example](https://github.com/stefftaelman/Xicor.jl/tree/main/examples/causal_direction_prediction.jl).