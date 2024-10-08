## MortalityBrazil

### Description 
`MortalityBrazil` provides functions and tools to create and visualize continuous
map of risk by combining spatially discrete mortality data with spatial covariates.

To achieve this aim, a dissagregation model based on the ELGM approach by Stringer and others 2021 (DOI: https://doi.org/10.1080/10618600.2022.2099403) is implemented. 
This package utilizes source code of the `disaggregation` package:  https://cran.r-project.org/web/packages/disaggregation/index.html



### **Table of Contents**
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Features](#features)

### **Installation**

```r
devtools::install_github("sofiar/MortalityBrazil")
```

### **Quick Start**
```r
library(MortalityBrazil)

# Load Brazil limits
data(brazil_limits)

# Load population offset
population_offset = load_population_offset()

# Mortality data
data(mortality_data)

```

### **Features**
- **Demographic Loading**: Easily load various demographic spatial covarites of Brazil based on the census 2010.
- **Dissagregation Analysis**: Functions for conducting disaggregation analysis using ELGM.
- **Visualization**: Create plots and maps using `ggplot2` and `mapmisc` - based functions



### **References**
Stringer, A., Brown, P., & Stafford, J. (2022). Fast, Scalable
Approximations to Posterior Distributions in Extended Latent Gaussian Models.
Journal of Computational and Graphical Statistics, 32(1), 84–98. https://doi.org/10.1080/10618600.2022.2099403

