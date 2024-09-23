## MortalityBrazil

### Description 
`MortalityBrazil` provides functions and tools to create and visualize continuous
map of risk by combining spatially discrete mortality data with spatial covariates.

To achieve this aim, a `dissagregation model` based on the ELGM approach by Stringer and others 2021 (DOI: https://doi.org/10.1080/10618600.2022.2099403) is implemented. 
This package utilizes source code of the disaggregation* package:  https://cran.r-project.org/web/packages/disaggregation/index.html



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
- **Dissagregation Analysis**: Functions for conduction disaggregation analysis using ELGM.
- **Visualization**: Create plots and maps using `ggplot2` and `mapmisc` -based functions
