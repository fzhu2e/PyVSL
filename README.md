# PyVSL
VS-Lite in Python.

It contains two modes:
1. A python wrapper of the original R package of VS-Lite by Suz Tolwinski-Ward by Suz Tolwinski-Ward (https://github.com/suztolwinskiward/VSLiteR)
2. VS-Lite in pure Python (TODO)

## How to install

For the wrapper mode, we need to install the original R package of VS-Lite, which can be done by executing the below lines in the R console:
```R
install.packages("devtools")
library(devtools)
install_github("fzhu2e/VSLiteR")
```

After that, we need to install the `rpy2` python package:
```bash
pip install rpy2
```

Then we are ready to install `PyVSL`:
```bash
pip install PyVSL
```
to install from PyPi, or
```bash
pip install -e .
```
in the directory of the downloaded repo to install from the local path.