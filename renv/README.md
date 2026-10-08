# renv/

`msna_analysis_main` does not use renv (packages are loaded with `pacman` from
the user library), so renv is **not activated** in this project: activating it
would create a private library and re-install ~80 packages that are already
installed.

`../renv.lock` records the exact package versions this project was developed,
validated and benchmarked with (data.table, srvyr, survey, analysistools,
writexl, readxl, testthat, ... and their dependencies), taken from the user
library on 2026-10-08.

To switch to an isolated renv library later (e.g. on another computer):

```r
setwd("C:/Users/User/Music/new_msna_analysis")
renv::init(bare = TRUE)   # creates renv/activate.R and .Rprofile
renv::restore()           # installs the versions recorded in renv.lock
```
