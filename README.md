# Developmental modularity drives convergence and diversity in Cetacea 🐬🐋🐟🧠🦴📈
### Testing the influence of different developemtnal patterns between the rostrum and braincase in baleen and toothed whales using 3D geometric morphometrics 

Authors: [Agnese Lanzetti](mailto:agnese.lanzetti@gmail.com?subject=[GitHub]%20Modularity%20Cetacea%20Paper%20Code),
Anjali Goswami

To cite the paper: 

Available at: https://github.com/AgneseLan/cetacea-mouth-brain

If using any of this code or data please cite the paper above and this repo.
To cite this repo: 


## Data :floppy_disk: 

The data are provided in the Data folder. Mesh files (PLY) needed to test postioning when importing landmarks are available at Phenome10k (https://www.phenome10k.org/).

- __Landmark data__: *pts folder* <br />
Text files with landmark coordinates for each specimen in PTS format. Unzip folder first.

- __Surface data__: *ply folder* <br />
Empty folder where mesh files from Phenome10k need to be saved to reproduce the code.

- __Specimens' classifiers, landmark/curves lists__: *absent_curves_all.csv, absent_LMs_all.csv, curves_all.csv, LMs_all.csv, info_lumbar.csv, specimens_all.csv* <br />
Spreadsheets with additional inforation for analyses: list of absent landmarks and curves, list of curves, list of landmarks, classifiers for specimens (specimen names, ID, group, family, genera, common name, length, bizygomatic width, age, growth stage).

- __Reference mesh for plotting__: *refmesh_all.zip* (*myst_adult.ply, myst_fetus.ply, odont_adult.ply, odont_fetus.ply, refmesh_all.ply*) <br />
Reduced meshes in PLY format used for plotting landmarks.

- __Silhouettes of taxa for plots__: *megaptera.png, stenella.png*

## Analysis :computer:
In this repository you will find raw data (.csv and data files) and code for analyses (code supplied as .R files)

📁 Data

As described above. Meshes used to collect and test landmarks available on Phenome10k (https://www.phenome10k.org/). 

⌨ Code for analyses - .R files

*1-Import-resample-slide.R, 1-Slider3d_2.R, 2-Absent_bones.R, 3-GPA.R, 4-Modularity.R, 5-PCA.R, 6-Disparity.R, 7-Trajectory.R, 8-Allometry.R*

Code files are numbered providing the order the analyses need to be performed in.
Before running analyses, save Data folder in the same directory as the R project. This will allow to import the data as detailed in the code provided. Also create an Output folder and one subfolder for each script to allow for the results to export correctly.

## License 📃
This project is licensed under the MIT License - see the LICENSE.md file for details

## Session Info 📋
For reproducibility purposes, here is the output of sessioninfo::session_info() used to perform the analyses in the publication.

```
─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.3.2 (2023-10-31 ucrt)
 os       Windows 10 x64 (build 19045)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  English_United Kingdom.utf8
 ctype    English_United Kingdom.utf8
 tz       Europe/London
 date     2024-01-23
 rstudio  2023.12.0+369 Ocean Storm (desktop)
 pandoc   3.1.1 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package           * version    date (UTC) lib source
 abind             * 1.4-5      2016-07-21 [1] CRAN (R 4.3.1)
 ape               * 5.7-1      2023-03-13 [1] CRAN (R 4.3.2)
 backports           1.4.1      2021-12-13 [1] CRAN (R 4.3.1)
 base64enc           0.1-3      2015-07-28 [1] CRAN (R 4.3.1)
 bezier              1.1.2      2018-12-14 [1] CRAN (R 4.3.1)
 borealis          * 2022.10.27 2024-01-12 [1] Github (aphanotus/borealis@b2d9fc1)
 broom               1.0.5      2023-06-09 [1] CRAN (R 4.3.2)
 cachem              1.0.8      2023-05-01 [1] CRAN (R 4.3.2)
 car               * 3.1-2      2023-03-30 [1] CRAN (R 4.3.2)
 carData           * 3.0-5      2022-01-06 [1] CRAN (R 4.3.2)
 checkmate           2.3.1      2023-12-04 [1] CRAN (R 4.3.2)
 cli                 3.6.2      2023-12-11 [1] CRAN (R 4.3.2)
 cluster             2.1.6      2023-12-01 [2] CRAN (R 4.3.2)
 clusterGeneration   1.3.8      2023-08-16 [1] CRAN (R 4.3.2)
 coda                0.19-4     2020-09-30 [1] CRAN (R 4.3.2)
 codetools           0.2-19     2023-02-01 [2] CRAN (R 4.3.2)
 colorRamps          2.3.1      2022-05-02 [1] CRAN (R 4.3.1)
 colorspace          2.1-0      2023-01-23 [1] CRAN (R 4.3.2)
 combinat            0.0-8      2012-10-29 [1] CRAN (R 4.3.1)
 corpcor             1.6.10     2021-09-16 [1] CRAN (R 4.3.1)
 cowplot             1.1.2      2023-12-15 [1] CRAN (R 4.3.2)
 curl                5.2.0      2023-12-08 [1] CRAN (R 4.3.2)
 data.table          1.14.10    2023-12-08 [1] CRAN (R 4.3.2)
 deSolve             1.40       2023-11-27 [1] CRAN (R 4.3.2)
 devtools          * 2.4.5      2022-10-11 [1] CRAN (R 4.3.2)
 digest              0.6.33     2023-07-07 [1] CRAN (R 4.3.2)
 doParallel          1.0.17     2022-02-07 [1] CRAN (R 4.3.2)
 dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.3.2)
 ellipsis            0.3.2      2021-04-29 [1] CRAN (R 4.3.2)
 EMMLi             * 0.0.3      2017-02-17 [1] CRAN (R 4.3.2)
 evaluate            0.23       2023-11-01 [1] CRAN (R 4.3.2)
 expm                0.999-8    2023-11-29 [1] CRAN (R 4.3.2)
 fansi               1.0.6      2023-12-08 [1] CRAN (R 4.3.2)
 fastmap             1.1.1      2023-02-24 [1] CRAN (R 4.3.2)
 fastmatch           1.1-4      2023-08-18 [1] CRAN (R 4.3.1)
 fdrtool             1.2.17     2021-11-13 [1] CRAN (R 4.3.1)
 forcats           * 1.0.0      2023-01-29 [1] CRAN (R 4.3.2)
 foreach             1.5.2      2022-02-02 [1] CRAN (R 4.3.2)
 foreign             0.8-86     2023-11-28 [2] CRAN (R 4.3.2)
 Formula             1.2-5      2023-02-24 [1] CRAN (R 4.3.1)
 fs                  1.6.3      2023-07-20 [1] CRAN (R 4.3.2)
 geiger            * 2.0.11     2023-04-03 [1] CRAN (R 4.3.2)
 generics            0.1.3      2022-07-05 [1] CRAN (R 4.3.2)
 geomorph          * 4.0.6      2023-08-31 [1] CRAN (R 4.3.2)
 ggfortify         * 0.4.16     2023-03-20 [1] CRAN (R 4.3.2)
 gginnards         * 0.1.2      2023-05-24 [1] CRAN (R 4.3.2)
 ggphylomorpho     * 0.1.0      2024-01-12 [1] Github (wabarr/ggphylomorpho@7e1228b)
 ggplot2           * 3.4.4      2023-10-12 [1] CRAN (R 4.3.2)
 ggplotify         * 0.1.2      2023-08-09 [1] CRAN (R 4.3.2)
 ggpubr            * 0.6.0      2023-02-10 [1] CRAN (R 4.3.2)
 ggrepel           * 0.9.5      2024-01-10 [1] CRAN (R 4.3.2)
 ggsignif            0.6.4      2022-10-13 [1] CRAN (R 4.3.2)
 ggthemes          * 5.0.0      2023-11-21 [1] CRAN (R 4.3.2)
 glasso              1.11       2019-10-01 [1] CRAN (R 4.3.1)
 glue                1.7.0      2024-01-09 [1] CRAN (R 4.3.2)
 gridExtra         * 2.3        2017-09-09 [1] CRAN (R 4.3.2)
 gridGraphics        0.5-1      2020-12-13 [1] CRAN (R 4.3.2)
 grImport2           0.3-1      2023-10-27 [1] CRAN (R 4.3.2)
 gtable              0.3.4      2023-08-21 [1] CRAN (R 4.3.2)
 gtools              3.9.5      2023-11-20 [1] CRAN (R 4.3.2)
 Hmisc               5.1-1      2023-09-12 [1] CRAN (R 4.3.2)
 hms                 1.1.3      2023-03-21 [1] CRAN (R 4.3.2)
 htmlTable           2.4.2      2023-10-29 [1] CRAN (R 4.3.2)
 htmltools           0.5.7      2023-11-03 [1] CRAN (R 4.3.2)
 htmlwidgets         1.6.4      2023-12-06 [1] CRAN (R 4.3.2)
 httpuv              1.6.13     2023-12-06 [1] CRAN (R 4.3.2)
 httr                1.4.7      2023-08-15 [1] CRAN (R 4.3.2)
 igraph              1.6.0      2023-12-11 [1] CRAN (R 4.3.2)
 iterators           1.0.14     2022-02-05 [1] CRAN (R 4.3.2)
 janeaustenr         1.0.0      2022-08-26 [1] CRAN (R 4.3.2)
 jpeg                0.1-10     2022-11-29 [1] CRAN (R 4.3.1)
 jsonlite            1.8.8      2023-12-04 [1] CRAN (R 4.3.2)
 knitr               1.45       2023-10-30 [1] CRAN (R 4.3.2)
 later               1.3.2      2023-12-06 [1] CRAN (R 4.3.2)
 lattice             0.22-5     2023-10-24 [2] CRAN (R 4.3.2)
 lavaan              0.6-17     2023-12-20 [1] CRAN (R 4.3.2)
 lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.3.2)
 loo                 2.6.0      2023-03-31 [1] CRAN (R 4.3.2)
 lubridate         * 1.9.3      2023-09-27 [1] CRAN (R 4.3.2)
 magick            * 2.8.2      2023-12-20 [1] CRAN (R 4.3.2)
 magrittr            2.0.3      2022-03-30 [1] CRAN (R 4.3.2)
 maps              * 3.4.2      2023-12-15 [1] CRAN (R 4.3.2)
 MASS                7.3-60     2023-05-04 [2] CRAN (R 4.3.2)
 Matrix            * 1.6-4      2023-11-30 [2] CRAN (R 4.3.2)
 matrixStats         1.2.0      2023-12-11 [1] CRAN (R 4.3.2)
 mcp               * 0.3.3      2023-03-22 [1] CRAN (R 4.3.2)
 memoise             2.0.1      2021-11-26 [1] CRAN (R 4.3.2)
 mime                0.12       2021-09-28 [1] CRAN (R 4.3.1)
 miniUI              0.1.1.1    2018-05-18 [1] CRAN (R 4.3.2)
 mnormt              2.1.1      2022-09-26 [1] CRAN (R 4.3.1)
 Morpho            * 2.12       2023-12-06 [1] CRAN (R 4.3.2)
 munsell             0.5.0      2018-06-12 [1] CRAN (R 4.3.2)
 mvtnorm             1.2-4      2023-11-27 [1] CRAN (R 4.3.2)
 nlme                3.1-164    2023-11-27 [2] CRAN (R 4.3.2)
 nnet                7.3-19     2023-05-03 [2] CRAN (R 4.3.2)
 numDeriv            2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.1)
 optimParallel       1.0-2      2021-02-11 [1] CRAN (R 4.3.2)
 paleomorph        * 0.1.4      2017-04-19 [1] CRAN (R 4.3.2)
 patchwork           1.2.0      2024-01-08 [1] CRAN (R 4.3.2)
 pbapply             1.7-2      2023-06-27 [1] CRAN (R 4.3.2)
 pbivnorm            0.6.0      2015-01-23 [1] CRAN (R 4.3.1)
 phangorn            2.11.1     2023-01-23 [1] CRAN (R 4.3.2)
 phytools          * 2.1-1      2024-01-09 [1] CRAN (R 4.3.2)
 pillar              1.9.0      2023-03-22 [1] CRAN (R 4.3.2)
 pkgbuild            1.4.3      2023-12-10 [1] CRAN (R 4.3.2)
 pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.3.2)
 pkgload             1.3.3      2023-09-22 [1] CRAN (R 4.3.2)
 plyr                1.8.9      2023-10-02 [1] CRAN (R 4.3.2)
 png               * 0.1-8      2022-11-29 [1] CRAN (R 4.3.1)
 profvis             0.3.8      2023-05-02 [1] CRAN (R 4.3.2)
 promises            1.2.1      2023-08-10 [1] CRAN (R 4.3.2)
 psych               2.3.12     2023-12-20 [1] CRAN (R 4.3.2)
 purrr             * 1.0.2      2023-08-10 [1] CRAN (R 4.3.2)
 qgraph            * 1.9.8      2023-11-03 [1] CRAN (R 4.3.2)
 quadprog            1.5-8      2019-11-20 [1] CRAN (R 4.3.1)
 R6                  2.5.1      2021-08-19 [1] CRAN (R 4.3.2)
 RColorBrewer      * 1.1-3      2022-04-03 [1] CRAN (R 4.3.1)
 Rcpp                1.0.12     2024-01-09 [1] CRAN (R 4.3.2)
 readr             * 2.1.5      2024-01-10 [1] CRAN (R 4.3.2)
 remotes             2.4.2.1    2023-07-18 [1] CRAN (R 4.3.2)
 reshape2          * 1.4.4      2020-04-09 [1] CRAN (R 4.3.2)
 rgl               * 1.2.8      2023-11-29 [1] CRAN (R 4.3.2)
 rlang               1.1.3      2024-01-10 [1] CRAN (R 4.3.2)
 rmarkdown           2.25       2023-09-18 [1] CRAN (R 4.3.2)
 rpart               4.1.23     2023-12-05 [2] CRAN (R 4.3.2)
 rphylopic         * 1.3.0      2023-12-20 [1] CRAN (R 4.3.2)
 RRPP              * 1.4.0      2023-08-15 [1] CRAN (R 4.3.2)
 rstatix             0.7.2      2023-02-01 [1] CRAN (R 4.3.2)
 rstudioapi          0.15.0     2023-07-07 [1] CRAN (R 4.3.2)
 rsvg                2.6.0      2023-10-08 [1] CRAN (R 4.3.2)
 Rvcg              * 0.22.2     2023-12-06 [1] CRAN (R 4.3.2)
 scales            * 1.3.0      2023-11-28 [1] CRAN (R 4.3.2)
 scatterplot3d       0.3-44     2023-05-05 [1] CRAN (R 4.3.1)
 sessioninfo         1.2.2      2021-12-06 [1] CRAN (R 4.3.2)
 shiny               1.8.0      2023-11-17 [1] CRAN (R 4.3.2)
 SnowballC           0.7.1      2023-04-25 [1] CRAN (R 4.3.1)
 stringi             1.8.3      2023-12-11 [1] CRAN (R 4.3.2)
 stringr           * 1.5.1      2023-11-14 [1] CRAN (R 4.3.2)
 subplex             1.8        2022-04-12 [1] CRAN (R 4.3.1)
 SURGE             * 0.1.0      2024-01-12 [1] Github (rnfelice/SURGE@225a842)
 tibble            * 3.2.1      2023-03-20 [1] CRAN (R 4.3.2)
 tidyr             * 1.3.0      2023-01-24 [1] CRAN (R 4.3.2)
 tidyselect          1.2.0      2022-10-10 [1] CRAN (R 4.3.2)
 tidytext          * 0.4.1      2023-01-07 [1] CRAN (R 4.3.2)
 tidyverse         * 2.0.0      2023-02-22 [1] CRAN (R 4.3.2)
 timechange          0.2.0      2023-01-11 [1] CRAN (R 4.3.2)
 tokenizers          0.3.0      2022-12-22 [1] CRAN (R 4.3.2)
 tzdb                0.4.0      2023-05-12 [1] CRAN (R 4.3.2)
 urlchecker          1.0.1      2021-11-30 [1] CRAN (R 4.3.2)
 usethis           * 2.2.2      2023-07-06 [1] CRAN (R 4.3.2)
 utf8                1.2.4      2023-10-22 [1] CRAN (R 4.3.2)
 vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.3.2)
 withr               2.5.2      2023-10-30 [1] CRAN (R 4.3.2)
 xfun                0.41       2023-11-01 [1] CRAN (R 4.3.2)
 XML                 3.99-0.16  2023-11-29 [1] CRAN (R 4.3.2)
 xtable              1.8-4      2019-04-21 [1] CRAN (R 4.3.2)
 yulab.utils         0.1.3      2024-01-08 [1] CRAN (R 4.3.2)

 [1] C:/Users/agnL/AppData/Local/R/win-library/4.3
 [2] C:/Program Files/R/R-4.3.2/library

───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
