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

- __Phylogenies by stage for mean shapes, all data and by group__: *tree_all_adult.txt*, *tree_all_adult_myst.txt*, *tree_all_adult_odont.txt*, *tree_all_early.txt*, *tree_all_early_myst.txt*, *tree_all_early_odont.txt*, *tree_all_immature.txt*, *tree_all_immature_myst.txt*, 
                                                                   *tree_all_immature_odont.txt*, *tree_all_late_new.txt*, *tree_all_late_new_myst.txt*, *tree_all_late_new_odont.txt*

## Analysis :computer:
In this repository you will find raw data (.csv and data files) and code for analyses (code supplied as .R files)

📁 Data

As described above. Meshes used to collect and test landmarks available on Phenome10k (https://www.phenome10k.org/). 

⌨ Code for analyses - .R files

*1-Import-resample-slide.R, 1-Slider3d_2.R, 2-Absent_bones.R, 3-GPA mean shapes.R, 4-Modularity means.R, 5-PCA means.R, 6-Disparity means.R, 7-Trajectory means.R, 8-Allometry means.R*

Code files are numbered providing the order the analyses need to be performed in.
Before running analyses, save Data folder in the same directory as the R project. This will allow to import the data as detailed in the code provided. Also create an Output folder and one subfolder for each script to allow for the results to export correctly.

## License 📃
This project is licensed under the MIT License - see the LICENSE.md file for details

## Session Info 📋
For reproducibility purposes, here is the output of sessioninfo::session_info() used to perform the analyses in the publication.

```
─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.4.1 (2024-06-14 ucrt)
 os       Windows 10 x64 (build 19045)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  English_United Kingdom.utf8
 ctype    English_United Kingdom.utf8
 tz       Europe/London
 date     2024-08-09
 rstudio  2024.04.2+764 Chocolate Cosmos (desktop)
 pandoc   3.1.11 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package           * version    date (UTC) lib source
 abind             * 1.4-5      2016-07-21 [1] CRAN (R 4.4.0)
 ape               * 5.8        2024-04-11 [1] CRAN (R 4.4.1)
 backports           1.5.0      2024-05-23 [1] CRAN (R 4.4.0)
 base64enc           0.1-3      2015-07-28 [1] CRAN (R 4.4.0)
 bezier              1.1.2      2018-12-14 [1] CRAN (R 4.4.0)
 borealis          * 2022.10.27 2024-07-12 [1] Github (aphanotus/borealis@b2d9fc1)
 broom               1.0.6      2024-05-17 [1] CRAN (R 4.4.1)
 cachem              1.1.0      2024-05-16 [1] CRAN (R 4.4.1)
 car               * 3.1-2      2023-03-30 [1] CRAN (R 4.4.1)
 carData           * 3.0-5      2022-01-06 [1] CRAN (R 4.4.1)
 checkmate           2.3.1      2023-12-04 [1] CRAN (R 4.4.1)
 cli                 3.6.3      2024-06-21 [2] CRAN (R 4.4.1)
 cluster             2.1.6      2023-12-01 [2] CRAN (R 4.4.1)
 clusterGeneration   1.3.8      2023-08-16 [1] CRAN (R 4.4.1)
 coda                0.19-4.1   2024-01-31 [1] CRAN (R 4.4.1)
 codetools           0.2-20     2024-03-31 [2] CRAN (R 4.4.1)
 colorRamps          2.3.4      2024-03-07 [1] CRAN (R 4.4.0)
 colorspace          2.1-0      2023-01-23 [1] CRAN (R 4.4.1)
 combinat            0.0-8      2012-10-29 [1] CRAN (R 4.4.0)
 corpcor             1.6.10     2021-09-16 [1] CRAN (R 4.4.0)
 cowplot             1.1.3      2024-01-22 [1] CRAN (R 4.4.1)
 curl                5.2.1      2024-03-01 [1] CRAN (R 4.4.1)
 data.table          1.15.4     2024-03-30 [1] CRAN (R 4.4.1)
 DEoptim             2.2-8      2022-11-11 [1] CRAN (R 4.4.1)
 deSolve             1.40       2023-11-27 [1] CRAN (R 4.4.1)
 devtools          * 2.4.5      2022-10-11 [1] CRAN (R 4.4.1)
 digest              0.6.36     2024-06-23 [1] CRAN (R 4.4.1)
 doParallel          1.0.17     2022-02-07 [1] CRAN (R 4.4.1)
 dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.4.1)
 ellipsis            0.3.2      2021-04-29 [1] CRAN (R 4.4.1)
 emmeans           * 1.10.3     2024-07-01 [1] CRAN (R 4.4.1)
 EMMLi             * 0.0.3      2017-02-17 [1] CRAN (R 4.4.1)
 estimability        1.5.1      2024-05-12 [1] CRAN (R 4.4.1)
 evaluate            0.24.0     2024-06-10 [1] CRAN (R 4.4.1)
 expm                0.999-9    2024-01-11 [1] CRAN (R 4.4.1)
 fansi               1.0.6      2023-12-08 [1] CRAN (R 4.4.1)
 fastmap             1.2.0      2024-05-15 [1] CRAN (R 4.4.1)
 fastmatch           1.1-4      2023-08-18 [1] CRAN (R 4.4.0)
 fdrtool             1.2.17     2021-11-13 [1] CRAN (R 4.4.0)
 forcats           * 1.0.0      2023-01-29 [1] CRAN (R 4.4.1)
 foreach             1.5.2      2022-02-02 [1] CRAN (R 4.4.1)
 foreign             0.8-87     2024-06-26 [2] CRAN (R 4.4.1)
 Formula             1.2-5      2023-02-24 [1] CRAN (R 4.4.0)
 fs                  1.6.4      2024-04-25 [1] CRAN (R 4.4.1)
 geiger            * 2.0.11     2023-04-03 [1] CRAN (R 4.4.1)
 generics            0.1.3      2022-07-05 [1] CRAN (R 4.4.1)
 geomorph          * 4.0.8      2024-07-12 [1] Github (geomorphR/geomorph@374c1b5)
 ggfortify         * 0.4.17     2024-04-17 [1] CRAN (R 4.4.1)
 gginnards         * 0.2.0      2024-05-01 [1] CRAN (R 4.4.1)
 ggphylomorpho     * 0.1.0      2024-07-12 [1] Github (wabarr/ggphylomorpho@7e1228b)
 ggplot2           * 3.5.1      2024-04-23 [1] CRAN (R 4.4.1)
 ggplotify         * 0.1.2      2023-08-09 [1] CRAN (R 4.4.1)
 ggpubr            * 0.6.0      2023-02-10 [1] CRAN (R 4.4.1)
 ggrepel           * 0.9.5      2024-01-10 [1] CRAN (R 4.4.1)
 ggsignif            0.6.4      2022-10-13 [1] CRAN (R 4.4.1)
 ggthemes          * 5.1.0      2024-02-10 [1] CRAN (R 4.4.1)
 glasso              1.11       2019-10-01 [1] CRAN (R 4.4.0)
 glue                1.7.0      2024-01-09 [2] CRAN (R 4.4.1)
 gridExtra         * 2.3        2017-09-09 [1] CRAN (R 4.4.1)
 gridGraphics        0.5-1      2020-12-13 [1] CRAN (R 4.4.1)
 grImport2           0.3-1      2023-10-27 [1] CRAN (R 4.4.1)
 gtable              0.3.5      2024-04-22 [1] CRAN (R 4.4.1)
 gtools              3.9.5      2023-11-20 [1] CRAN (R 4.4.1)
 Hmisc               5.1-3      2024-05-28 [1] CRAN (R 4.4.1)
 hms                 1.1.3      2023-03-21 [1] CRAN (R 4.4.1)
 htmlTable           2.4.2      2023-10-29 [1] CRAN (R 4.4.1)
 htmltools           0.5.8.1    2024-04-04 [1] CRAN (R 4.4.1)
 htmlwidgets         1.6.4      2023-12-06 [1] CRAN (R 4.4.1)
 httpuv              1.6.15     2024-03-26 [1] CRAN (R 4.4.1)
 httr                1.4.7      2023-08-15 [1] CRAN (R 4.4.1)
 igraph              2.0.3      2024-03-13 [1] CRAN (R 4.4.1)
 iterators           1.0.14     2022-02-05 [1] CRAN (R 4.4.1)
 janeaustenr         1.0.0      2022-08-26 [1] CRAN (R 4.4.1)
 jpeg                0.1-10     2022-11-29 [1] CRAN (R 4.4.0)
 jsonlite            1.8.8      2023-12-04 [1] CRAN (R 4.4.1)
 knitr               1.48       2024-07-07 [1] CRAN (R 4.4.1)
 later               1.3.2      2023-12-06 [1] CRAN (R 4.4.1)
 lattice             0.22-6     2024-03-20 [2] CRAN (R 4.4.1)
 lavaan              0.6-18     2024-06-07 [1] CRAN (R 4.4.1)
 lifecycle           1.0.4      2023-11-07 [2] CRAN (R 4.4.1)
 loo                 2.8.0      2024-07-03 [1] CRAN (R 4.4.1)
 lubridate         * 1.9.3      2023-09-27 [1] CRAN (R 4.4.1)
 magick            * 2.8.4      2024-07-14 [1] CRAN (R 4.4.1)
 magrittr            2.0.3      2022-03-30 [2] CRAN (R 4.4.1)
 maps              * 3.4.2      2023-12-15 [1] CRAN (R 4.4.1)
 MASS                7.3-61     2024-06-13 [2] CRAN (R 4.4.1)
 Matrix            * 1.7-0      2024-04-26 [2] CRAN (R 4.4.1)
 matrixStats         1.3.0      2024-04-11 [1] CRAN (R 4.4.1)
 mcp               * 0.3.4      2024-03-17 [1] CRAN (R 4.4.1)
 memoise             2.0.1      2021-11-26 [1] CRAN (R 4.4.1)
 mime                0.12       2021-09-28 [1] CRAN (R 4.4.0)
 miniUI              0.1.1.1    2018-05-18 [1] CRAN (R 4.4.1)
 mnormt              2.1.1      2022-09-26 [1] CRAN (R 4.4.0)
 Morpho            * 2.12       2023-12-06 [1] CRAN (R 4.4.1)
 munsell             0.5.1      2024-04-01 [1] CRAN (R 4.4.1)
 mvtnorm             1.2-5      2024-05-21 [1] CRAN (R 4.4.1)
 nlme                3.1-165    2024-06-06 [2] CRAN (R 4.4.1)
 nnet                7.3-19     2023-05-03 [2] CRAN (R 4.4.1)
 numDeriv            2016.8-1.1 2019-06-06 [1] CRAN (R 4.4.0)
 optimParallel       1.0-2      2021-02-11 [1] CRAN (R 4.4.1)
 paleomorph        * 0.1.4      2017-04-19 [1] CRAN (R 4.4.1)
 patchwork           1.2.0      2024-01-08 [1] CRAN (R 4.4.1)
 pbapply             1.7-2      2023-06-27 [1] CRAN (R 4.4.1)
 pbivnorm            0.6.0      2015-01-23 [1] CRAN (R 4.4.0)
 phangorn            2.11.1     2023-01-23 [1] CRAN (R 4.4.1)
 phytools          * 2.3-0      2024-06-13 [1] CRAN (R 4.4.1)
 pillar              1.9.0      2023-03-22 [1] CRAN (R 4.4.1)
 pkgbuild            1.4.4      2024-03-17 [1] CRAN (R 4.4.1)
 pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.4.1)
 pkgload             1.4.0      2024-06-28 [1] CRAN (R 4.4.1)
 plyr                1.8.9      2023-10-02 [1] CRAN (R 4.4.1)
 png               * 0.1-8      2022-11-29 [1] CRAN (R 4.4.0)
 profvis             0.3.8      2023-05-02 [1] CRAN (R 4.4.1)
 promises            1.3.0      2024-04-05 [1] CRAN (R 4.4.1)
 psych               2.4.6.26   2024-06-27 [1] CRAN (R 4.4.1)
 purrr             * 1.0.2      2023-08-10 [1] CRAN (R 4.4.1)
 qgraph            * 1.9.8      2023-11-03 [1] CRAN (R 4.4.1)
 quadprog            1.5-8      2019-11-20 [1] CRAN (R 4.4.0)
 R6                  2.5.1      2021-08-19 [1] CRAN (R 4.4.1)
 RColorBrewer      * 1.1-3      2022-04-03 [1] CRAN (R 4.4.0)
 Rcpp                1.0.12     2024-01-09 [1] CRAN (R 4.4.1)
 readr             * 2.1.5      2024-01-10 [1] CRAN (R 4.4.1)
 remotes             2.5.0      2024-03-17 [1] CRAN (R 4.4.1)
 reshape2          * 1.4.4      2020-04-09 [1] CRAN (R 4.4.1)
 rgl               * 1.3.1      2024-03-05 [1] CRAN (R 4.4.1)
 rlang               1.1.4      2024-06-04 [2] CRAN (R 4.4.1)
 rmarkdown           2.27       2024-05-17 [1] CRAN (R 4.4.1)
 rpart               4.1.23     2023-12-05 [2] CRAN (R 4.4.1)
 rphylopic         * 1.4.0      2024-04-23 [1] CRAN (R 4.4.1)
 RRPP              * 2.0.3      2024-06-22 [1] CRAN (R 4.4.1)
 rstatix             0.7.2      2023-02-01 [1] CRAN (R 4.4.1)
 rstudioapi          0.16.0     2024-03-24 [1] CRAN (R 4.4.1)
 rsvg                2.6.0      2023-10-08 [1] CRAN (R 4.4.1)
 Rvcg              * 0.23       2024-06-27 [1] CRAN (R 4.4.1)
 scales            * 1.3.0      2023-11-28 [1] CRAN (R 4.4.1)
 scatterplot3d       0.3-44     2023-05-05 [1] CRAN (R 4.4.0)
 sessioninfo         1.2.2      2021-12-06 [1] CRAN (R 4.4.1)
 shiny               1.8.1.1    2024-04-02 [1] CRAN (R 4.4.1)
 SnowballC           0.7.1      2023-04-25 [1] CRAN (R 4.4.0)
 stringi             1.8.4      2024-05-06 [2] CRAN (R 4.4.0)
 stringr           * 1.5.1      2023-11-14 [2] CRAN (R 4.4.1)
 subplex             1.9        2024-07-05 [1] CRAN (R 4.4.1)
 SURGE             * 0.1.0      2024-07-16 [1] Github (rnfelice/SURGE@225a842)
 tibble            * 3.2.1      2023-03-20 [1] CRAN (R 4.4.1)
 tidyr             * 1.3.1      2024-01-24 [1] CRAN (R 4.4.1)
 tidyselect          1.2.1      2024-03-11 [1] CRAN (R 4.4.1)
 tidytext          * 0.4.2      2024-04-10 [1] CRAN (R 4.4.1)
 tidyverse         * 2.0.0      2023-02-22 [1] CRAN (R 4.4.1)
 timechange          0.3.0      2024-01-18 [1] CRAN (R 4.4.1)
 tokenizers          0.3.0      2022-12-22 [1] CRAN (R 4.4.1)
 tzdb                0.4.0      2023-05-12 [1] CRAN (R 4.4.1)
 urlchecker          1.0.1      2021-11-30 [1] CRAN (R 4.4.1)
 usethis           * 2.2.3      2024-02-19 [1] CRAN (R 4.4.1)
 utf8                1.2.4      2023-10-22 [1] CRAN (R 4.4.1)
 vctrs               0.6.5      2023-12-01 [2] CRAN (R 4.4.1)
 withr               3.0.0      2024-01-16 [1] CRAN (R 4.4.1)
 xfun                0.45       2024-06-16 [1] CRAN (R 4.4.1)
 XML                 3.99-0.17  2024-06-25 [1] CRAN (R 4.4.1)
 xtable              1.8-4      2019-04-21 [1] CRAN (R 4.4.1)
 yulab.utils         0.1.4      2024-01-28 [1] CRAN (R 4.4.1)

 [1] C:/Users/agnL/AppData/Local/R/win-library/4.4
 [2] C:/Program Files/R/R-4.4.1/library

────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

```
