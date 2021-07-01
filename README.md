# FlowSOM protocol
R code to demonstrate the FlowSOM analysis pipeline.   
The protocol, including installing the necessary packages and downloading the 
used dataset, can be found in [R/FlowSOM_protocol.R](R/FlowSOM_protocol.R). 
Typically, the installation of the packages takes less than ten minutes. An 
average FlowSOM analysis takes one to three hours to complete, quality issues 
can increase the time considerably.  
This automated analysis pipeline will result in a FlowSOM object 
and multiple plots that will help in understanding, visualizing and analyzing 
high-dimensional data.

### Key steps in the protocol
* Prepare the data, including quality control (steps 1-16)
* Train FlowSOM model (steps 17-19)
* Check the quality of the FlowSOM model (step 20)
* Use the FlowSOM model for further discovery (steps 21-26)
![FlowSOM analysis pipeline](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41596-021-00550-0/MediaObjects/41596_2021_550_Fig1_HTML.png?as=webp)  

### References
Quintelier, K. et al. Analyzing high-dimensional cytometry data using FlowSOM. 
Nat Protoc 1–27 (2021) doi:10.1038/s41596-021-00550-0. (https://rdcu.be/cndgZ)  
(Additional code to reproduce the figures and timing is listed in 
[R/Figures_and_timing.R](R/Figures_and_timing.R).)  
  
Van Gassen, S. et al. FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data. Cytometry A 87, 636–645 (2015) doi:10.1002/cyto.a.22625.



### Session information
```{r}
sessionInfo()
```

```
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.3.3  FlowSOM_2.1.10 igraph_1.2.6   flowCore_2.2.0

loaded via a namespace (and not attached):
  [1] Rtsne_0.15                  ggnewscale_0.4.5            colorspace_2.0-0            ggsignif_0.6.1              rjson_0.2.20               
  [6] ellipsis_0.3.1              rio_0.5.16                  colorRamps_2.3              rsconnect_0.8.16            rprojroot_2.0.2            
 [11] circlize_0.4.12             cytolib_2.2.1               GlobalOptions_0.1.2         base64enc_0.1-3             clue_0.3-58                
 [16] rstudioapi_0.13             ggpubr_0.4.0                farver_2.0.3                hexbin_1.28.2               CytoML_2.2.1               
 [21] ggrepel_0.9.1               fansi_0.4.2                 xml2_1.3.2                  knitr_1.31                  polyclip_1.10-0            
 [26] jsonlite_1.7.2              Cairo_1.5-12.2              broom_0.7.5                 cluster_2.1.0               png_0.1-7                  
 [31] pheatmap_1.0.12             graph_1.68.0                ggforce_0.3.2               compiler_4.0.3              httr_1.4.2                 
 [36] backports_1.2.1             assertthat_0.2.1            tweenr_1.0.1                htmltools_0.5.1.1           tools_4.0.3                
 [41] ncdfFlow_2.36.0             gtable_0.3.0                glue_1.4.2                  flowWorkspace_4.2.0         dplyr_1.0.4                
 [46] ggcyto_1.18.0               Rcpp_1.0.6                  scattermore_0.7             carData_3.0-4               Biobase_2.50.0             
 [51] cellranger_1.1.0            vctrs_0.3.6                 xfun_0.21                   openxlsx_4.2.3              lifecycle_1.0.0            
 [56] rstatix_0.7.0               XML_3.99-0.5                zlibbioc_1.36.0             MASS_7.3-53                 scales_1.1.1               
 [61] RProtoBufLib_2.2.0          hms_1.0.0                   parallel_4.0.3              RBGL_1.66.0                 RColorBrewer_1.1-2         
 [66] ComplexHeatmap_2.6.2        yaml_2.2.1                  curl_4.3                    aws.signature_0.6.0         gridExtra_2.3              
 [71] latticeExtra_0.6-29         stringi_1.5.3               S4Vectors_0.28.1            PeacoQC_1.4.0               BiocGenerics_0.36.0        
 [76] zip_2.1.1                   shape_1.4.5                 rlang_0.4.10                pkgconfig_2.0.3             matrixStats_0.58.0         
 [81] bitops_1.0-6                evaluate_0.14               lattice_0.20-41             purrr_0.3.4                 labeling_0.4.2             
 [86] cowplot_1.1.1               tidyselect_1.1.0            plyr_1.8.6                  magrittr_2.0.1              FlowRepositoryR_1.20.0     
 [91] R6_2.5.0                    IRanges_2.24.1              generics_0.1.0              DBI_1.1.0                   pillar_1.5.0               
 [96] haven_2.3.1                 foreign_0.8-80              withr_2.4.1                 abind_1.4-5                 RCurl_1.98-1.2             
[101] tibble_3.0.6                crayon_1.4.1                car_3.0-10                  utf8_1.1.4                  ggpointdensity_0.1.0       
[106] rmarkdown_2.6               aws.s3_0.3.21               jpeg_0.1-8.1                GetoptLong_1.0.5            grid_4.0.3                 
[111] readxl_1.3.1                data.table_1.13.6           Rgraphviz_2.34.0            ConsensusClusterPlus_1.54.0 forcats_0.5.1              
[116] digest_0.6.27               tidyr_1.1.2                 RcppParallel_5.0.2          stats4_4.0.3                munsell_0.5.0    
```
