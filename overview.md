This is the accompanying interactive webserver for the Dnmt1<sup>-/-</sup>, Dnmt3a<sup>-/-</sup> and Dnmt3b<sup>-/-</sup> KO analysis of our study [Single-cell multi-omics profiling links dynamic DNA methylation to cell fate decisions during mouse early organogenesis](XXX). It covers most of Figure 1-2 as well as Supplementary Figures 1-6.


### Experimental design

We generated *Dnmt1*<sup>-/-</sup>, *Dnmt3a*<sup>-/-</sup> and *Dnmt3b*<sup>-/-</sup> embryos together with matching wildtypes from heterozygous matings. We collected embryos at E8.5, when progenitor cells for all major organs have formed and methylation mutants are not yet lethal, and performed scRNA-seq. To increase the statistical power of our analysis we combined our data set of KO embryos with a [published data set](https://www.nature.com/articles/s41586-020-2552-x) where *Dnmt1*, *Dnmt3a* and *Dnmt3b* were disrupted using CRISPR-Cas9 and also profiled using scRNA-seq at E8.5. In total, our analysis comprises 51,811 cells from 17 WT embryos, 45,579 cells from 14 *Dnmt3a*<sup>-/-</sup> embryos, 55,237 cells from 12 *Dnmt3b*<sup>-/-</sup> embryos and 25,185 cells from 15 *Dnmt1*<sup>-/-</sup> embryos. We assigned celltype labels by mapping the RNA expression profiles to a [ reference atlas that spans E6.5 to E8.5] (https://www.nature.com/articles/s41586-019-0933-9)

### Key results

- *Dnmt3a*<sup>-/-</sup> and *Dnmt3b*<sup>-/-</sup> embryos show relatively minor defects in cell type proportions. In contrast, *Dnmt1*<sup>-/-</sup> embryos show widespread defects in cell type proportions, including a relative overrepresentation of ExE ectoderm and immature embryonic cell types such as rostral neuroectoderm and caudal epiblast. We also observe a relative underrepresentation of some mature embryonic cell types, including Neural crest, NMPs, Brain, Spinal cord and Gut cells.  

- We find a small number of DE genes when comparing *Dnmt3a*<sup>-/-</sup> and *Dnmt3b*<sup>-/-</sup> to WT samples. In contrast, we find a large number of DE genes in the *Dnmt1*<sup>-/-</sup> KO across most cell types, but particularly in the Neural crest, Caudal mesoderm and Blood progenitors.  

- Posterior Hox genes (for example *Hoxc9*, *Hoxc8*, *Hoxb9* and *Hoxa9*) are downregulated in posterior cell types of *Dnmt1*<sup>-/-</sup> KO cells, including NMPs, somitic mesoderm, intermediate mesoderm and ExE mesoderm.  

- Among the genes that are upregulated in the *Dnmt1*<sup>-/-</sup> KO we observe primed pluripotency markers (*Pou5f1*, *Utf1*, *Slc7a3*, *Fgf5* and *Pim2*) and ExE marker genes (*Rhox5*, *Krt8*, *Apoe*, *Ascl2*, *Trap1a* and *Xlr3a*) across most cell types.


### Code availability

Analysis code is available at https://github.com/rargelaguet/10x_gastrulation_DNMTs.  
Code for the Shiny a is available at https://github.com/rargelaguet/shiny_dnmt.

### Data availability

Raw data is available at [GEO: GSE204908](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204908).  
Links to the parsed data objects is available in the github repository above.


### Contact

For questions on the experimental work, please reach Tim Lohoff (tlohoff431@gmail.com) or Stephen Clark (Stephen.Clark@babraham.ac.uk). For questions on the computational analysis, please reach Ricard Argelaguet (ricard.argelaguet@gmail.com).