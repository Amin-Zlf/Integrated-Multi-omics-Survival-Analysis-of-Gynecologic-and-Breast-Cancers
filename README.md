# Integrated Multi-omics Survival Analysis of Gynecologic and Breast Cancers
This study aimed to conduct a multi-omics survival analysis, incorporating gene expression, methylation, copy number variation (CNV), and mutation data, on Breast Carcinoma (BRCA) and Gynecologic Cancers (Ovarian Serous Cystadenocarcinoma (OV), Uterine Corpus Endometrial Carcinoma (UCEC), Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESC), and Uterine Carcinosarcoma (UCS)). The initial objective was to identify pathways significantly associated with patients’ survival. Subsequent analyses sought to determine which modules within each significant pathway and which genes within each significant module were significantly associated with survival. The MOSClip R package was used to perform these analyses.
MOSClip is an R package that facilitates integrated multi-omic survival analysis through the use of multivariate models and dimensionality reduction of multi-omics data (https://doi.org/10.1093/nar/gkz324). This topological pathway analysis tool can identify significant pathways, modules, and genes in survival analysis. MOSClip was chosen for this study because it is currently the only available tool capable of performing survival analysis using multi-omics data while accounting for interactions among genes. To be more specific, MOSClip is a versatile approach that employs omic-specific dimensionality reduction techniques to conduct multi-omic survival tests on individual pathways or modules derived from a pathway’s graph structure. These tests can be performed independently at the pathway or module level using a multivariate survival model to identify associations with patient survival. 
Cytoscape was then used to visualize the topology of significant genes within each module. Through this analysis, 33 genes were identified as being common among different types of cancers. Afterward, heatmaps were created for each cancer type to illustrate the effect of significant genes on patient survival. In these heatmaps, the genes are carefully ordered based on their relative importance in determining patient survival. The score values for each cancer type are then meticulously divided into three distinct ranges, representing the low, medium, or high risk of death for each individual patient. Then, Kaplan-Meier plots were compared among different types of cancers to provide valuable insights into the survival rates and differences among cancer types. The underlying rationale for this approach was that if the previous tests and analyses had been performed correctly, the resulting Kaplan-Meier plot for each cancer type should accurately represent the relative survival rates of the different risk groups. The Kaplan-Meier plots for all five cancer types demonstrated that the separation of patients based on their death risk was performed correctly. An additional test was performed to assess the accuracy of survival prediction, with rates of 69.56522% for BRCA, 94.73684% for CESC, 94.44444% for OV, 81.81818% for UCEC, and 66.66667% for UCS.
In summary, this research provides valuable insights that may inform future research and treatment strategies. We hope that these findings will serve as a foundation for further investigation into the underlying mechanisms and potential therapeutic targets for these cancers. 

