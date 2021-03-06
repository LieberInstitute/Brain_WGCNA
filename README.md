# Brain_WGCNA
## Scripts for the WGCNA paper: 
## *A consensus molecular environment of schizophrenia risk genes in shifting co-expression networks across brain development, ageing and brain regions*

### For the "Regional-coexpression" and "Age-period" studies (NC samples only): 
* [Code1(NC).R](preprocess/Code1(NC).R): Preprocess tissue rse data
* [Code2(NC).R](preprocess/Code2(NC).R): Remove confounders effect
* [Code3(NC).R](preprocess/Code3(NC).R): Estimate beta for WGCNA via connectivity match
* [Code4(NC).R](preprocess/Code4(NC).R): Perform WGCNA to obtain age-parsed and non-parsed networks 
* [prepare_WGCNAnetwork_script.R](combine/prepare_WGCNAnetwork_script.R): Combine WGCNA network output (module assignments, connectivity etc) in a single file
* [prepare_OTHERnetwork_script.R](combine/prepare_OTHERnetwork_script.R): Preprocess and collate network data from various published networks
* [prepare_ALL_WGCNA_network_script.R](combine/prepare_ALL_WGCNA_network_script.R): Combine WGCNA output from our networks and other previous published networks
* [new_wide_form_data_script.R](combine/new_wide_form_data_script.R): Create a combined wide_data_file with WGCNA output from our and other networks, plus genestats, MAGMA stats etc
* [get_gene_module_list script.R](combine/get_gene_module_list%20script.R): Get a gene_module_list for all networks from the wide_data_file to be used for enrichment analysis
* [Enrichments.R](enrich_and_plot/Enrichments.R):
  - Script to run various enrichments on gene_module_list: SCZ, DEGs, DMGs, TWAS, GO, MAGMA, Cell Specificity etc
  - Script to make sankey plots for DLPFC, HP and CN SCZ risk genes
  - Script to make binwise PGC3 enrichment scatterplots for networks
  - Prepare and format data for other visualisations
* [Enrichment Visualisation.R](enrich_and_plot/Enrichment%20Visualisation.R): 
  - Script to plot main enrichment plots for SCZ risk modules for all networks
  - Script to plot GO enrichment plots for SCZ risk modules for all networks
  - Script to plot TF enrichment plots for SCZ risk modules for all networks
  - Script to plot other Pathology enrichment plots for SCZ risk modules for all networks


### For the "Cell population enrichment" study:
* [matched-sample-HP-DG(Code1).R](cell_population_enrichment/matched-sample-HP-DG(Code1).R): Preprocess samples match rse data for HP-DG
* [matched-sample-HP-DG(Code2,noqsva).R](cell_population_enrichment/matched-sample-HP-DG(Code2%2Cnoqsva).R): Remove confounders effect
* [matched-sample-HP-DG(Code2,qsva).R](cell_population_enrichment/matched-sample-HP-DG(Code2%2Cqsva).R):  Remove confounders effect + QSV
* [matched-sample-HP-DG(Code3,noqsva).R](cell_population_enrichment/matched-sample-HP-DG(Code3%2Cnoqsva).R): Estimate beta for WGCNA via connectivity match 
* [matched-sample-HP-DG(Code3,qsva).R](cell_population_enrichment/matched-sample-HP-DG(Code3%2Cqsva).R): Estimate beta for WGCNA via connectivity match [QSVA pipeline]
* [matched-sample-HP-DG(Code4,noqsva).R](cell_population_enrichment/matched-sample-HP-DG(Code4%2Cnoqsva).R): Perform WGCNA to obtain samples match HP-DG networks
* [matched-sample-HP-DG(Code4,qsva).R](cell_population_enrichment/matched-sample-HP-DG(Code4%2Cqsva).R): Perform WGCNA to obtain samples match HP-DG networks [QSVA pipeline]


### For the "Sliding age windows" study:
* [sliding_window(Code1,SCZ).R](sliding_windows/sliding_window(Code1%2CSCZ).R): Preprocess tissue rse data (SCZ samples only)
* [sliding_window(Code2,SCZ).R](sliding_windows/sliding_window(Code2%2CSCZ).R): Remove confounders effect
* [sliding_window(Code3,NC).R](sliding_windows/sliding_window(Code3%2CNC).R): Arrange samples in sliding age windows and estimate beta for WGCNA via connectivity match (NC samples only)
* [sliding_window(Code3,SCZ).R](sliding_windows/sliding_window(Code3%2CSCZ).R): Arrange samples in sliding age windows and estimate beta for WGCNA via connectivity match (SCZ samples only)
* [sliding_window(Code4,NC).R](sliding_windows/sliding_window(Code4%2CNC).R): Perform WGCNA to obtain sliding window networks on the NC samples
* [sliding_window(Code4,SCZ).R](sliding_windows/sliding_window(Code4%2CSCZ).R): Perform WGCNA to obtain sliding window networks on the SCZ samples
* [Enrichments and Visualisation(Sliding Window).R](sliding_windows/enrich_and_plot/Enrichments%20and%20Visualisation(Sliding%20Window).R): Script to plot SCZ enrichment measures for sliding window networks (NC + SCZ)


### For the iPSC data analysis:
* [stemcell(Code1).R](stemcells/stemcell(Code1).R): Prepare and preprocess iPSC data
* [stemcell(Code2).R](stemcells/stemcell(Code2).R): Remove confounders effect
* [stemcell(Code3)[noKO].R](stemcells/stemcell(Code3)%5BnoKO%5D.R): Estimate beta for WGCNA via connectivity match
* [stemcell(Code3)[shuffleKO].R](stemcells/stemcell(Code3)%5BshuffleKO%5D.R): Make shuffled/KO expression and estimate beta for WGCNA via connectivity match
* [stemcell(Code4)[noKO].R](stemcells/stemcell(Code4)%5BnoKO%5D.R): Perform WGCNA to obtain iPSC networks [noKO]
* [stemcell(Code4)[shuffleKO].R](stemcells/stemcell(Code4)%5BshuffleKO%5D.R): Perform WGCNA to obtain iPSC networks [shuffled/KO]
* [Enrichment and Visualisation (stemcells).R](stemcells/enrich_and_plot/Enrichment%20and%20Visualisation%20(stemcells).R): Script to plot enrichment in PGC lists vs enrichment in consensus list for the iPSC networks [noKO and shuffled/KO]



### Additional data files:
* [External (Zenodo) link](https://zenodo.org/record/5706100): Interactive Sankey files, wide_form data files etc
