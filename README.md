# Brain_WGCNA
## Scripts for the WGCNA paper: 
## *Schizophrenia risk genes converge into shifting co-expression networks across brain development, ageing and brain regions*

### For the "Regional-coexpression" and "Age-period" studies (NC samples only): 
* [Code1(NC).R](Code1(NC).R): Preprocess tissue rse data
* [Code2(NC).R](Code2(NC).R): Remove confounders effect
* [Code3(NC).R](Code3(NC).R): Estimate beta for WGCNA via connectivity match
* [Code4(NC).R](Code4(NC).R): Perform WGCNA to obtain age-parsed and non-parsed networks 
* [prepare_WGCNAnetwork_script.R](combine/prepare_WGCNAnetwork_script.R): Combine WGCNA network output (module assignments, connectivity etc) in a single file
* [prepare_OTHERnetwork_script.R](combine/prepare_OTHERnetwork_script.R): Preprocess and collate network data from various published networks
* [prepare_ALL_WGCNA_network_script.R](combine/prepare_ALL_WGCNA_network_script.R): Combine WGCNA output from our networks and other previous published networks
* [new_wide_form_data_script.R](combine/new_wide_form_data_script.R): Create a combined wide_data_file with WGCNA output from our and other networks, plus genestats, MAGMA stats etc
* [get_gene_module_list script.R](combine/get_gene_module_list%20script.R): Get a gene_module_list for all networks from the wide_data_file to be used for enrichment analysis
* [Enrichments.R](enrich_and_plot/Enrichments.R): Script to run various enrichments on gene_module_list: SCZ, DEGs, DMGs, TWAS, GO, MAGMA, Cell Specificity etc. Script to make sankey plots for DLPFC, HP and CN SCZ risk genes. Script to make binwise PGC3 enrichment scatterplots for networks. Prepare and format data for other visualisations
* [Enrichment Visualisation.R](enrich_and_plot/Enrichment%20Visualisation.R): Script to plot main enrichment plots for SCZ risk modules for all networks. Script to plot GO enrichment plots for SCZ risk modules for all networks. Script to plot TF enrichment plots for SCZ risk modules for all networks. Script to plot other Pathology enrichment plots for SCZ risk modules for all networks

### For the "Cell population enrichment" study:
* [matched-sample-HP-DG(Code1).R](cell_population_enrichment/matched-sample-HP-DG(Code1).R): Preprocess sample match rse data for HP-DG
* [matched-sample-HP-DG(Code2,noqsva).R](cell_population_enrichment/matched-sample-HP-DG(Code2%2Cnoqsva).R): Remove confounders effect
* [matched-sample-HP-DG(Code2,qsva).R](cell_population_enrichment/matched-sample-HP-DG(Code2%2Cqsva).R):  Remove confounders effect + QSV
* [matched-sample-HP-DG(Code3,noqsva).R](cell_population_enrichment/matched-sample-HP-DG(Code3%2Cnoqsva).R): Estimate beta for WGCNA via connectivity match 
* [matched-sample-HP-DG(Code3,qsva).R](cell_population_enrichment/matched-sample-HP-DG(Code3%2Cqsva).R): Estimate beta for WGCNA via connectivity match [QSVA pipeline]
* [matched-sample-HP-DG(Code4,noqsva).R](cell_population_enrichment/matched-sample-HP-DG(Code4%2Cnoqsva).R): Perform WGCNA to obtain age-parsed and non-parsed networks
* [matched-sample-HP-DG(Code4,qsva).R](cell_population_enrichment/matched-sample-HP-DG(Code4%2Cqsva).R): Perform WGCNA to obtain age-parsed and non-parsed networks [QSVA pipeline]
