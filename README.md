# Brain_WGCNA
## Scripts for the manuscript: 
## *Consensus molecular environment of schizophrenia risk genes in co-expression networks shifting across age and brain regions*
"Giulio Pergola et al. ,Consensus molecular environment of schizophrenia risk genes in coexpression networks shifting across age and brain regions.Sci. Adv.9, eade2812(2023).DOI:10.1126/sciadv.ade2812"

### Scripts for the "Fixed-windows" and "Sliding-windows" studies (NC samples only):
* [Code1(NC).R](revision/preprocess/revision_script1.R): Preprocess tissue rse data for NC samples
* [Code2(NC).R](revision/preprocess/revision_script2.R): Remove confounders effect
* [Code3(NC).R](revision/preprocess/revision_script3.R): Estimate beta for WGCNA via connectivity match
* [Code4(NC).R](revision/preprocess/revision_script4.R): Perform WGCNA to obtain age-parsed and non-parsed networks 
* [prepare_WGCNAnetwork_script.R](combine/prepare_WGCNAnetwork_script.R): Combine WGCNA network output (module assignments, connectivity etc) in a single file
* [prepare_OTHERnetwork_script.R](combine/prepare_OTHERnetwork_script.R): Preprocess and collate network data from various published networks
* [prepare_ALL_WGCNA_network_script.R](combine/prepare_ALL_WGCNA_network_script.R): Combine WGCNA output from our networks and other previous published networks
* [new_wide_form_data_script.R](combine/new_wide_form_data_script.R): Create a combined wide_data_file with WGCNA output from our and other networks, plus genestats, MAGMA stats etc
* [get_gene_module_list script.R](combine/get_gene_module_list%20script.R): Get a gene_module_list for all networks from the wide_data_file to be used for enrichment analysis
* [Enrichments.online.script.R](revision/enrich_and_plot/Enrichments.online.script.R):
  - Script to run various enrichments on gene_module_list: SCZ, DEGs, DMGs, TWAS, GO, MAGMA, Cell Specificity etc
  - Script to make sankey plots for DLPFC, HP and CN SCZ risk genes
  - Prepare and format data for other visualisations
* [Enrichment Visualisation.online.script.R](revision/enrich_and_plot/Enrichment%20Visualisation.online.script.R): 
  - Script to plot main enrichment plots for SCZ risk modules for all networks
  - Script to plot GO enrichment plots for SCZ risk modules for all networks
  - Script to plot TF enrichment plots for SCZ risk modules for all networks
  - Script to plot other Pathology enrichment plots for SCZ risk modules for all networks

### Scripts for the "MAGMA prediction model":
* [MAGMA_prediction.py](MLcode/MAGMA_prediction.py): Code to perfrom ML MAGMA prediction for Age-parsed and non-parsed networks

### Scripts for the "Cell population enrichment" study:
* [matched-sample-HP-DG(Code1).R](cell_population_enrichment/matched-sample-HP-DG(Code1).R): Preprocess samples match rse data for HP-DG
* [matched-sample-HP-DG(Code2,noqsva).R](cell_population_enrichment/matched-sample-HP-DG(Code2%2Cnoqsva).R): Remove confounders effect
* [matched-sample-HP-DG(Code2,qsva).R](cell_population_enrichment/matched-sample-HP-DG(Code2%2Cqsva).R):  Remove confounders effect + QSV
* [matched-sample-HP-DG(Code3,noqsva).R](cell_population_enrichment/matched-sample-HP-DG(Code3%2Cnoqsva).R): Estimate beta for WGCNA via connectivity match 
* [matched-sample-HP-DG(Code3,qsva).R](cell_population_enrichment/matched-sample-HP-DG(Code3%2Cqsva).R): Estimate beta for WGCNA via connectivity match [QSVA pipeline]
* [matched-sample-HP-DG(Code4,noqsva).R](cell_population_enrichment/matched-sample-HP-DG(Code4%2Cnoqsva).R): Perform WGCNA to obtain samples match HP-DG networks
* [matched-sample-HP-DG(Code4,qsva).R](cell_population_enrichment/matched-sample-HP-DG(Code4%2Cqsva).R): Perform WGCNA to obtain samples match HP-DG networks [QSVA pipeline]


### Scripts for the "Sliding-windows" study:
* [sliding_window(Code1,SCZ).R](sliding_windows/sliding_window(Code1%2CSCZ).R): Preprocess tissue rse data (SCZ samples only)
* [sliding_window(Code2,SCZ).R](sliding_windows/sliding_window(Code2%2CSCZ).R): Remove confounders effect
* [sliding_window(Code3,NC).R](sliding_windows/sliding_window(Code3%2CNC).R): Arrange samples in sliding age windows and estimate beta for WGCNA via connectivity match (NC samples only)
* [sliding_window(Code3,SCZ).R](sliding_windows/sliding_window(Code3%2CSCZ).R): Arrange samples in sliding age windows and estimate beta for WGCNA via connectivity match (SCZ samples only)
* [sliding_window(Code4,NC).R](sliding_windows/sliding_window(Code4%2CNC).R): Perform WGCNA to obtain sliding window networks on the NC samples
* [sliding_window(Code4,SCZ).R](sliding_windows/sliding_window(Code4%2CSCZ).R): Perform WGCNA to obtain sliding window networks on the SCZ samples
* [Enrichments and Visualisation(Sliding Window).R](sliding_windows/enrich_and_plot/Enrichments%20and%20Visualisation(Sliding%20Window).R): Script to plot SCZ enrichment measures for sliding window networks (NC + SCZ)
* [Visualisation(Sliding Window).R](sliding_windows/enrich_and_plot/Visualisation(Sliding%20Window).R): Script for additional sliding window plots (max.Fold change vs median age with rugs)


### Scripts for the "Sliding-windows enrichment":
* [module_enrichments_MAGMA.R](sliding_windows_enrichment/module_enrichments_MAGMA.R): Calculate enrichment of MAGMA per module using genesettest 
* [module_enrichments_DRONC.R](sliding_windows_enrichment/module_enrichments_DRONC.R): Calculate enrichment of DRONC celltypes per module using genesettest 
* [module_enrichments_GO.R](sliding_windows_enrichment/module_enrichments_GO.R): Calculate GO term ratio per module
* [combine_MAGMA.R](sliding_windows_enrichment/combine_MAGMA.R): Add MAGMA module enrichments to GO term ratio and DRONC celltype enrichment dataframes 
* [module_enrichments_p_values_DRONCtoMAGMA.R](sliding_windows_enrichment/module_enrichments_p_values_DRONCtoMAGMA.R): Calculate across module association of MAGMA enrichment to cell type enrichment
* [module_enrichments_p_values_GOtoMAGMA.R](sliding_windows_enrichment/module_enrichments_p_values_GOtoMAGMA.R): Calculate across module association of MAGMA enrichment to GO term ratio


### Scripts for the "replication in human iPSC data":
* [stemcell(Code1).R](stemcells/stemcell(Code1).R): Prepare and preprocess iPSC data
* [stemcell(Code2).R](stemcells/stemcell(Code2).R): Remove confounders effect
* [stemcell(Code3)[noKO].R](stemcells/stemcell(Code3)%5BnoKO%5D.R): Estimate beta for WGCNA via connectivity match
* [stemcell(Code4)[noKO].R](stemcells/stemcell(Code4)%5BnoKO%5D.R): Perform WGCNA to obtain iPSC networks [noKO]
* [Enrichment and Visualisation (stemcells).R](stemcells/enrich_and_plot/Enrichment%20and%20Visualisation%20(stemcells).R): Script to plot enrichment in PGC lists vs enrichment in consensus list for the iPSC networks [noKO and shuffled/KO]

### Scripts for the "Consensus environment":
* [Consensus environment enrichment for GO and SCZ.R](Consensus%20environment/Consensus%20environment%20enrichment%20for%20GO%20and%20SCZ.R): Run GO and SCZ enrichments on consensus gene environment
* [Consensus environment enrichment for KEGG.R](Consensus%20environment/Consensus%20environment%20enrichment%20for%20KEGG.R): Run KEGG enrichment on consensus gene environment

### Uncategorized scripts:
* [Miscellenous scripts](misc/scripts/): Consensus genes computation, Jaccard index calculations etc

### Additional links:
* [Zenodo](https://doi.org/10.5281/zenodo.5676480): Interactive Sankey files, wide_form data files, modulewise SCZ enrichment results, preprocessed data files etc
* [NETS@LIBD](https://nets.libd.org/age_wgcna/): Future project updates and general updates of our research group at Lieber Institute for Brain Development will be available here.

For any data or code inquiries please contact Giulio Pergola: [Giulio.Pergola@libd.org]
