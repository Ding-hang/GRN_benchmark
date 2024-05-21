include { base_GRN } from './base_GRN/base_GRN.nf'


process genome_download {

    input:
        val genome_name
    
    output:
        tuple val(genome_name), path(genome_name)

    script:
        """
        #!/usr/bin/env python3

        import genomepy

        genomepy.install_genome(name="$genome_name", provider="UCSC",genomes_dir="./genome/")

        """ 

}

process RNA_pre {
    label 'py'
    
    input:
        tuple val(id), val(genome_name), path(refgenome_filefolder),  path(RNA_expression)
        
    output:
        tuple val(id), path("RNA_etal_15.h5ad")

    script:
        """
        #!/usr/bin/env python3

        import scanpy as sc

        # loading dataset.
        adata = sc.read_h5ad("$RNA_expression")

        # Only consider genes with more than 1 count
        sc.pp.filter_genes(adata, min_counts=1)

        # Normalize gene expression matrix with total UMI count per cell
        sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

        # Select top 2000 highly-variable genes
        filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                                    flavor='cell_ranger',
                                                    n_top_genes=2000,
                                                    log=False)

        # Subset the genes
        adata = adata[:, filter_result.gene_subset].copy()

        # Renormalize after filtering
        sc.pp.normalize_per_cell(adata)

        # keep raw cont data before log transformation
        adata.raw = adata
        adata.layers["raw_count"] = adata.raw.X.copy()


        # Log transformation and scaling
        sc.pp.log1p(adata)
        sc.pp.scale(adata)

        # PCA
        sc.tl.pca(adata, svd_solver='arpack')

        # Diffusion map
        sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)

        sc.tl.diffmap(adata)
        # Calculate neihbors again based on diffusionmap 
        sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

        sc.tl.louvain(adata, resolution=0.8)

        # PAGA graph construction
        sc.tl.paga(adata, groups='louvain')
        sc.pl.paga(adata)
        sc.tl.draw_graph(adata, init_pos='paga', random_state=123)
        sc.pl.draw_graph(adata, color='louvain', legend_loc='on data')
        adata.write_h5ad("RNA_etal_15.h5ad")
        
        """

}

process GRN_pruing_weight{
    label 'py'

    input:
        tuple val(id), path (baseNet), path (RNA_adata)
    
    output:
        tuple val(id), path("GRN.celloracle.oracle")

    script:
        """
        #!/usr/bin/env python3
        
        import scanpy as sc
        import pandas as pd
        import numpy as np

        import celloracle as co

        adata = sc.read_h5ad("$RNA_adata")

        '''
        n_cells_downsample = 30000 
        if adata.shape[0] > n_cells_downsample:
            # Let's dowmsample into 30K cells
            sc.pp.subsample(adata, n_obs=n_cells_downsample, random_state=123)
        '''
        base_GRN = pd.read_parquet("$baseNet")
        #base_GRN = pd.read_parquet("/storage/zhangkaiLab/dingyihang/benchmark/cellOracle/data_folder/mm9_mouse_atac_atlas_data_TSS_and_cicero_0.9_accum_threshold_10.5_DF_peaks_by_TFs_v202204.parquet")
        #base_GRN = tfi.to_dataframe()

        #Object
        oracle = co.Oracle()
        adata.X = adata.layers["raw_count"].copy()

        # Instantiate Oracle object.
        oracle.import_anndata_as_raw_count(adata=adata,
                                        cluster_column_name="louvain",
                                        embedding_name="X_draw_graph_fr")
        oracle.import_TF_data(TF_info_matrix=base_GRN)

        '''
        Add TF-target gene pair manually
        '''

        #KNN imputation
        oracle.perform_PCA()
        n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
        n_comps = min(n_comps, 50)
        n_cell = oracle.adata.shape[0]
        k = int(0.025*n_cell)
        oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                            b_maxl=k*4, n_jobs=4)

        oracle.to_hdf5("GRN.celloracle.oracle")

        #oracle = co.load_hdf5("Paul_15_data.celloracle.oracle")

        #get GRN
        #links = oracle.get_links(cluster_name_for_GRN_unit="louvain_annot", alpha=10,
        #                         verbose_level=10)
        #links.links_dict.keys()
        #links.links_dict["Ery_0"]
        # Set cluster name
        #cluster = "Ery_0"
        # Save as csv
        #links.links_dict[cluster].to_csv(f"raw_GRN_for_{cluster}.csv")
        
        """

}

workflow {

    download = false
    if (download) {
        genomes = Channel.fromList(["hg19"]) | genome_download
    }
    else {
        genomes = Channel.of(["hg19", file("./genome/hg19")])
    }

    rawdata = Channel.of(
        [
            "id": "Buenrostro",
            "genome_name": "hg19",
            "atac_matrix": file("./data/ATAC_GSE96769_Buenrostro/ATAC_spares_matrix.mtx"), 
            "cell_labels": file("./data/ATAC_GSE96769_Buenrostro/ATAC_Cell_labels.txt"), 
            "peaks": file("./data/ATAC_GSE96769_Buenrostro/peaks.bed"), 
            "RNA_expression":file("./data/Paul_raw_15.h5ad")
        ]
    )   
        | map { [it.genome_name,it] }
    
    dataset = genomes.join(rawdata) | map { genome_name, genome_dir, rawdata -> rawdata['genome_filefolder'] = genome_dir; rawdata }

    baseNet = dataset | base_GRN

    RNA_adata = dataset 
        | map { [it.id, it.genome_name, it.genome_filefolder, it.RNA_expression] } 
        | RNA_pre
    
    baseNet.join(RNA_adata) | GRN_pruing_weight

}
