//This workflow is used for the celloracle pipeline to build base_GRN based with the scATAC-seq data

process GRN_base_R{
    container 'kaizhang/cicero'
    
    input:
        tuple val(id), val(genome_name), path(refgenome_filefolder), path(atac_matrix), path(cell_labels), path(peaks)
    
    output:
        tuple val(id), val(genome_name), path(refgenome_filefolder), path("all_peaks.csv"), path("cicero_connections.csv")

    script:
    """
    #!/usr/bin/env Rscript

    library(cicero)
    library(monocle3)

    # Create a folder to save results
    output_folder <- "cicero_output"
    dir.create(output_folder)

    # Read in matrix data using the Matrix package
    indata <- Matrix::readMM("$atac_matrix")
    # Binarize the matrix
    indata@x[indata@x > 0] <- 1

    # Format cell info
    cellinfo <- read.table("$cell_labels")
    row.names(cellinfo) <- cellinfo\$V1
    names(cellinfo) <- "cells"

    # Format peak info
    peakinfo <- read.table("$peaks")
    names(peakinfo) <- c("chr", "bp1", "bp2")
    peakinfo\$site_name <- paste(peakinfo\$chr, peakinfo\$bp1, peakinfo\$bp2, sep="_")
    row.names(peakinfo) <- peakinfo\$site_name

    row.names(indata) <- row.names(peakinfo)
    colnames(indata) <- row.names(cellinfo)

    # Make CDS
    input_cds <-  suppressWarnings(new_cell_data_set(indata,
    cell_metadata = cellinfo,
    gene_metadata = peakinfo))

    input_cds <- monocle3::detect_genes(input_cds)

    #Ensure there are no peaks included with zero reads
    input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

    # Filter cells by peak_count
    max_count <-  15000 
    min_count <- 2000
    input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) >= min_count]
    input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) <= max_count]

    set.seed(2017)

    input_cds <- detect_genes(input_cds)
    input_cds <- estimate_size_factors(input_cds)
    input_cds <- preprocess_cds(input_cds, method = "LSI")

    # Dimensional reduction with umap
    input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")    
    umap_coords <- reducedDims(input_cds)\$UMAP


    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
    chromosome_length <- read.table(paste0("$refgenome_filefolder", "/", "$genome_name", ".fa.sizes"), header = F)
    #length <- list.files("$refgenome_filefolder", pattern = "\\\\.sizes\$", full.names = TRUE)
    #chromosome_length <- read.table(length[1], header = FALSE)

    conns <- run_cicero(cicero_cds, chromosome_length)
    all_peaks <- row.names(exprs(input_cds))
    write.csv(x = all_peaks, file = "all_peaks.csv")
    write.csv(x = conns, file = "cicero_connections.csv")
    """
}

process GRN_base_py{
    label 'py'

    input:
        tuple val(id), val(genome_name), path(refgenome_filefolder), path(all_peaks), path(cicero_connections)

    output:
        tuple val(id), path("base_GRN_dataframe.parquet")
    
    script:
    """
    #!/usr/bin/env python3
    import os
    os.environ['HOME'] = os.getcwd()
    
    import pandas as pd
    from celloracle import motif_analysis as ma
    import os

    ref_genome = "$genome_name"

    # Load scATAC-seq peak list.
    peaks = pd.read_csv("$all_peaks", index_col=0)
    peaks = peaks.x.values
    cicero_connections = pd.read_csv("$cicero_connections", index_col=0)

    ##!! Please make sure to specify the correct reference genome here
    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=ref_genome)

    integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, 
                                                cicero_connections=cicero_connections)
    peak = integrated[integrated.coaccess >= 0.8]
    peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)

    peak = pd.DataFrame(peak, columns=["peak_id", "gene_short_name"])

    #regdir = os.path.dirname($refgenome_filefolder)

    # Instantiate TFinfo object
    tfi = ma.TFinfo(peak_data_frame=peak, ref_genome=ref_genome,genomes_dir="$refgenome_filefolder") 
    # Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
    tfi.scan(fpr=0.02, motifs=None, verbose=True)
    # If you enter None, default motifs will be loaded.

    # Save tfinfo object
    tfi.to_hdf5(file_path="tfi.celloracle.tfinfo")

    # Reset filtering 
    tfi.reset_filtering()

    # Do filtering
    tfi.filter_motifs_by_score(threshold=10)

    # Format post-filtering results.
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
    df = tfi.to_dataframe()
    df.to_parquet("base_GRN_dataframe.parquet")
    """

}

workflow base_GRN{
    take:
        ATAC_dataset
    
    main:
        base_GRN=
            ATAC_dataset 
            | map { [it.id, it.genome_name, it.genome_filefolder, it.atac_matrix, it.cell_labels, it.peaks] }
            | GRN_base_R 
            | GRN_base_py

    emit:
        base_GRN
    
}