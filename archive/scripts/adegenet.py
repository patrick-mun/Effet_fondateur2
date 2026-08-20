# adegenet.py
"""
Module adegenet : génération et exécution d'un script R pour analyse DAPC, structure, CA et arbre phylogénétique via adegenet
"""
import os
import subprocess
import logging

def generate_adegenet_script(raw_path: str, group_file: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    script_path = os.path.join(output_dir, "adegenet_dapc.R")
    pdf_output = os.path.join(output_dir, "dapc_plot.pdf")
    tree_output = os.path.join(output_dir, "dapc_tree.pdf")
    hexp_output = os.path.join(output_dir, "dapc_hexp.pdf")
    ca_output = os.path.join(output_dir, "genpop_CA_plot.pdf")
    tree_group_output = os.path.join(output_dir, "genpop_group_trees.pdf")
    genin_path = os.path.join(output_dir, "genin.csv")
    genpop_path = os.path.join(output_dir, "genpop.txt")
    debug_path = os.path.join(output_dir, "debug_adegenet.log")

    script_content = f"""
    library(adegenet)
    library(ggplot2)
    library(ape)
    library(ade4)

    log <- file("{debug_path}", open="wt")
    writeLines("[DEBUG] Début de l'analyse DAPC", log)

    geno <- read.table('{raw_path}', header=TRUE)
    groups <- read.table('{group_file}', header=FALSE)
    data <- geno[, 7:ncol(geno)]
    rownames(data) <- geno$IID
    writeLines(paste("Colonnes initiales:", ncol(data)), log)

    colnames(data) <- sub("_.*", "", colnames(data))

    grps <- groups$V2[match(geno$IID, groups$V1)]
    valid_rows <- which(!is.na(grps))
    data <- data[valid_rows, ]
    grps <- as.factor(grps[valid_rows])
    rownames(data) <- geno$IID[valid_rows]
    writeLines(paste("Lignes valides (groupes):", nrow(data)), log)

    complete_idx <- complete.cases(data)
    data <- data[complete_idx, ]
    grps <- grps[complete_idx]
    writeLines(paste("Lignes complètes:", nrow(data)), log)

    if (ncol(data) > 0) {{
        data <- data[, colSums(is.na(data)) < nrow(data), drop=FALSE]
    }}
    writeLines(paste("Colonnes après filtre NA:", ncol(data)), log)

    if (ncol(data) > 0) {{
        colnames(data) <- paste0("locus", seq_len(ncol(data)))
    }} else {{
        writeLines("[ERREUR] Aucune colonne de génotypes valide disponible après filtrage.", log)
        close(log)
        stop("Aucune colonne de génotypes valide disponible après filtrage.")
    }}

    write.table(data, file="{genin_path}", sep=",", quote=FALSE, col.names=NA)
    write.table(grps, file="{genpop_path}", row.names=FALSE, col.names=FALSE, quote=FALSE)

    genind_obj <- df2genind(data, ploidy=2, sep="/", ind.names=rownames(data))
    pop(genind_obj) <- grps

    safe_pca <- min(50, floor(nrow(data) / 3))
    dapc_res <- dapc(genind_obj, n.pca=safe_pca, n.da=2)
    pdf('{pdf_output}')
    scatter(dapc_res)
    dev.off()

    pdf('{tree_output}')
    myTree <- nj(dist(genind_obj))
    plot(myTree, main="Neighbour-Joining Tree")
    dev.off()

    pdf('{hexp_output}')
    basic_stats <- summary(genind_obj)
    barplot(basic_stats$Hobs, beside=TRUE, col="skyblue", main="Hobs par population")
    barplot(basic_stats$Hexp, beside=TRUE, col="orange", main="Hexp par population")
    dev.off()

    genpop_obj <- genind2genpop(genind_obj)
    pdf('{ca_output}')
    ca_res <- dudi.coa(tab(genpop_obj), scannf=FALSE, nf=3)
    nf <- ncol(ca_res$li)
    writeLines(paste("Composantes disponibles pour CA:", nf), log)
    if (nf >= 2) {{
        fac_vector <- factor(rownames(ca_res$li))
        s.class(dfxy = ca_res$li, fac = fac_vector, xax = 1, yax = 2,
                col = rainbow(length(levels(fac_vector))), cstar = 0, grid = TRUE,
                clabel=1.2, sub="Analyse de correspondance (genpop)", cpoint=2)
        legend("topright", legend=levels(fac_vector), fill=rainbow(length(levels(fac_vector))), border=NA, bty="n")
    }} else {{
        writeLines("[AVERTISSEMENT] Moins de deux composantes disponibles pour la CA. Un graphique alternatif sera généré.", log)
        barplot(ca_res$li[,1], col=rainbow(length(ca_res$li[,1])), main="Coordonnées sur la première composante", ylab="Coordonnée CA1")
    }}
    dev.off()

    pdf('{tree_group_output}')
    for (grp in levels(grps)) {{
        sel <- which(grps == grp)
        if (length(sel) > 2) {{
            sub_gen <- genind_obj[sel, ]
            plot(nj(dist(sub_gen)), main=paste("Arbre groupe:", grp))
        }}
    }}
    dev.off()

    pdf(file.path(dirname('{ca_output}'), "dapc_scatterplot_fallback.pdf"))
    temp_pca <- dudi.pca(tab(genind_obj), scannf=FALSE, nf=2)
    fac <- pop(genind_obj)
    s.class(dfxy=temp_pca$li, fac=fac, xax=1, yax=2,
            col=rainbow(length(levels(fac))), cstar=0, grid=TRUE,
            clabel=1.2, cpoint=2,
            xlim=range(temp_pca$li[,1])*1.1, ylim=range(temp_pca$li[,2])*1.1)
    title("Projection PCA des individus selon les groupes", cex.main=1.1, line=1)
    legend("topleft", legend=levels(fac), fill=rainbow(length(levels(fac))), border=NA, bty="n", cex=0.7, inset=c(0.01, 0.01))
    cas_ids <- rownames(tab(genind_obj))[which(fac == "cas")]
    if (length(cas_ids) > 0) {{
        text(temp_pca$li[cas_ids, 1], temp_pca$li[cas_ids, 2], labels=cas_ids,
             cex=0.6, pos=3, col="firebrick")
    }}
    dev.off()

    writeLines("[INFO] Fallback scatterplot PCA généré.", log)
    writeLines("[DEBUG] Fin de l'analyse DAPC", log)
    close(log)
    """

    with open(script_path, "w") as f:
        f.write(script_content)

    logging.info(f"[adegenet] Script R généré : {script_path}")
    return script_path

def run_adegenet_script(script_path: str):
    cmd = f"Rscript {script_path}"
    try:
        logging.info(f"[adegenet] Exécution R : {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        logging.info("[adegenet] Analyse R terminée.")
    except subprocess.CalledProcessError as e:
        logging.error(f"[adegenet] Erreur R : {e}")
        raise
