library(tidyverse)
library(biomaRt)
library(httr)
library(jsonlite)
library(rlang)
# the select function from dplyr is masked by biomaRt, so we need to specify it explicitly
# dplyr::select

# ------------------------------------------------------------
# Load SHaRe variant annotations and gene info
# ------------------------------------------------------------

HCM_Gen_Variant_Annotation_2026Q1 <- read_delim("HCM_Gen_Variant_Annotation_2026Q1.csv", delim = ",")

# identical(HCM_Gen_Variant_Annotation_2026Q1$SHaRe_VarClass...10,
#           HCM_Gen_Variant_Annotation_2026Q1$SHaRe_VarClass...21)
# TRUE

# problems(HCM_Gen_Variant_Annotation_2026Q1)


HCM_gene_info <- read_csv("HCM_gene_info.csv")


# ------------------------------------------------------------
# Get Transcripts and MANE mapping for SHaRe genes
# ------------------------------------------------------------

genes <- HCM_gene_info$Gene
txs <- HCM_gene_info|>separate_longer_delim(TranscriptIDs,delim = "; ")|>rename("ensembl_transcript_id"=TranscriptIDs)

mane <- read_delim(
  "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz",
  delim = "\t"
) |>
  dplyr::select(-1) |>
  separate(Ensembl_Gene,sep = "\\.", into = c("ensembl_gene_id",NA), remove = F) |>
  separate(Ensembl_nuc,sep = "\\.", into = c("ensembl_transcript_id",NA), remove = F)

txs_mane <- txs |>
  left_join(mane, by = "ensembl_transcript_id")|>
  # filter(MANE_status == "MANE Select")|>
  # filter(symbol %in% genes)|>
  arrange(Gene)|>
  dplyr::select(-starts_with("chr"))

txs_mane <- txs_mane |>
  rowwise()|>
  mutate(
    "Ensembl_Gene"   = if_else(ensembl_transcript_id == "ENST00000710581","ENSG00000151067.24", Ensembl_Gene),
    "ensembl_gene_id"= if_else(ensembl_transcript_id == "ENST00000710581","ENSG00000151067", ensembl_gene_id),
    "HGNC_ID"        = if_else(ensembl_transcript_id == "ENST00000710581","HGNC:1390", HGNC_ID),
    "symbol"         = if_else(ensembl_transcript_id == "ENST00000710581","CACNA1C",symbol),
    "name"           = if_else(ensembl_transcript_id == "ENST00000710581","calcium voltage-gated channel subunit alpha1 C",name),
    "RefSeq_nuc"     = if_else(ensembl_transcript_id == "ENST00000710581",NA,RefSeq_nuc),
    "RefSeq_prot"    = if_else(ensembl_transcript_id == "ENST00000710581",NA,RefSeq_prot),
    "Ensembl_nuc"    = if_else(ensembl_transcript_id == "ENST00000710581","ENST00000710581.1",Ensembl_nuc),
    "Ensembl_prot"   = if_else(ensembl_transcript_id == "ENST00000710581","ENSP00000518352.1",Ensembl_prot),
    "MANE_status"    = if_else(ensembl_transcript_id == "ENST00000710581","CANONICAL",MANE_status),
    "GRCh38_chr"     = if_else(ensembl_transcript_id == "ENST00000710581","NC_000012.12",GRCh38_chr)) |>
  mutate(
    "Ensembl_Gene"   = if_else(ensembl_transcript_id == "ENST00000367318","ENSG00000118194.24", Ensembl_Gene),
    "ensembl_gene_id"= if_else(ensembl_transcript_id == "ENST00000367318","ENSG00000118194", ensembl_gene_id),
    "HGNC_ID"        = if_else(ensembl_transcript_id == "ENST00000367318","HGNC:11949", HGNC_ID),
    "symbol"         = if_else(ensembl_transcript_id == "ENST00000367318","TNNT2",symbol),
    "name"           = if_else(ensembl_transcript_id == "ENST00000367318","troponin T2, cardiac type",name),
    "RefSeq_nuc"     = if_else(ensembl_transcript_id == "ENST00000367318",NA,RefSeq_nuc),
    "RefSeq_prot"    = if_else(ensembl_transcript_id == "ENST00000367318",NA,RefSeq_prot),
    "Ensembl_nuc"    = if_else(ensembl_transcript_id == "ENST00000367318","ENST00000367318.10",Ensembl_nuc),
    "Ensembl_prot"   = if_else(ensembl_transcript_id == "ENST00000367318","ENSP00000356287.5",Ensembl_prot),
    "MANE_status"    = if_else(ensembl_transcript_id == "ENST00000367318","Cardiac_Expression",MANE_status),
    "GRCh38_chr"     = if_else(ensembl_transcript_id == "ENST00000367318","NC_000001.11",GRCh38_chr))




txs_ids <- txs_mane$ensembl_transcript_id


# ------------------------------------------------------------
# Connect to Ensembl (gene dataset)
# ------------------------------------------------------------
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

attrs_structure <- c(
  "external_gene_name",
  "ensembl_transcript_id",
  "ensembl_gene_id",
  "transcript_start",
  "transcript_end",
  "start_position",
  "end_position",
  "strand",
  "5_utr_start",
  "5_utr_end",
  "3_utr_start",
  "3_utr_end",
  "cds_length",
  "cds_start",
  "cds_end",
  "ensembl_exon_id",
  "exon_chrom_start",
  "exon_chrom_end"
)

attrs_feature <- c(
  "ensembl_transcript_id",
  "ensembl_gene_id",
  # "external_gene_name",
  "ensembl_peptide_id",
  "refseq_mrna",
  # "refseq_mrna_predicted",
  "pfam",
  "pfam_start",
  "pfam_end"
)



# ------------------------------------------------------------
# Query biomaRt (structure + feature)
# ------------------------------------------------------------

df_structure <- getBM(
  attributes = attrs_structure,
  filters = "ensembl_transcript_id",
  values = txs_ids,
  mart = mart
) 

df_feature <- getBM(
  attributes = attrs_feature,
  filters = "ensembl_transcript_id",
  values = txs_ids,
  mart = mart
) 



# Merge Transcript → structure → feature → MANE
txs_out1 <- txs_mane %>%
  left_join(df_structure, by = "ensembl_transcript_id", relationship = "many-to-many") |>
  left_join(df_feature, by = "ensembl_transcript_id", relationship = "many-to-many")






# ------------------------------------------------------------
# Query biomaRt (Uniprot)
# ------------------------------------------------------------

ensp_ids <- txs_out1$ensembl_peptide_id


df_UniProt <- getBM(
  attributes = c("ensembl_peptide_id", "uniprotswissprot"),
  filters = "ensembl_peptide_id",
  values = ensp_ids,
  mart = mart
)


txs_out2 <- txs_out1 |>
  left_join(df_UniProt, by = "ensembl_peptide_id")
  




# Step 1: Define mapping function
get_genomic_pos <- function(exons, cds_pos, strand) {
  # Order exons depending on strand
  exons <- exons %>%
    arrange(if (strand == 1) exon_chrom_start else desc(exon_chrom_start))
  
  # Compute exon lengths
  exon_lengths <- exons$exon_chrom_end - exons$exon_chrom_start + 1
  cum_lengths <- cumsum(exon_lengths)
  start_offsets <- lag(cum_lengths, default = 0) + 1
  end_offsets <- cum_lengths
  
  # Find exon containing the CDS position
  i <- which(cds_pos >= start_offsets & cds_pos <= end_offsets)[1]
  if (is.na(i)) return(NA)
  
  offset_within_exon <- cds_pos - start_offsets[i]
  if (strand == 1) {
    genomic_pos <- exons$exon_chrom_start[i] + offset_within_exon
  } else {
    genomic_pos <- exons$exon_chrom_end[i] - offset_within_exon
  }
  return(genomic_pos)
}

safe_pfam_bounds <- function(pfam, raw_start, raw_end, bound = c("start", "end")) {
  bound <- match.arg(bound)

  if (is.na(pfam) || pfam == "" || all(is.na(c(raw_start, raw_end)))) {
    return(NA_real_)
  }

  vals <- c(raw_start, raw_end)
  vals <- vals[!is.na(vals)]

  if (length(vals) == 0) {
    return(NA_real_)
  }

  if (bound == "start") min(vals) else max(vals)
}

# Step 2: Apply per transcript
share_gene_info <- txs_out2 %>%
  mutate(
    cds_pos_start = (pfam_start - 1) * 3 + 1,
    cds_pos_end   = pfam_end * 3
  ) %>% 
  group_by(ensembl_transcript_id) %>%
  group_modify(~ {
    strand <- .x$strand[1]
    
    # Get all exons for this transcript
    exons <- .x %>%
      dplyr::select(exon_chrom_start, exon_chrom_end) %>%
      distinct()
    
    .x %>%
      dplyr::rowwise() %>%
      mutate(
        raw_start = get_genomic_pos(exons, cds_pos_start, strand),
        raw_end   = get_genomic_pos(exons, cds_pos_end, strand),
        pfam_start_g = safe_pfam_bounds(pfam, raw_start, raw_end, "start"),
        pfam_end_g   = safe_pfam_bounds(pfam, raw_start, raw_end, "end")
      ) %>%
      dplyr::select(-raw_start, -raw_end, -starts_with("cds")) %>%
      ungroup()
  }) %>%
  ungroup() %>%
  rename(
    transcript_stable_id_version = Ensembl_nuc,
    refseq_match_transcript_mane_select = RefSeq_nuc,
    gene_stable_id_version = Ensembl_Gene,
    exon_region_start_bp = exon_chrom_start,
    exon_region_end_bp = exon_chrom_end,
    pfam_id = pfam
  ) %>% 
  mutate(hg38_chr = dplyr::case_when(
    GRCh38_chr == "NC_000023.11" ~ "chrX",
    GRCh38_chr == "NC_000024.10" ~ "chrY",
    GRCh38_chr == "NC_012920.1" ~ "chrM",
    TRUE ~ paste0("chr", as.integer(sub("NC_0+([0-9]+)\\..*", "\\1", GRCh38_chr)))
  ), .after = GRCh38_chr)






# Step 3: Get Pfam name and description


get_pfam_name <- function(df, pfam_col) {
  
  pf_col <- ensym(pfam_col)
  
  # Unique Pfam IDs
  unique_ids <- df %>% distinct(!!pf_col) %>% pull(!!pf_col)
  unique_ids <- unique_ids[!is.na(unique_ids)]
  
  fetch_one <- function(pf) {
    url <- paste0(
      "https://www.ebi.ac.uk/interpro/api/entry/pfam/",
      pf
    )
    
    res <- GET(url, accept("application/json"))
    
    if (res$status_code != 200) {
      return(tibble(pfam_id = pf, pfam_name = NA_character_))
    }
    
    x <- fromJSON(rawToChar(res$content))
    
    # Extract name safely
    name_raw <- x$metadata$name
    
    # Force to character scalar
    name_clean <- if (is.null(name_raw)) {
      NA_character_
    } else if (is.list(name_raw)) {
      as.character(name_raw[[1]])
    } else {
      as.character(name_raw)
    }
    
    tibble(
      pfam_id = pf,
      pfam_name = name_clean
    )
  }
  
  lookup <- map_dfr(unique_ids, fetch_one)
  
  df %>% left_join(lookup, by = setNames("pfam_id", as_string(pf_col)))
}





share_genes <- share_gene_info %>%
  get_pfam_name(pfam_id) |>
  relocate("uniprot_id" = uniprotswissprot,.after = last_col())




# ------------------------------------------------------------
# Write HGVSc inputs for local VEP annotations of SHaRe variants
# ------------------------------------------------------------

# SHaRe annotations
share_2026_Q1 <- HCM_Gen_Variant_Annotation_2026Q1 |>
  mutate(VarClass_date = if_else(is.na(VarClass_date), as_date("2025-06-04"), VarClass_date))|>
  dplyr::select(-c(10,13,16,19,21,35,36))|>
  filter(VarClass_status != "nonHCM gene")|>
  filter(VarClass != "nonHCM VUS")|>
  filter(!str_starts(VariantID, "H")) |>
  filter(!str_starts(external_gene_name,"Intronic")) |>
  rename(
    Impact = VEP_impact
  )


vep_hgvsc_input <- share_2026_Q1 %>%
  filter(!is.na(HGVSc), HGVSc != "") %>%
  separate(VariantID, into = c("chromosome", "position", "ref", "alt"), sep = "-", remove = FALSE, fill = "right") %>%
  # mutate(HGVSc = if_else(VariantID=="11-47331882-CT-<GSPLICEACCEPTOR>", "ENST00000545968.6:c.3815-2_3815-1del", HGVSc)) %>% 
  # mutate(HGVSc = if_else(VariantID=="11-47351506-CT-<GSPLICEACCEPTOR>", "ENST00000545968.6:c.26-3_26-2del", HGVSc)) %>% 
  mutate(
    chromosome = str_remove(chromosome, "^chr"),
    chromosome = case_when(
      chromosome %in% as.character(1:22) ~ chromosome,
      str_to_upper(chromosome) == "X" ~ "X",
      str_to_upper(chromosome) == "Y" ~ "Y",
      TRUE ~ chromosome
    ),
    chromosome_rank = match(chromosome, c(as.character(1:22), "X", "Y")),
    chromosome_rank = coalesce(chromosome_rank, 99L),
    position = suppressWarnings(as.numeric(position))
  ) %>%
  distinct(VariantID, HGVSc, .keep_all = TRUE) %>%
  arrange(chromosome_rank, position, HGVSc) %>%
  dplyr::select(VariantID, HGVSc)


  

# 11:47351506-47351508 - ENST00000545968.6:c.26-3_26-2del
# 11:47331881-47331883 - ENST00000545968.6:c.3815-2_3815-1del

# ------------------------------------------------------------
# Write VEP input lines in HGVSc
# ------------------------------------------------------------

# write_lines(vep_hgvsc_input$HGVSc, "../data/share_2026Q1_vep_hgvsc.txt")
# write_csv(vep_hgvsc_input, "../data/share_2026Q1_vep_hgvsc_map.csv")


# ------------------------------------------------------------
# Read in VEP115
# ------------------------------------------------------------

share_2026Q1_vep_hgvsc_map <- read_csv("../data/share_2026Q1_vep_hgvsc_map.csv")

share_2026Q1_VEP115_in <- read_delim("../data/share_2026Q1_vep_hgvsc.VEP115.txt", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)



share_2026Q1 <- share_2026Q1_vep_hgvsc_map |>
  left_join(share_2026Q1_VEP115_in, by = "HGVSc")

which(duplicated(share_2026Q1$VariantID))

share_2026Q1_uniq <- distinct(share_2026Q1) |>
  filter(!is.na(CHROM)) |>
  filter(!duplicated(HGVSc))

which(duplicated(share_2026Q1_uniq$HGVSc))

share_2026Q1_VEPout <- share_2026Q1_uniq |>
  left_join(share_2026_Q1) |>
  filter(!is.na(VarClass))







# ------------------------------------------------------------
# Prepare app-ready objects for the SHaRe browser
# ------------------------------------------------------------

genomic_classifications <- c("P/LP", "VUS", "B/LB")
size_map <- c("0" = 1, "1-5" = 10, "6-10" = 15, "11-20" = 20, "21-50" = 30, ">50" = 40)

bin_AC <- function(AC) {
  AC <- coalesce(AC, 0)
  cut(
    AC,
    breaks = c(0, 1, 5, 10, 20, 50, 200),
    labels = c("0", "1-5", "6-10", "11-20", "21-50", ">50"),
    right = TRUE,
    include.lowest = TRUE
  )
}

normalize_variant_data <- function(df) {
  df <- tibble::as_tibble(df)
  
  if (!"Clinvar_LastEvaluated" %in% names(df)) df$Clinvar_LastEvaluated <- NA_character_
  if (!"Impact" %in% names(df)) df$Impact <- NA_character_
  if (!"VEP_consequence" %in% names(df)) df$VEP_consequence <- NA_character_
  if (!"SpliceAI_GENE" %in% names(df)) df$SpliceAI_GENE <- NA_character_
  
  df %>%
    mutate(
      Gene = dplyr::coalesce(.data$SYMBOL, .data$Gene),
      Impact = dplyr::coalesce(.data$Impact, .data$IMPACT),
      VEP_consequence = dplyr::coalesce(.data$VEP_consequence, .data$Consequence),
      SpliceAI_GENE = dplyr::coalesce(.data$SpliceAI_GENE, .data$SpliceAI_pred_SYMBOL),
      SpliceAI_max = pmax(
        as.numeric(.data$SpliceAI_pred_DS_AG),
        as.numeric(.data$SpliceAI_pred_DS_AL),
        as.numeric(.data$SpliceAI_pred_DS_DG),
        as.numeric(.data$SpliceAI_pred_DS_DL),
        na.rm = TRUE
      ),
      SpliceAI_max = na_if(SpliceAI_max, -Inf),
      SpliceAI_pred = SpliceAI_max,
      Flag_PTV = str_detect(
        str_to_lower(coalesce(.data$Consequence, "")),
        "stop_gained|frameshift|splice_acceptor|splice_donor"
      ),
      Flag_PAV = str_detect(
        str_to_lower(coalesce(.data$Consequence, "")),
        "missense|inframe_insertion|inframe_deletion|start_lost|stop_lost"
      )
    ) %>%
    separate(VariantID, into = c("CHR", "POS"), extra = "drop", remove = FALSE, fill = "right") %>%
    mutate(POS = suppressWarnings(as.numeric(POS))) %>%
    filter(!is.na(POS), !is.na(HGVSc), !is.na(AC)) %>%
    distinct(VariantID, .keep_all = TRUE) %>%
    mutate(
      VariantID = as.character(VariantID),
      classification = factor(VarClass, levels = genomic_classifications),
      AC_bin = bin_AC(AC),
      ypos = as.numeric(classification),
      size = unname(size_map[as.character(AC_bin)])
    ) %>%
    separate(HGVSc, into = c("transcript_ID", "cdot"), sep = ":", remove = FALSE, fill = "right") %>%
    separate(HGVSp, into = c("protein_ID", "pdot"), sep = ":", remove = FALSE, fill = "right") %>%
    dplyr::select(-any_of(c("ensembl_transcript_id", "external_gene_name")))
}

normalize_gene_data <- function(df) {
  tibble::as_tibble(df) %>%
    mutate(
      gene_name = Gene,
      exon_stable_id = ensembl_exon_id,
      gene_start_bp = start_position,
      gene_end_bp = end_position
    )
}

extract_exons <- function(df) {
  df %>%
    dplyr::select(-contains(c("pfam", "cds"))) %>%
    distinct(exon_stable_id, .keep_all = TRUE)
}

extract_domains <- function(df) {
  df %>%
    dplyr::select(contains(c("gene_name", "pfam"))) %>%
    distinct(pfam_name, .keep_all = TRUE) %>%
    arrange(pfam_start_g)
}

extract_meta <- function(df) {
  df %>%
    dplyr::select(gene_start_bp, gene_end_bp, strand) %>%
    distinct() %>%
    slice_head(n = 1)
}

share <- normalize_variant_data(share_2026Q1_VEPout)
gene_info_all <- normalize_gene_data(share_genes)
gene_info_filtered <- gene_info_all %>%
  filter(MANE_status == "MANE Select")
hcm_genes <- sort(unique(gene_info_all$Gene[!is.na(gene_info_all$Gene)]))
multi_tx_genes <- c("CACNA1C", "FHL1", "TNNT2")
multi_tx_exclude <- list(
  CACNA1C = c("ENST00000710581.1")
)

share_by_gene <- split(share, share$Gene)

gene_structures <- lapply(split(gene_info_filtered, gene_info_filtered$gene_name), function(df) {
  list(
    exons = extract_exons(df),
    domains = extract_domains(df),
    meta = extract_meta(df)
  )
})

gene_regions <- gene_info_filtered %>%
  group_by(gene_name, hg38_chr) %>%
  summarise(
    gene_start_bp = min(gene_start_bp, na.rm = TRUE),
    gene_end_bp = max(gene_end_bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(gene_start_bp), is.finite(gene_end_bp))

empty_exons <- gene_info_filtered[0, ] %>% extract_exons()
empty_domains <- gene_info_filtered[0, ] %>% extract_domains()
empty_meta <- gene_info_filtered[0, ] %>% extract_meta()

gene_clinvar_counts <- lapply(share_by_gene, function(df) {
  df %>%
    count(VarClass) %>%
    mutate(VarClass = factor(VarClass, levels = genomic_classifications)) %>%
    arrange(desc(VarClass))
})








# ------------------------------------------------------------
# Save objects to .RData file
# ------------------------------------------------------------

save(
  share,
  gene_info_filtered,
  gene_info_all,
  hcm_genes,
  multi_tx_genes,
  multi_tx_exclude,
  share_by_gene,
  gene_structures,
  gene_regions,
  empty_exons,
  empty_domains,
  empty_meta,
  gene_clinvar_counts,
  genomic_classifications,
  size_map,
  file = "share_info_2026Q1.RData"
)


##################
# Later, load them back into your environment
# load("share_info_2026Q1.RData")
