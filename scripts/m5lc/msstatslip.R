# MSstatsLiP 1.8.1 on the Vulcan D7 export (scripts/m5lc/lip_quant.py): F1 (alpha-synuclein fibrils) vs M1 (monomer).
# Rscript scripts/m5lc/msstatslip.R OUTDIR FASTA
# Inputs OUTDIR/msstatslip_{LiP,TrP}.csv (MSstats format; LiP adds FULL_PEPTIDE). Defaults of the package, except
# MBimpute.LiP = FALSE (package default) and equalizeMedians normalisation (package default).
suppressMessages({library(MSstatsLiP); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]; fasta <- args[2]
t0 <- Sys.time()
lip <- fread(file.path(out, "msstatslip_LiP.csv"))
trp <- fread(file.path(out, "msstatslip_TrP.csv"))
lip$Intensity <- as.numeric(lip$Intensity); trp$Intensity <- as.numeric(trp$Intensity)
s <- dataSummarizationLiP(list(LiP = lip, TrP = trp), use_log_file = FALSE, verbose = FALSE)
t1 <- Sys.time()
m <- groupComparisonLiP(s, fasta.path = fasta, use_log_file = FALSE, verbose = FALSE)
t2 <- Sys.time()
fwrite(as.data.table(m$LiP.Model), file.path(out, "msstatslip_LiP_model.csv"))
fwrite(as.data.table(m$TrP.Model), file.path(out, "msstatslip_TrP_model.csv"))
if (!is.null(m$Adjusted.LiP.Model)) fwrite(as.data.table(m$Adjusted.LiP.Model), file.path(out, "msstatslip_Adjusted_LiP_model.csv"))
writeLines(jsonlite::toJSON(list(summarization_s = as.numeric(difftime(t1, t0, units = "secs")),
                                 comparison_s = as.numeric(difftime(t2, t1, units = "secs")),
                                 msstatslip = as.character(packageVersion("MSstatsLiP")),
                                 msstats = as.character(packageVersion("MSstats"))), auto_unbox = TRUE),
           file.path(out, "msstatslip_run.json"))
cat("done\n")
