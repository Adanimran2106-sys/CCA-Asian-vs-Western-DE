library(DESeq2); library(limma)

# ===== ASIAN (DESeq2) =====
cnt <- as.matrix(asian[, 8:ncol(asian)])
ord <- order(rowSums(cnt), decreasing = TRUE)
cnt <- cnt[ord, ]; sym <- asian$Genesymbol[ord]
keep_g <- !duplicated(sym)
cnt <- cnt[keep_g, ]; rownames(cnt) <- sym[keep_g]

patient <- sub("[NT]$", "", colnames(cnt))
type    <- sub("^.*([NT])$", "\\1", colnames(cnt))
both    <- intersect(patient[type == "N"], patient[type == "T"])
use     <- both[1:18]
k       <- patient %in% use
cnt <- cnt[, k]; patient <- patient[k]; type <- type[k]

coldata <- data.frame(patient = factor(patient),
                      type = factor(type, levels = c("N", "T")))
dds <- DESeqDataSetFromMatrix(round(cnt), coldata, ~ patient + type)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res_a <- as.data.frame(results(dds, contrast = c("type", "T", "N")))
res_a$gene <- rownames(res_a)
deg_a <- subset(res_a, padj < 0.05 & abs(log2FoldChange) > 1)

# ===== WESTERN (limma) =====
tp <- substr(colnames(west), 14, 15)
pw <- substr(colnames(west), 1, 12)
sel <- tp == "11" | (tp == "01" & pw %in% pw[tp == "11"])
ex  <- as.matrix(west[, sel])
grp <- factor(ifelse(tp[sel] == "11", "N", "T"), levels = c("N", "T"))
pat <- factor(pw[sel])
print(table(grp))
design <- model.matrix(~ pat + grp)
fit <- eBayes(lmFit(ex, design))
res_w <- topTable(fit, coef = "grpT", number = Inf)
res_w$gene <- rownames(res_w)
deg_w <- subset(res_w, adj.P.Val < 0.05 & abs(logFC) > 1)

# ===== COMPARISON =====
common <- intersect(deg_a$gene, deg_w$gene)
cat("Asian DEGs:", nrow(deg_a), "\nWestern DEGs:", nrow(deg_w),
    "\nCommon:", length(common), "\n")
m <- merge(deg_a[, c("gene", "log2FoldChange")],
           deg_w[, c("gene", "logFC")], by = "gene")
m$same_direction <- sign(m$log2FoldChange) == sign(m$logFC)
head(m)
write.csv(m, "common_DEGs.csv", row.names = FALSE)

#  1) Volcano plot: Asian 
png("volcano_asian.png", width = 1600, height = 1400, res = 200)
plot(res_a$log2FoldChange, -log10(res_a$padj), pch = 20, col = "grey70",
     xlab = "log2 Fold Change", ylab = "-log10 adj. p-value",
     main = "Asian: Tumor vs Normal")
sig <- res_a$padj < 0.05 & abs(res_a$log2FoldChange) > 1
points(res_a$log2FoldChange[sig], -log10(res_a$padj[sig]), pch = 20, col = "red")
abline(v = c(-1, 1), h = -log10(0.05), lty = 2)
dev.off()

#  2) Volcano plot: Western 
png("volcano_western.png", width = 1600, height = 1400, res = 200)
plot(res_w$logFC, -log10(res_w$adj.P.Val), pch = 20, col = "grey70",
     xlab = "log2 Fold Change", ylab = "-log10 adj. p-value",
     main = "Western: Tumor vs Normal")
sigw <- res_w$adj.P.Val < 0.05 & abs(res_w$logFC) > 1
points(res_w$logFC[sigw], -log10(res_w$adj.P.Val[sigw]), pch = 20, col = "red")
abline(v = c(-1, 1), h = -log10(0.05), lty = 2)
dev.off()

# 3) Venn diagram (base R) 
a <- nrow(deg_a); w <- nrow(deg_w); cm <- length(common)
png("venn.png", width = 1600, height = 1200, res = 200)
plot(NA, xlim = c(0, 10), ylim = c(0, 6), asp = 1, axes = FALSE,
     xlab = "", ylab = "", main = "Overlap of DEGs")
symbols(c(3.8, 6.2), c(3, 3), circles = c(2.2, 2.2), inches = FALSE, add = TRUE,
        bg = c(rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.3)))
text(2.7, 3, a - cm); text(7.3, 3, w - cm); text(5, 3, cm)
text(2.7, 5.6, "Asian"); text(7.3, 5.6, "Western")
dev.off()

getwd()   # plots isi folder mein save hue hain