#R 4.4.0

library(edgeR)
library(clusterProfiler) 
library(ggplot2)

matrix <- read.table("summary_unique.tsv",sep = "\t",quote = NULL,header = T,
                     row.names = 2,check.names = FALSE)


size <- colSums(matrix)
size

y <- DGEList(counts=matrix)

group <- as.factor(c("KO", "KO", "KO", "KO","KO_treated", "KO_treated", "KO_treated","KO_treated"))
y$samples$group <- group
y$samples

#filtering 
table(rowSums(y$counts==0)==8)
keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)

#pre-normalization 
logcpm_before <- cpm(y, log=TRUE)
boxplot(logcpm_before)

#normalization 
y <- calcNormFactors(y, method = "TMM")

logcpm <- cpm(y, log=TRUE, prior.count = 1)
boxplot(logcpm)

#design matrix 
design <- model.matrix(~group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design 


plotMDS(logcpm, labels=group)

y <- estimateDisp(y, design)
plotBCV(y)
y
y$common.dispersion

#differential expression analysis - KO treated vs KO
fit <- glmQLFit(y, design)
qlf.2vs1 <- glmQLFTest(fit, coef=2) 
topTags(qlf.2vs1)

FDR <- p.adjust(qlf.2vs1$table$PValue, method="BH")

summary(decideTests(qlf.2vs1, p.value=0.05, lfc=1)) 
results <- topTags(qlf.2vs1, n=20000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)$table
results <-(results[which(results$FDR < 0.05 & (results$logFC >1 | results$logFC < -1)), ])
write.table(results, file= 'results_KO_treated.txt',sep='\t')

#enrichment analysis - KO treated vs KO
treated <- rownames(results)

GO_bp_treated <-enrichGO(treated, keyType = "ENSEMBL", 'org.Hs.eg.db', ont="BP",
                       pvalueCutoff=0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH", universe = matrix$gene_id)

dotplot(GO_bp_treated, showCategory=20)

results_GO_treated <- GO_bp_treated@result[which(GO_bp_treated@result$p.adjust < 0.05),]
write.table(results_GO_treated, file="enrichment_treated.txt", sep="\t")

