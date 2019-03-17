library(ggplot2)
library(reshape2)
library(dplyr)

# --------------------------------------------------
# land marker: differential expression at change points for each gene
# --------------------------------------------------

# --------------------------------------------------
# I/O
# --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
infile = args[1]
outprefix = args[2]
pmax = as.numeric(args[3])
dp_dif_min = as.numeric(args[4])

# --------------------------------------------------
# regression: get expected variance
# --------------------------------------------------
df <- read.table(infile, header=T, sep = "\t", as.is=T)

# split mean1 & mean2 into separate lines
df_reg1 = df[,c("seg_dstl", "mean_a", "cv2_a", "nsamples_a")]
df_reg2 = df[,c("seg_dstl", "mean_b", "cv2_b", "nsamples_b")]
colnames(df_reg1) = c("seg", "mean", "cv2", "n")
colnames(df_reg2) = c("seg", "mean", "cv2", "n")
df_reg = unique(rbind(df_reg1, df_reg2))

# loess regression
myloess <- loess(df_reg$cv2 ~ log2(df_reg$mean), data = df_reg)
predY <- predict(myloess, log2(df_reg$mean), se = T)
df_reg$exp_std <- sqrt(df_reg$mean * df_reg$mean * predY$fit)

# plot
pdf(paste(outprefix, "_regression.pdf", sep=""), width=2.5, height=2.5)
print( ggplot(df_reg, aes(log2(mean), cv2)) +
         geom_point(alpha=0.2) + geom_line(aes(y=predY$fit), colour="red") +
         labs(x = "log2(mean)", y="var / mean^2" ) + theme_bw() +
         ggtitle(paste("# change points = ", nrow(df), sep="")) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               aspect.ratio = 1) )
invisible(dev.off())

# --------------------------------------------------
# test
# --------------------------------------------------
df_pred = df_reg %>% group_by(seg) %>%
  summarise(mean1 = min(mean),
            mean2 = max(mean),
            n = n[which.max(mean)],
            exp_std2 = exp_std[which.max(mean)])
df_pred$pval = apply(df_pred, 1, function(x) pnorm(as.numeric(x["mean1"]), as.numeric(x["mean2"]),
                                                  sqrt(as.numeric(x["exp_std2"])**2 / as.numeric(x["n"]))))
df_new = merge(df, df_pred, by.x = "seg_dstl", by.y = "seg")
df_new$pBH = p.adjust(df_new$pval, "BH")

# proximal vs. distal fold change: color p-value
df_new$color = apply(df_new, 1, function(x) {
  this_pbh = as.numeric(as.character(x["pBH"]))
  this_rudif = as.numeric(as.character(x["ru_dif"]))
  return(ifelse(is.na(this_pbh),
                "insignificant",
                ifelse(this_pbh <= pmax,
                       ifelse(is.na(this_rudif),
                              "insignificant",
                              ifelse(abs(this_rudif) >= dp_dif_min,
                                     paste("p<=", pmax, ", abs(RUa-RUb)>=", dp_dif_min, sep=""),
                                     "insignificant")),
                       "insignificant")))
})

pdf(paste(outprefix, "_volcano.pdf", sep=""), width=3, height=4)
print( ggplot(df_new, aes(ru_dif, -log10(pBH))) + geom_point(alpha=1/3, aes(colour=color)) +
         scale_colour_manual(values=c("grey", "blue")) +
         labs(y = "-log10(BH p)", x="RUa - RUb" ) + theme_bw() +
         ggtitle(paste("# change points = ", nrow(df), sep="")) +
         guides(colour=guide_legend(nrow=2)) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               legend.position = "bottom", aspect.ratio = 1) )
invisible(dev.off())

# write output
write.table(df_new[,c("cp", "gene", "test_type", "side", "seg_prxl", "seg_dstl", "ru_dif", "ru_prxla", "ru_prxlb", "ru_dstla", "ru_dstlb", "pval", "pBH")], paste(outprefix, ".txt", sep=""), quote=F, row.names=F, sep="\t")
print(paste(nrow(df_new[!is.na(df_new$pBH) & df_new$pBH <= pmax,]), "/", nrow(df_new), "p <=", pmax, sep=" "))
print(paste(nrow(df_new[!is.na(df_new$pBH) & df_new$pBH <= pmax & !is.na(df_new$ru_dif) & abs(df_new$ru_dif) >= dp_dif_min,]), "/", nrow(df_new), "p <=", pmax, "& abs(RUa - RUb) >=", dp_dif_min, sep=" "))
