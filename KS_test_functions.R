
# condsider changing min_p to .Machine$double.xmin
run.ks.test <- function(df, group_var, rank_var, ties_method = "random", alt = "two.sided", 
	flip_rnk = T, min_p = 2.2e-16){
	ks_res <- data.frame(matrix(nrow = length(unique(df[[group_var]])), ncol = 4))
	colnames(ks_res) <- c("stat", "p.value", "D+", "D-")
	if(flip_rnk){
		rnk <- rank(-df[[rank_var]], ties = ties_method)
	} else {
		rnk <- rank(df[[rank_var]], ties = ties_method)
	}
	for(i in 1:length(unique(df[[group_var]]))){
		sub_rnk <- rnk[df[[group_var]] == unique(df[[group_var]])[i]]
		ks_g <- ks.test(sub_rnk, "punif", 1, nrow(df), alternative = "greater")
		ks_l <- ks.test(sub_rnk, "punif", 1, nrow(df), alternative = "less")
		if(alt == "two.sided"){
			ks_two <- ks.test(sub_rnk, "punif", 1, nrow(df), alternative = "two.sided")
			ks_res[i,1:2] <- c(ks_two$statistic, ks_two$p.value)
		} else if(alt %in% c("g", "greater")){
			ks_res[i,1:2] <- c(ks_g$statistic, ks_g$p.value)
		} else if(alt %in% c("l", "less")){
			ks_res[i,1:2] <- c(ks_l$statistic, ks_l$p.value)
		}
		ks_res[i,3] <- ks_g$stat
		ks_res[i,4] <- ks_l$stat
	}
	rownames(ks_res) <- unique(df[[group_var]])
	rnk <- rank(-df[[rank_var]], ties = "random")[df[[group_var]] == "TAG"]
	ks.test(rnk, "punif", 1, 904, alternative = "two.sided")

	ks_res[ks_res == 0] <- min_p
	if(alt == "two.sided"){
		ks_res$signed.p.score <- -log10(ks_res$p.value) * ifelse(ks_res[["D+"]] < ks_res[["D-"]], -1, 1)
	}
	return(ks_res)
}

source("C:/Users/Nick/Desktop/GBM/scripts/FGSEA_functions.R", chdir = T)

ks.test.and.plot <- function(cor_df, prefix){
	mean_df <- data.frame(cor_df %>% 
	group_by(class) %>%
	summarize(mean = mean(p.score)))

	cor_df$class <- factor(cor_df$class, levels = mean_df[order(mean_df[,2], decreasing = T),1])

	gg <- ggplot(cor_df, aes(x = rank(-p.score, ties.method = "random"), y = p.score)) + geom_bar(aes(color = sig), size = 0.25, width = 0.25, stat = "identity") + 
	facet_wrap(~ class, ncol = 1) + geom_hline(yintercept = 0, size = 0.25) + scale_color_manual(values = c("black", "red")) + theme_classic() +
	theme(strip.background = element_blank(), strip.text.x = element_blank(), 
		axis.line.x = element_blank(), axis.text.y = element_blank(), legend.position = "none",
		axis.title = element_text(size = 9), axis.text.x = element_text(size = 8, color = "black"),
		axis.line.y = element_line(size = 0.33), axis.ticks = element_line(size = 0.33)) +
	labs(x = "Signed p-value rank", y = "Signed log10 p-value") + scale_x_continuous(expand = c(0.01,0.01))

	ggsave(paste0("Lipid Correlation w ", prefix, " by Class.png"), gg, type = "cairo", dpi = 300, height = 4, width = 2)

	ks_res <- run.ks.test(cor_df, "class", "p.score")
	ks_res$class <- rownames(ks_res)
	ks_res$class <- factor(ks_res$class, levels = mean_df[order(mean_df[,2], decreasing = F),1])
	gg <- ggplot(ks_res, aes(x = signed.p.score, y = class, fill = signed.p.score > 0)) + 
	geom_bar(stat = "identity", color = "black", width = 0.75) + 
	scale_fill_manual(limits = c(TRUE, FALSE), values = c("red", "blue")) + 
	theme_classic() + geom_vline(xintercept = 0) + geom_vline(xintercept = c(log10(0.05), -log10(0.05)), lty = 2, lwd = 0.25) +
	theme(legend.position = "none", axis.line.y = element_blank(), axis.text.y = element_text(size = 8, color = "black"),
			axis.text.x = element_text(size = 8, color = "black"), axis.ticks.y = element_blank(), axis.title.x = element_text(size = 9)) + 
	labs(x = "Signed log10 p-value", y= "", fill = "")

	ggsave(paste0(prefix, " - lipid class KS test barchart.png"), gg, type = "cairo", dpi = 300, height = 2, width = 2.5)

	rnk <- cor_df$p.score
	names(rnk) <- rownames(cor_df)
	lipid_sets <- split(rownames(cor_df), cor_df$class)
	gsea_res <- data.frame(run.FGSEA(rnk, genesets = lipid_sets, minP = 1e-30))
	gsea_res$signed.p.score <- -log10(gsea_res$pval) * sign(gsea_res$NES)
	gsea_res$pathway <- factor(gsea_res$pathway, levels = mean_df[order(mean_df[,2], decreasing = F),1])
	gg <- ggplot(gsea_res, aes(x = signed.p.score, y = pathway, fill = signed.p.score > 0)) + 
	geom_bar(stat = "identity", color = "black", width = 0.75) + 
	scale_fill_manual(limits = c(TRUE, FALSE), values = c("red", "blue")) + 
	theme_classic() + geom_vline(xintercept = 0) + geom_vline(xintercept = c(log10(0.05), -log10(0.05)), lty = 2, lwd = 0.25) +
	theme(legend.position = "none", axis.line.y = element_blank(), axis.text.y = element_text(size = 8, color = "black"),
		axis.text.x = element_text(size = 8, color = "black"), axis.ticks.y = element_blank(), axis.title.x = element_text(size = 9)) + 
	labs(x = "Signed log10 p-value", y= "", fill = "")

	ggsave(paste0(prefix, " - lipid class GSEA test barchart.png"), gg, type = "cairo", dpi = 300, height = 2, width = 2.5)
	return(list(mean = mean_df, ks = ks_res, gsea = gsea_res))
}
