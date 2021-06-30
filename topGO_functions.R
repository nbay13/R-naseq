
list2genes <- function(sig_list, background_list){
	background_list <- setdiff(background_list, sig_list)
	genes <- c(rep(1, length(sig_list)), rep(0, length(background_list)))
	names(genes) <- c(sig_list, background_list)
	return(genes)
}

code2genes <- function(DE_df, code){
	x <- T
	for(i in 1:length(code)){
		x <- DE_df[,i] == code[i] & x
	}
	sig <- rownames(DE_df)[x]
	#sig <- rownames(DE_df)[DE_df[,1] == code[1] & DE_df[,2] == code[2] & DE_df[,3] == code[3]]
	set <- setdiff(rownames(DE_df), sig)
	genes <- c(rep(1, length(sig)), rep(0, length(set)))
	names(genes) <- c(sig, set)
	return(genes)
}

run.topGO <- function(geneList, method = "classic", ontology = "BP", size = 1000){
	library(topGO)
	library(org.Hs.eg.db)
	## prepare data
	tgd <- new("topGOdata", ontology = ontology, allGenes = factor(geneList), nodeSize=2,
		annot=annFUN.org, mapping="org.Hs.eg.db", ID = "alias")
	## run tests
	if(method == "classic"){
		tst <- runTest(tgd, algorithm = "classic", statistic = "Fisher")
		tbl <- GenTable(tgd, classic = tst, orderBy = "classic", topNodes = size, numChar = 1000)
	} else if(method == "weight01"){
		tst <- runTest(tgd, algorithm = "weight01", statistic = "Fisher" )
		tbl <- GenTable(tgd, weight01 = tst, orderBy = "weight01", topNodes = size, numChar = 1000)
	} else if(method == "parentchild"){
		tst <- runTest(tgd, algorithm = "parentchild", statistic = "Fisher")
		tbl <- GenTable(tgd, parentchild = tst, orderBy = "parentchild", topNodes = size, numChar = 1000)
		} else {
		return("no matching method")
	}
	## create table
	tbl$pval <- score(tst)[tbl$GO.ID]
	return(list(tgd = tgd, tbl = tbl))
}

search.tbl <- function(tbl, searchList, p.value = 0.05){
	inSearch <- grepl(paste(searchList,collapse="|"), tbl$Term)
	sig <- tbl[inSearch & tbl$pval < p.value,]
	return(sig)
}

GO2gene.signature <- function(tgd, GOList, geneList){
	mygenes <- genesInTerm(tgd, GOList)
	temp <- unique(unlist(mygenes))
	final <- temp[temp %in% names(geneList)[geneList == 1]]
	return(final)
}

GO2DE.metrics <- function(tgd, tbl, geneList, DETable, absolute = F, metric = "log2FoldChange"){
	fold_changes <- vector(length = nrow(tbl))
	sig <- names(geneList)[geneList == 1]
	for(i in 1:nrow(tbl)){
		mygenes <- genesInTerm(tgd, tbl$GO.ID[i])
		temp <- unique(unlist(mygenes))
		final <- temp[temp %in% sig]
		if(metric == "signed logp"){
			vals <- -log10(DETable[["padj"]][rownames(DETable) %in% final]) * sign(DETable[["stat"]][rownames(DETable) %in% final])
		} else {
			vals <- DETable[[metric]][rownames(DETable) %in% final]
		}
		if(absolute){
			val <- mean(abs(vals))
		} else {
			val <- mean(vals)
		}
		fold_changes[i] <- val
	}
	return(fold_changes)
}

## from plotGODESeq.R on github
## Core wrapping function
wrap.it <- function(x, len){ 
	sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
## Call this function with a list or vector
wrap.labels <- function(x, len){
	if (is.list(x)){
		lapply(x, wrap.it, len)
	} else {
		wrap.it(x, len)
	}
}

centered.breaks <- function(pal, vals){
	len <- length(pal)
	breaks <- c(
		seq(min(vals), 0, length.out=ceiling(len/2) + 1),
		seq(max(vals)/len, max(vals), length.out=floor(len/2))
		)
	return(breaks)
}

linMap <- function(x){(x-min(x))/(max(x)-min(x))}

bubbleGO.plot <- function(tbl, GOList, wrap = 25, x_max = 3.5, pval_min = 5*10^-10, x_label = expression('Avg. gene log'[2]~'Fold Change'),
	x_thresh = 1, pval_thresh = 0.005, metric = "l2fc", force = 2, xlim = c(-3.5,3.5), ylim = c(0,10), legend_label = expression('Avg. gene\nlog'[2]~'Fold Change')){
	library(ggplot2)
	library(ggrepel)
	library(RColorBrewer)
	color_pal <- brewer.pal("RdBu", n = 9)
	tbl[[metric]][is.na(tbl[[metric]])] <- 0
	tbl[[metric]][abs(tbl[[metric]]) > x_max] <- x_max * sign(tbl[[metric]][abs(tbl[[metric]]) > x_max])
	tbl$pval[tbl$pval < pval_min] <- pval_min
	sig_tbl <- tbl[tbl$GO.ID %in% GOList,]
	back_tbl <- tbl[!tbl$GO.ID %in% GOList,]
	color_breaks <- centered.breaks(color_pal, sig_tbl[[metric]])
	label_tbl <- sig_tbl[abs(sig_tbl[[metric]]) > x_thresh & sig_tbl$pval < pval_thresh,]
	labels <- wrap.labels(label_tbl$Term, wrap)
	gg_bubble <- ggplot() + 
	geom_point(aes(x = as.numeric(back_tbl[[metric]]), 
		y = -log10(as.numeric(back_tbl$pval)),
		size = as.numeric(back_tbl$Significant)),
		fill = "grey90", col = "grey90", shape = 21) +
	geom_hline(yintercept = -log10(0.05), lty = 2) +
	geom_point(aes(x = as.numeric(sig_tbl[[metric]]), 
		y = -log10(as.numeric(sig_tbl$pval)), 
		size = as.numeric(sig_tbl$Significant),
		fill = as.numeric(sig_tbl[[metric]])), col = "white", shape = 21) +
	scale_size_binned(range=c(2,8), breaks = c(10,25,50,100)) + ylim(0,-log10(min(tbl$pval))) + xlim(-max(abs(tbl[[metric]])),max(abs(tbl[[metric]]))) + theme_classic() +
	geom_vline(xintercept = 0, lty = 2) +
	geom_text_repel(aes(x = as.numeric(label_tbl[[metric]]), 
		y = -log10(as.numeric(label_tbl$pval)),
		label = labels), force = force, xlim = xlim, ylim = ylim) +
	scale_fill_gradientn(colors = color_pal, limits = c(min(sig_tbl[[metric]]),max(sig_tbl[[metric]])), 
		breaks = color_breaks, labels = c(round(color_breaks[1],1), rep("", length(color_breaks)-2), round(color_breaks[length(color_breaks)],1)), values = linMap(color_breaks)) +
	labs(x = x_label, y = expression('Geneset -log'[10]~'p-value'), fill = legend_label, size = "No. of genes") +
	theme(legend.position = "right", axis.text = element_text(color = "black"))
	plot(gg_bubble)
	return(gg_bubble)
}


list.to.gmx <- function(list, wd, filename){
	old_dir <- getwd()
	n.obs <- sapply(list, length)
	seq.max <- seq_len(max(n.obs))
	mat <- sapply(list, "[", i = seq.max)
	setwd(wd)
	write.table(mat, filename, col.names = T, row.names = F, quote = F, sep = "\t")
	setwd(old_dir)
}