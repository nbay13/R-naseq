get.group.labels <- function(x){
	group_labels <- gsub("[\\.]*[0-9].*", "", x)
	return(group_labels)
}

get.id.labels <- function(x){
	temp <- gsub("[A-Z]*", "", x)
	id <- gsub("\\..*", "",temp)
	final_ids <- substr(id, 1, 3)
	return(final_ids)
}

run.generalized.anova <- function(data, vars){
	vars <- apply(vars, 2, function(x){
		if(!is.factor(x)){
			factor(x)
		}
	})
	res_list <- apply(data, 1, function(x){
		temp <- data.frame(vals = x, vars)
		design <- as.formula(paste("vals", paste(colnames(vars), collapse = " + "), sep = " ~ "))
		aov(data = temp, design)
	})
	names(res_list) <- rownames(data)
	return(res_list)
}

run.two.way.anova <- function(data, between_var, within_var){
	if(!is.factor(between_var)){
		between_var <- factor(between_var)
	}
	if(!is.factor(within_var)){
		within_var <- factor(within_var)
	}
	res_list <- apply(data, 1, function(x){
		aov(as.numeric(x) ~ between_var + within_var)
	})
	names(res_list) <- rownames(data)
	return(res_list)
}

get.summary <- function(res_list, var_name){
	summary(res_list[[var_name]])
}

get.post.hoc.summary <- function(res_list, var_name){
	TukeyHSD(res_list[[var_name]])
}

get.p.values <- function(res_list, n_vars){
	p_vals <- lapply(res_list, function(x){
		summary(x)[[1]][["Pr(>F)"]][1:n_vars]
		})
	df <- do.call(rbind, p_vals)
	colnames(df) <- c("Between.Var", "Within.Var")
	return(df)
}

add.var.names <- function(res_df, between_var_name, within_var_name){
	colnames(res_df) <- c(between_var_name, within_var_name)
	return(res_df)
}

adjust.p.values <- function(res_df, adj_method = "fdr"){
	adj <- apply(res_df, 2, function(x){
		p.adjust(x, method = adj_method)
	})
	colnames(adj) <- paste0("adj.",colnames(res_df))
	new_df <- data.frame(cbind(res_df, adj))
	return(new_df)
}

anova.from.file <- function(data_filename, anno_filename, output_filename){
	data_tabl <- read.table(data_filename, sep = "\t", stringsAsFactors = F, row.names = 1, header = T)
	anno_tabl <- read.table(anno_filename, sep = "\t", stringsAsFactors = F, header = T)
	final_tabl <- na.omit(data_tabl)
	anova_list <- run.generalized.anova(final_tabl, anno_tabl)
	pvalue_df <- get.p.values(anova_list, ncol(anno_tabl))
	colnames(pvalue_df) <- colnames(anno_tabl)
	final_df <- adjust.p.values(pvalue_df)
	write.table(final_df, output_filename, sep = "\t", row.names = T, col.names = NA, quote = F)
	cat(paste0("\nOutput file: '", output_filename, "'\nSaved at: ", getwd(), "\n"))
	return(anova_list)
}
