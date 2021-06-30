
ssgsea.norm <- function(x, min, max){
	return(x / (max - min))
}

#linMap <- function(x){(x-min(x))/(max(x)-min(x))}
#qinMap <- function(x){(x-quantile(x,0.05))/(quantile(x, 0.95)-quantile(x, 0.05))}
qinProj <- function(x,lo,hi){(x-lo)/(hi-lo)}


call.cellular.state <- function(state_comp){
	call <- apply(state_comp, 1, function(x){
		names(sort(x, decreasing = T))[1:2]
		})
	return(t(call))
}

eucDist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

get.simplicity.score <- function(state_comp){
	scores <- apply(state_comp, 1, function(x){
		eucDist(x, c(0.25,0.25,0.25,0.25)) / eucDist(c(1,0,0,0), c(0.25, 0.25, 0.25, 0.25))
		})
	return(scores)
}

get.simplicity.score.subtypes <- function(state_comp){
	scores <- apply(state_comp, 1, function(x){
		eucDist(x, c(0.25,0.25,0.25)) / eucDist(c(1,0,0), c(0.25, 0.25, 0.25))
		})
	return(scores)
}