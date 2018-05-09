m.chi <- function(tm, tm.h0){
	rs2 <- rowSums(tm.h0)
	rs1 <- rowSums(tm)
	rs2nz <- rs2 > 0
	rs1nz <- rs1 > 0
	dof1 <- sum(rs1nz)
	dof2 <- sum(rs2nz)
	rs2 <- rs2 + (rs2 == 0)
	dof <- (dof1 - 1) * (dof2 - 1)
	p <- diag(1/rs2) %*% tm.h0
	E <- diag(rs1) %*% p
	num <- (tm - E)^2
	E <- E + (E == 0)
	chi2 <- sum(num/E)
	pvalue <- 1-pchisq(chi2,dof)
	output <-c("Chi2"=chi2, "p-value"=pvalue, "d.o.f"=dof)
	return(output)
	}

