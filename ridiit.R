riditize=function(x, latent=FALSE, generate=FALSE)
{
	# THIS FUNCTION DOES RIDIT SCORING (AGRESTI, 2010) FOR A SET OF ORDERED FACTORS.
	# THE SCORES ASSIGNED ARE BASED ON CUMULATIVE PROBABILITIES. 
	# THAT IS, THE PROBABILITY EXHAUSTED BY THE FIRST K CATEGORIES PLUS HALF OF THE PROBABILITY OF K+1. 
	# IT ALSO CAN REPLACE THE VALUES OF X WITH THESE SCORES. 
	# ALTERNATIVELY IT CAN RETURN THE RESPECTIVE QUANTILES OF THE STANDARD NORMAL CDF THAT THE RIDIT VALUES CORRESPOND TO.
	# CHECKS TO SEE IF FUNCTION TAKES ORDERED FACTOR AS ARGUMENT
	if (is.ordered(x)==FALSE)
	{
		stop("must take ordered factor as argument")	
	}
	xlevels=levels(x); length_levels=length(xlevels)
	
	# CHECKS TO SEE IF FUNCTION IS AT LEAST THREE CATEGORIES.
	if (length_levels<3)
	{
		stop("ridit scoring requires 3 or more categories")	
	}
	
	# COMPUTES RIDIT SCORES
	lex=length(x)
	counts=NULL
	for (i in 1:length_levels)
	{
		sump=0
		for (j in 1:lex)
		{
			if (xlevels[i]==x[j])
			{
			sump=sump+1		
			}
		}
		counts[i]=sump
	}
	pk=counts/sum(counts)
	ridits=NULL
	for (i in 1:length_levels)
	{
		if (i==1)
		{
			ridits[i]=(pk[i]/2)
		}
		else
		{
			ridits[i]=sum(pk[1:(i-1)])+(pk[i]/2)
		}
	}
	# PRODUCES VALUES FROM NORMAL CDF FROM RIDITS.
	latentnormalridits=qnorm(ridits)
	# REPLACES ORIGINAL VECTOR VALUES WITH RIDIT SCORES
	if (generate==TRUE)
	{
		riditx=NULL
		for (i in 1:lex)
		{
			for (j in 1:length_levels)
			{
				if (xlevels[j]==x[i] && latent==FALSE)
				{
					riditx[i]=ridits[j]
				}
				else if (xlevels[j]==x[i] && latent==TRUE)
				{
					riditx[i]=latentnormalridits[j]
				}
			}
		}
		return(riditx)
	}
	# RETURNS RIDIT VALUES.
	else
	{
		if (latent==TRUE)
		{	
			return(latentnormalridits)	
		}
		else
		{
			return(ridits)	
		}
	}
}
