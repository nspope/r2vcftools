#' @include vcflink.R utils.R
NULL

#' Fst calculation from VCF file 
#'
#' Given a "vcfLink" object and a list of population IDs, returns a data frame containing pairwise Fst 
#'
#' @param object an object of class "vcfLink"
#' @param by a vector containing population assignments for each sample
#' @param mean_Fst return an estimate of Fst averaged over loci, or separate Fst values for each locus?
#' @rdname Fstats
#' @export
setGeneric(name="Fstats", function(object, ...) standardGeneric("Fstats") )

#' @rdname Fstats
#' @export
setMethod("Fstats", signature="vcfLink", function(object, by, mean_Fst=TRUE){
	# checks
	if( length(by) != length(object@sample_id) )
		stop("'by' must be of same length as sample id")
	if( length(unique(by)) == 1 )
		stop("'by' must have at least two unique values")
	P <- length(unique(by))
	indList <- split( object@sample_id, by ) 
	indList <- lapply( indList, paste, collapse = "\n" )
	fnames <- sapply(1:P, function(i) tempfile_vcfLink() )
	for(i in 1:P)
		cat(indList[[i]], file=fnames[i])

	#TODO: fix the stupid way of counting loci
	loci <- nrow(Query(object, "site-depth", verbose=FALSE))
	fstMat <- matrix(0, loci, P*(P-1)/2)
	colnames(fstMat) <- rep("blank", P*(P-1)/2 )
	z = 1
	for(i in 1:(P-1)){
		for(j in (i+1):P){
			mo <- system2( "vcftools", paste0("--vcf ", object@vcf_scratch, 
				" --weir-fst-pop ", fnames[i],
				" --weir-fst-pop ", fnames[j] 
				), stderr = TRUE, stdout = TRUE)
			me <- read.table("out.weir.fst", header=TRUE)
			fstMat[,z] <- me[,3]
      cnms <- ifelse(is.numeric(by), paste0("pop", unique(by)[i], "-", "pop", unique(by)[j]),
                     paste0(unique(by)[i], "-", unique(by)[j]))
			colnames(fstMat)[z] <- cnms
			z = z + 1
		}
	}
	rownames(fstMat) <- me[,1]
	if( mean_Fst ){
		fstMat <- apply(fstMat, 2, mean, na.rm=TRUE)
		nms <- matrix(Reduce(rbind, strsplit(names(fstMat), split="-", fixed=TRUE)), ncol=2)
		names(fstMat) <- NULL
		colnames(nms) <- c("Population1", "Population2")
		fstMat <- data.frame(nms, mean_fst=fstMat)
	}
	fstMat
})

#' Relatedness calculation from VCF file 
#'
#' Given a "vcfLink" object, returns a data frame containing pairwise relatedness coefficients between samples. 
#'
#' @param object an object of class "vcfLink"
#' @param type the type of relatedness coefficient to return; see vcftools manual
#' @param verbose print stdout from vcftools? 
#' @rdname Relatedness 
#' @export
setGeneric("Relatedness", function(object, ...) standardGeneric("Relatedness") )

#' @rdname Relatedness 
#' @export
setMethod("Relatedness", "vcfLink", function(object, type = c("yang", "manichaikul"), verbose=TRUE)
	  {
      type <- match.arg(type)
		  tmpout <- tempfile_vcfLink(); tmperr <- tempfile_vcfLink();
		  argval <- paste0("--relatedness", ifelse(type=="yang", "", "2"))
		  argval <- paste("--vcf", object@vcf_scratch, argval, "--stdout")
		  out <- system2("vcftools", argval, stdout = tmpout, stderr = tmperr)
		  if(verbose)
			  writeLines(readLines(tmperr))
		  read.table(tmpout, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	  }
)
