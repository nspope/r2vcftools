#' @include vcflink.R utils.R
NULL

#' Filter VCF files 
#'
#' Given a "vcfLink" object, removes loci/samples corresponding to a set of filtering options. 
#'
#' @param object an object of class "vcfLink"
#' @param arguments optional arguments for vcftools as returned by \code{filterOptions}, see vcftools manual.
#' @param indels what to do with indels; "keep" leaves them in, "remove" removes them, and "only" discards everything but indels.
#' @param verbose print stdout from vcftools? 
#' @rdname Filter 
#' @export
setGeneric("Filter", function(object, ...) useMethod("Filter"))

#' @rdname Filter 
#' @export
setMethod("Filter", "vcfLink", function(object, arguments = filterOptions(), indels = c("keep", "remove", "only"), verbose=TRUE ){
  indels <- match.arg(indels)
	vcf_scratch <- tempfile_vcfLink()
	argvec <- arguments
	argvec <- paste0("--", names(argvec), " ", argvec)
  if(indels=="remove")
  {
    argvec <- c("--remove-indels", argvec)
  } 
  else if(indels=="only")
  {
    argvec <- c("--keep-only-indels")
  }
	argvec <- paste("--vcf", object@vcf_scratch, paste(argvec, collapse=" "), "--recode --out", vcf_scratch)
	res <- system2("vcftools", argvec, stdout = TRUE, stderr = TRUE) 
	if(verbose) cat(res, sep="\n")
	rename <- file.rename(paste0(vcf_scratch,".recode.vcf"), vcf_scratch)
	if(!rename) stop("Problem creating new temporary file")
	object <- UpdateLink(object, vcf_scratch)
	object@site_id <- .getRowIds(vcf_scratch) 
	return(object)
}) 

#' @rdname Filter 
#' @export
filterOptions <- function(min.meanDP = NULL,
			  max.meanDP = NULL,
			  minDP = NULL,
			  maxDP = NULL,
			  min.alleles = NULL,
			  max.alleles = NULL,
			  maf = NULL,
			  max.maf = NULL,
			  minQ = NULL,
			  max.missing = NULL,
			  hwe = NULL,
        thin = NULL){
	c("min-meanDP"=min.meanDP, "max-meanDP"=max.meanDP,
	  "minDP"=minDP, "maxDP"=maxDP,
	  "min-alleles"=min.alleles, "max-alleles"=max.alleles,
	  "maf"=maf, "max-maf"=max.maf,
	  "minQ"=minQ, "max-missing"=max.missing, "hwe"=hwe,
    "thin"=thin)
}

#' Subset VCF files 
#'
#' Given a "vcfLink" object, subsets loci/samples to a given vector of sample or site IDs. 
#'
#' @param object an object of class "vcfLink"
#' @param samples a character vector containing either positions or names of samples to keep
#' @param site a character vector containing either positions or names of sites to keep
#' @param verbose print stdout from vcftools? 
#' @rdname Subset 
#' @export
setGeneric("Subset", function(object, ...) standardGeneric("Subset") )

#' @rdname Subset 
#' @export
setMethod("Subset", "vcfLink", function(object, samples=NULL, sites=NULL, verbose=TRUE){
		if( !is.null(samples) ){
			# standardize names, order
			if( is.numeric(samples) ) {
				if( max(samples) > length(object@sample_id ) )
					stop("Max index exceeds number of samples")
				samples <- object@sample_id[samples]
			}
			if( !all( samples %in% object@sample_id ) )
				stop("Not all 'samples' found in sample ids")
			norder <- order(match(samples, object@sample_id))
			samples <- samples[norder]

			# temporary files and subsetting 
			tmpkeep <- tempfile_vcfLink()
			vcf_scratch <- tempfile_vcfLink()
			cat( samples, file=tmpkeep, sep="\n" )
			argkeep <- paste("--keep", tmpkeep)	
			argkeep <- paste("--vcf", object@vcf_scratch, argkeep, "--recode --out", vcf_scratch)
			out <- system2("vcftools", argkeep, stdout=TRUE, stderr=TRUE)
			if(verbose) cat(out, sep="\n")
			rename <- file.rename(paste0(vcf_scratch,".recode.vcf"), vcf_scratch)
			object <- UpdateLink(object, vcf_scratch)
			
			object@meta <- object@meta[samples, ]
			object@sample_id <- samples
		}
		if( !is.null(sites) ){
			if( is.numeric(sites) ){
				if( max(sites) > length(object@site_id ) )
					stop("Max index exceeds number of sites")
				sites <- object@site_id[sites]
			}
			if( !all( sites %in% object@site_id ) )
				stop("Not all 'sites' found in site ids")
			norder <- order(match(sites, object@site_id))
			sites <- sites[norder]

			# temporary files and subsetting 
			tmpkeep <- tempfile_vcfLink()
			vcf_scratch <- tempfile_vcfLink()
			cat( sites, file=tmpkeep, sep="\n" )
			argkeep <- paste("--snps", tmpkeep)	
			argkeep <- paste("--vcf", object@vcf_scratch, argkeep, "--recode --out", vcf_scratch)
			out <- system2("vcftools", argkeep, stdout=TRUE, stderr=TRUE)
			if(verbose) cat(out, sep="\n")
			rename <- file.rename(paste0(vcf_scratch,".recode.vcf"), vcf_scratch)
			object <- UpdateLink(object, vcf_scratch)
	
			object@site_id <- sites 
		}
		return(object)
}
)
