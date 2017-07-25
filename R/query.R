#' @include vcflink.R utils.R
NULL

#' Extract SNP ids and positions
#'
#' Given a "vcfLink" object, returns a data frame with the chromosome (contig), position, and ID of each SNP.
#'
#' @param object an object of class "vcfLink"
#' @rdname Chrom
#' @export
setGeneric(name="Chrom", function(object, ...) standardGeneric("Chrom") )

#' @rdname Chrom
#' @export
setMethod("Chrom", "vcfLink", function(object)
          {
            arg <- paste("cat", object@vcf_scratch, "| perl -pe 's/^#.+\n//g' | awk '{print $1}'")
            chrom <- system(arg, intern=TRUE)
            arg <- paste("cat", object@vcf_scratch, "| perl -pe 's/^#.+\n//g' | awk '{print $2}'")
            pos <- system(arg, intern=TRUE)
            cbind("CHROM"=chrom, "POS"=pos, "ID"=object@site_id)
          }
)

#' Query attributes of a VCF file
#'
#' Use vcftools to extract information from a VCF file.
#'
#' @param object an object of class "vcfLink"
#' @param type the type of information to return; see the vcftools manual for details
#' @param verbose boolean; if true, return stdout from vcftools
#' @rdname Query
#' @export
setGeneric(name="Query", function(object, ...) standardGeneric("Query") )

#' @rdname Query
#' @export
setMethod("Query", "vcfLink", function(object, type = c("counts2", "freq2", "site-depth", "site-mean-depth", "geno-depth", "depth", "site-quality", "het", "hardy", "site-pi"), verbose = TRUE)
	  {
    type <- match.arg(type)
    tmpout <- tempfile_vcfLink(); tmperr <- tempfile_vcfLink()
		argvec <- paste("--vcf", object@vcf_scratch, paste0("--", type), paste0("--stdout") )
		out <- system2("vcftools", argvec, stdout=tmpout, stderr=tmperr)
		if( verbose )
			writeLines(readLines(tmperr))
		colnms <- strsplit( readLines(tmpout, 1) , split="\t" )[[1]]
		ncols <- max( count.fields(tmpout, sep="\t") )
		if( length(colnms) < ncols )
			colnms <- c(colnms, paste0(colnms[length(colnms)], 1:(ncols-length(colnms)) ) )
		read.table(tmpout, skip=1, sep="\t", fill=TRUE, col.names=colnms, stringsAsFactors=FALSE)
	  }
)

#' Calculate linkage disequillibrium between loci 
#'
#' Given a "vcfLink" object, returns a data frame with pairwise p-values between loci, under the null hypothesis that the loci are not in linkage disequillibrium. 
#'
#' @param object an object of class "vcfLink"
#' @param type the type of comparison to make. See vcftools manual. 
#' @param arguments optional arguments for vcftools as returned by \code{linkageOptions}, see vcftools manual.
#' @param verbose print stdout from vcftools? 
#' @rdname Linkage 
#' @export
setGeneric(name="Linkage", function(object, ...) standardGeneric("Linkage") )

#' @rdname Linkage 
#' @export
setMethod("Linkage", "vcfLink", function(object, type=c("geno-r2", "geno-chisq", "interchrom-geno-r2"), arguments=linkageOptions(), verbose=TRUE)
          {
            type <- match.arg(type)
            tmpout <- tempfile_vcfLink(); tmperr <- tempfile_vcfLink()
            argvec <- arguments
            if( !is.null(argvec) )
              argvec <- paste0("--", names(argvec), " ", argvec )
            argvec <- paste("--vcf", object@vcf_scratch,  paste0("--", type), paste(argvec, collapse=" "), paste0("--stdout") )
            out <- system2("vcftools", argvec, stdout=tmpout, stderr=tmperr)
            if(verbose)
              writeLines(readLines(tmperr))
            ld <- read.table(tmpout, header=TRUE)
            pos <- Chrom(object)
            pos_unique <- paste(pos[,1], pos[,2], sep="-")
            if("CHR" %in% colnames(ld))
              ld$CHR1 = ld$CHR2 = ld$CHR
            ld$ID1 <- pos[match(paste(ld$CHR1, ld$POS1, sep="-"), pos_unique), "ID"]
            ld$ID2 <- pos[match(paste(ld$CHR2, ld$POS2, sep="-"), pos_unique), "ID"]
            ld <- ld[,c("CHR1","POS1","ID1","CHR2","POS2","ID2","N_INDV","R.2")]
            attr(ld, "snpid") <- unique(c(ld$ID1, ld$ID2))
            return(ld)
          }
)

#' @rdname Linkage
#' @export
linkageOptions <- function(min.r2=NULL,
                           ld.window=NULL,
                           ld.window.bp=NULL,
                           ld.window.min=NULL,
                           ld.window.bp.min=NULL)
{
  c("min-r2"=min.r2,
              "ld-window"=ld.window,
              "ld-window-bp"=ld.window.bp,
              "ld-window-min"=ld.window.min,
              "ld-window-bp-min"=ld.window.bp.min)
}
