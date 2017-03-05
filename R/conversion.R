#' @include vcflink.R
NULL

#' Convert VCF to Plink 
#'
#' Given a "vcfLink" object, converts the VCF file to Plink format and returns the locations of the .ped and .map files.
#'
#' @param object an object of class "vcfLink"
#' @param output.prefix the prefix to use for output files
#' @param verbose print stdout from vcftools?
#' @rdname Plink 
#' @export
setGeneric("Plink", function(object, ...) standardGeneric("Plink") )

#' @rdname Plink 
#' @export
setMethod("Plink", "vcfLink", function(object, output.prefix=NULL, verbose=TRUE) {
    # TODO: Admixture doesn't seem to like these, might have to use different software
		if(is.null(output.prefix))
			output.prefix <- sub(object@vcf_path, pattern="\\.vcf", replacement="")
		argvar <- paste("--vcf", object@vcf_scratch, "--plink", "--out", output.prefix)
		me <- system2("vcftools", argvar, stdout=TRUE, stderr=TRUE)
		if(verbose){
			writeLines(me)
			cat("Returning location of plink files.", sep="\n")
		}
		return(paste0(output.prefix, c(".ped", ".map")))
	}
)

#' Convert VCF to genotype matrix 
#'
#' Given a "vcfLink" object, converts the VCF file to genotype matrix format. Doing so involves dropping indels and multiallelic SNPs, so this fuction returns a new vcfLink object that has been appropriately subset.
#'
#' @param object an object of class "vcfLink"
#' @param output.file the filename to save the genotype matrix to
#' @param verbose print stdout from vcftools?
#' @rdname Geno
#' @export
setGeneric("Geno", function(object, ...) standardGeneric("Geno") )

#' @rdname Geno
#' @export
setMethod("Geno", "vcfLink", function(object, output.file, verbose=TRUE) {
  filename <- paste0(object@vcf_scratch, ".vcf")
  file.copy(object@vcf_scratch, filename, overwrite=TRUE)
  geno <-  vcf2geno(filename, output.file=output.file, force=TRUE)
  if(!grepl("\\.geno", output.file))
  {
    output.file = paste0(output.file, ".geno")
    warning(paste("Adding suffix '.geno' to output.file: now\n",output.file))
  }
  removed <- sub("\\.geno", "\\.removed", output.file)
  vcfsnp <- sub("\\.geno", "\\.vcfsnp", output.file)
  #if(loadFile){
  #    geno <- read.geno(geno)
  #    removed <- read.table(removed, header=FALSE)
  #    vcfsnp <- read.table(vcfsnp, header=FALSE)
  #}
  file.remove(filename)

  # subsetting based on removed SNPs
  if(file.info(removed)$size > 0)
  {
    removed_snps <- read.table(removed, sep=" ")$V3
    warnings("Subsetting vcf to biallelic SNPs to match geno file")
    snps_geno <- object@site_id[!(object@site_id %in% removed_snps)] #
    object <- Subset(object, sites=snps_geno, verbose=verbose)
  }

  if(verbose)
    cat("Output files:", output.file, removed, vcfsnp, sep="\n")
  return(object)
})

#' Convert VCF to genotype matrix and return
#'
#' Given a "vcfLink" object, converts the VCF file to genotype matrix format and returns this matrix. See the "--012" flag in the vcftools documentation.
#'
#' @param object an object of class "vcfLink"
#' @param output.file the filename to save the genotype matrix to, if you want to write it outside of R
#' @param verbose print stdout from vcftools?
#' @rdname GenotypeMatrix 
#' @export
setGeneric("GenotypeMatrix", function(object, ...) standardGeneric("GenotypeMatrix") )

#' @rdname GenotypeMatrix 
#' @export
setMethod("GenotypeMatrix", "vcfLink", function(object, filename=NULL, verbose=TRUE) {
		if(is.null(filename))
			filename <- tempfile_vcfLink() 
    object <- Filter(object, filterOptions(max.alleles=2, min.alleles=2), verbose=FALSE)
		argvar <- paste("--vcf", object@vcf_scratch, "--012", "--out", filename)
		me <- system2("vcftools", argvar, stdout=TRUE, stderr=TRUE)
		if( verbose )
			writeLines(me)
		out <- as.matrix(read.table(paste0(filename,".012"), sep="\t", row.names=object@sample_id))
    colnames(out) = c("sample", object@site_id)
		out[,-c(1)]
})

#' Convert VCF to LFMM format
#'
#' Given a "vcfLink" object, converts the VCF file to the format used by LFMM in package LEA. Involves intermediate creation of a genotype matrix, so only biallelic SNPs are retained.
#'
#' @param object an object of class "vcfLink"
#' @param output.file the filename to save the LFMM file to 
#' @param verbose print stdout from vcftools?
#' @rdname Lfmm 
#' @export
setGeneric("Lfmm", function(object, ...) standardGeneric("Lfmm") )

#' @rdname Lfmm 
#' @export
setMethod("Lfmm", "vcfLink", function(object, output.file, verbose=TRUE) {
  tmpfile <- paste0(tempfile_vcfLink(), ".geno")
  if(verbose)
    cat("\nConverting to '.geno'", sep="\n")
  object <- Geno(object, output.file=tmpfile, verbose=verbose)
  if(verbose)
    cat("\nConverting to '.lfmm'", sep="\n")
  output.file <- geno2lfmm(tmpfile, output.file=output.file)

  if(verbose)
    cat("Output files:", output.file, sep="\n")
  return( object )
})

#' Create environmenal variable file for LFMM 
#'
#' Given a "vcfLink" object, writes given variables in the metadata to a .env file used by LFMM in package LEA. 
#'
#' @param object an object of class "vcfLink"
#' @param varnames a character vector of the variables in the metadata to write to the .env file
#' @param output.file the filename to save the ENV file to
#' @rdname LfmmEnv
#' @export
setGeneric("LfmmEnv", function(object, ...) standardGeneric("LfmmEnv") )

#' @rdname LfmmEnv
#' @export
setMethod("LfmmEnv", "vcfLink", function(object, varnames, output.file = NULL) {
  if( !all( varnames %in% colnames(object@meta) ) )
    stop("Some variables not found in metadata")
  tempenv <- output.file 
  if(is.null(output.file))
    tempenv <- paste0(tempfile_vcfLink(), ".env")
  write.table(object@meta[,varnames], tempenv, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  return(tempenv)
})
