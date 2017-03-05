.setRowIds <- function(filename){
  arg <- paste("cat", filename, "| perl -pe 's/^#.+\n//g' | awk '{OFS=\"\\t\"; $3=\"snp\"NR; print $0}' > vlinesplit")
  system(arg, intern=FALSE)
  system( paste( "grep -E '^#'", filename, "> header") )
  system( paste("cat header vlinesplit >", filename) )
  system( "rm vlinesplit header")
}

.getRowIds <- function(filename){
 	arg <- paste("cat", filename, "| perl -pe 's/^#.+\n//g' | awk '{print $3}'")
	out <- system(arg, intern=TRUE)
	out
}

tempfile_vcfLink <- function(){
  tempfile(pattern = "vcfLink_file")
}

#' An S4 class containing a link to a VCF file which vcftools interacts with.
#'
#' @slot vcf_path The path to the original VCF file
#' @slot vcf_scratch The path to the temporary (working) VCF file
#' @slot meta_path The path to a csv file containing metadata describing the samples
#' @slot sample_id A character vector containing the names of samples
#' @slot meta A data frame containing metadata describing the samples
#' @slot site_id A character vector containing the names of SNPs
#' @export
setClass(
	"vcfLink",
	## members
	representation(
		vcf_path = "character",
		vcf_scratch = "character",
		meta_path = "character",
		sample_id = "character",
		meta = "data.frame",
		site_id = "character"
	)
)

#' Create a vcfLink object
#'
#' Creates an object of class vcfLink, which acts as an interface to a temporary VCF file. Also loads metadata (if present) associated with the samples.
#'
#' @param vcf_path Character, the path to a VCF file
#' @param meta Character, the path to a csv containing metadata associated with each sample
#' @param overwriteID Boolean, overwrite SNP ids with numbers? Only set to TRUE if the VCF file does not already contain SNP ids 
#' @export
vcfLink <- function(vcf_path, meta=NULL, overwriteID=FALSE) {
	# if meta already exists and is not supplied, load it and issue warning
	meta_path <- strsplit(vcf_path, ".", fixed=TRUE)[[1]]
	meta_path <- paste(meta_path[1:length(meta_path)], collapse=".")
	meta_path <- paste0( meta_path, ".meta")
	if( !file.exists(meta_path) & is.null(meta) ){
		warning("Meta not supplied and does not already exist") 
	} else if( file.exists(meta_path) & is.null(meta) ){
		# load current meta
		meta <- read.table(meta_path, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
	} else if( exists(meta_path) & !is.null(meta) ){
			warning("Meta for VCF already exists but is being ignored")
	} else {}

	# create temporary file
	vcf_scratch <- tempfile_vcfLink()
	file.copy(from = vcf_path, to = vcf_scratch) 

	# parse sample names and match to meta
	snames <- system2("grep", paste("-Eo '^#CHROM.+'", vcf_scratch), stdout=TRUE, stderr=TRUE)
	snames <- strsplit(snames, split="\t")[[1]][-c(1:9)]
	if( is.null(meta) ){
		meta <- data.frame(sample_num = 1:length(snames), sample_name = snames, stringsAsFactors=FALSE)
		rownames(meta) <- meta$sample_name
	}
	sample_id <- rownames(meta)
	if( is.null(sample_id) | all(sample_id == 1:nrow(meta)) | any(duplicated(sample_id)) )
		stop("Meta must have row names that are unique sample IDs")
	if( !all(snames %in% sample_id) )
		stop("Some VCF entries not found in 'sample_id'")
	if( !all(sample_id %in% snames) )
		stop("Some of 'sample_id' not found in VCF")

	# reorder metadata 
	norder <- match(sample_id, snames)
	sample_id <- sample_id[norder]
	if( !all(sample_id == snames) )
		stop("Problem matching samples -- check for duplicates")
	meta <- meta[norder,]

	# if overwriting site ids
	if(overwriteID){
		warning("Overwriting SNP ids")
		.setRowIds(vcf_scratch)
	}
	site_id = .getRowIds(vcf_scratch)

	new("vcfLink", meta=meta, vcf_path=vcf_path, vcf_scratch=vcf_scratch, 
	    meta_path=meta_path, sample_id=sample_id, site_id=site_id)
}

