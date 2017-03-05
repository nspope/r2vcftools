#' @include vcflink.R
NULL

#' @title Clean unused temporary VCF files 
#'
#' @description
#' Sweeps the temporary directory and clears temporary VCF files that are not currently linked to an object.
#'
#' @param object an object of class "vcfLink"
#' @rdname Clean
#' @export
setGeneric("Clean", function(object, ...) standardGeneric("Clean") )

#' @rdname Clean
#' @export
setMethod("Clean", "vcfLink", function(object){
	envi <- environment() 
	oln <- c()
	while( !(identical(envi, .GlobalEnv)) ){
		ol <- ls(envir=envi)
		olc <- sapply( ol, function(i) class(get(i, envir=envi) ))
		oln <- c(oln, sapply( ol[olc == "vcfLink"], function(i) get(i, envir=envi)@vcf_scratch ) )
		envi <- parent.env(envi)
	}
	ol <- ls(envir=.GlobalEnv)
	olc <- sapply( ol, function(i) class(get(i, envir=.GlobalEnv) ))
	oln <- c(oln, sapply( ol[olc == "vcfLink"], function(i) get(i, envir=.GlobalEnv)@vcf_scratch) )
	tf <- paste0(tempdir(), "/", list.files(tempdir()))
	tfd <- tf %in% c(oln)
	sapply(tf[!tfd], function(del) if(grepl("vcfLink_file",del)) file.remove(del))
}
)

#' @title Update link to temporary VCF
#'
#' @description
#' Links a vcfLink to a new temporary path
#'
#' @param object an object of class "vcfLink"
#' @param newlink character; new path to update object with
#' @param clean boolean; clean temporary directory? If true, the old file that the object is linked to will be deleted 
#' @rdname UpdateLink 
#' @export
setGeneric("UpdateLink", function(object, newlink, ...) standardGeneric("UpdateLink") )

#' @rdname UpdateLink 
#' @export
setMethod("UpdateLink", "vcfLink", function(object, newlink, clean=TRUE){
	  object@vcf_scratch <- newlink
	  if(clean)
		  Clean(object)
	  return( object )
}
)

#' @title Copy vcfLink object 
#'
#' @description
#' Copies a vcfLink object, and creates a copy of the temporary (working) VCF file for it to link to
#'
#' @param object an object of class "vcfLink"
#' @rdname Copy 
#' @export
setGeneric("Copy", function(object, ...) standardGeneric("Copy") )

#' @rdname Copy 
#' @export
setMethod("Copy", "vcfLink", function(object){
		  newlink <- tempfile_vcfLink()
		  file.copy(object@vcf_scratch, newlink)
		  object <- UpdateLink(object, newlink)
		  return(object) 
})

#' @title Save vcfLink object 
#'
#' @description
#' Saves the temporary VCF linked to by a vcfLink object to a given filename. Also, saves any sample metadata which is associated with the VCF file.
#'
#' @param object an object of class "vcfLink"
#' @param filename character; the path to which to save the VCF file to
#' @rdname Save 
#' @export
setGeneric("Save", function(object, ...) standardGeneric("Save") )

#' @rdname Save 
#' @export
setMethod("Save", "vcfLink", function(object, filename){
	me <- system2("cp", paste(object@vcf_scratch, filename), stdout=TRUE, stderr=TRUE)
	meta_path <- strsplit(filename, ".", fixed=TRUE)[[1]]
	meta_path <- paste(meta_path[1:length(meta_path)], collapse=".")
	meta_path <- paste0(meta_path, ".meta")
	write.table(object@meta, meta_path, sep="\t", quote=FALSE) 
	object@vcf_path <- filename
	object@meta_path <- meta_path
	return(object)
})

#' @title Attach metadata to vcfLink object 
#'
#' @description
#' Attaches a dataframe containing sample metadata to a vcfLink object. 
#'
#' @param object an object of class "vcfLink"
#' @param meta a "data.frame" containing sample metadata. Columns correspond to variables, rows correspond to samples.
#' @rdname LoadMeta 
#' @export
setGeneric("LoadMeta", function(object, ...) standardGeneric("LoadMeta"))

#' @rdname LoadMeta 
#' @export
setMethod("LoadMeta", "vcfLink", function(object, meta){
  sample_id <- rownames(meta)
  snames <- object@sample_id
  if( is.null(sample_id) | all(sample_id == 1:nrow(meta)) | any(duplicated(sample_id)) )
    stop("Meta must have row names that are unique sample IDs")
  if( !all(snames %in% sample_id) )
    stop("Some VCF entries not found in 'sample_id'")
  if( !all(sample_id %in% snames) ){
    warning("Some of 'sample_id' not found in VCF, subsetting to only those found in VCF")
    meta <- meta[sample_id %in% snames,]
    sample_id <- sample_id[sample_id %in% snames]
  }
  # reorder metadata 
  norder <- match(sample_id, snames)
  sample_id <- sample_id[norder]
  if( !all(sample_id == snames) )
    stop("Problem matching samples -- check for duplicates")
  meta <- meta[norder,]
  object@meta <- data.frame(object@meta, meta)
  return(object)
})
