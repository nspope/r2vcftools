.checkDependencies <- function()
{
  # I'm sure there's a better way to do this.
  sout <- tempfile()
  serr <- tempfile()

  cat("Checking for awk ... ")
  awk <- system2("awk", paste0("'{print $0}' ", sout), stderr=serr, stdout=sout)
  if(awk!=0)
    stop("awk not found!", call. = FALSE)
  cat("OK\n")

  cat("Checking for perl ... ")
  perl <- system2("perl", "-h", stderr=serr, stdout=sout)
  if(perl!=0)
    stop("perl not found!", call. = FALSE)
  cat("OK\n")

  cat("Checking for vcftools ... ")
  vcftools <- system2("vcftools", "-h", stderr=serr, stdout=sout)
  if(vcftools!=0)
    stop("vcftools not found!", call. = FALSE)
  cat("OK\n")

}

.checkDependencies()
