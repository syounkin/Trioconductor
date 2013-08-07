THISPKG <- "Trioconductor"
.onAttach <- function(libname, pkgname) {
	version <- packageDescription("Trioconductor", fields="Version")
	packageStartupMessage(paste("
Welcome to Trioconductor version ", version, "\n",
"~~~~~~~~~~~~~~~~~", "\n",
">()_  >()_  >()_ ", "\n",
" (__)  (__)  (__)", "\n",
"~~~~~~~~~~~~~~~~~", sep = "" ) )
}

.onUnload <- function(libpath){
	library.dynam.unload(THISPKG, libpath)
}
