THISPKG <- "trioClasses"
.onAttach <- function(libname, pkgname) {
	version <- packageDescription("trioClasses", field="Version")
	packageStartupMessage(paste("
Welcome to trioClasses version ", version, "\n",
"~~~~~~~~~~~~~~~~~", "\n",
">()_  >()_  >()_ ", "\n",
" (__)  (__)  (__)", "\n",
"~~~~~~~~~~~~~~~~~", sep = "" ) )
}

.onUnload <- function(libpath){
	library.dynam.unload(THISPKG, libpath)
}
