THISPKG <- "trioClasses"
.onAttach <- function(libname, pkgname) {
	version <- packageDescription("trioClasses", field="Version")
	packageStartupMessage(paste("
==================================
Welcome to trioClasses version ", version, "\n",
"by Scharpf and Younkin
==================================", sep = "" ) )
}

.onUnload <- function(libpath){
	library.dynam.unload(THISPKG, libpath)
}
