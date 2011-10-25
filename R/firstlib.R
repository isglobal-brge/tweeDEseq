
############ First.lib ###############

.onLoad <- function(lib, pkg){
   library.dynam("tweeDEseq", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("tweeDEseq", libpath)


############ End of .First.lib ###############


