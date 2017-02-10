.onUnload <- function (libpath) {
  library.dynam.unload("FarmCPUpp", libpath)
}
