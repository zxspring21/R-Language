quit <- function() {
  .Tcl("set ok_to_quit 1")
  tkdestroy(tt)
  tkdestroy(QI_interface)
}

SAQCgui=paste(getwd(), "/PROGRAM/SAQC_interactive.r",sep="")
source(SAQCgui)

tkwait.window(QI_interface)

	# put below the "RUN" button
	# tt <<- tktoplevel()
	# OK.but <- tkbutton(tt,text="OK",command=quit)
	# tkgrid(OK.but)
	# tkfocus(tt) 