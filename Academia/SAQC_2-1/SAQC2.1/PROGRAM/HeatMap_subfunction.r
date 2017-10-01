
PolyGonPlot_interactive=function(a,ub,ind_name,QIabline)
{
	#a=z;ub=UBQI ;ind_name=indname;QIabline=QIline
	ub=t(ub)
	ind_name=as.character(ind_name)
	nind=nrow(a)
	nchip=ncol(a)

	ub_matrix=rep(QIabline,nind);
	ub_matrix=matrix(ub_matrix,nind,nchip,byrow=T)
	z=a-ub_matrix
	z[z>=0]=1;z[z<0]=0
	z=as.matrix(z)
	z=matrix(as.numeric(z),nind,nchip)
	y=c(1:nchip);x=c(1:nind)
	nY=3
	y.lim = if(index.2chip){
		3
	} else 1
	#if(option=="2chip"){nY=3}
	#if(option=="1chip"){nY=1}

	if(sum(z)==nind*nchip){
		#image(x,y,z,col=c(terrain.colors(18, alpha = 1)[18]),axes=F,ylab="",xlab="",ylim=c(0.5,nY+0.5),cex.main=3,font.main=4)
		image(x,y,z,col=c(terrain.colors(18, alpha = 1)[18]),axes=F,ylab="",xlab="",ylim=c(0.5,y.lim+0.5),cex.main=3,font.main=4)
	}
	if(sum(z)==0){
		#image(x,y,z,col=c(terrain.colors(18, alpha = 1)[1]),axes=F,ylab="",xlab="",ylim=c(0.5,nY+0.5),cex.main=3,font.main=4)
		image(x,y,z,col=c(terrain.colors(18, alpha = 1)[1]),axes=F,ylab="",xlab="",ylim=c(0.5,y.lim+0.5),cex.main=3,font.main=4)
	}
	if(sum(z)<nind*nchip & sum(z)>0){
		#image(x,y,z,col=c(terrain.colors(18, alpha = 1)[1],terrain.colors(18, alpha = 1)[18]),axes=F,ylab="",xlab="",ylim=c(0.5,nY+0.5),cex.main=3,font.main=4)
		image(x,y,z,col=c(terrain.colors(18, alpha = 1)[1],terrain.colors(18, alpha = 1)[18]),axes=F,ylab="",xlab="",ylim=c(0.5,y.lim+0.5),cex.main=3,font.main=4)
	}

	abline(v=x+0.5,col="skyblue",lty=2)
	abline(h=y+0.5,col="skyblue")
	colquantile=c("skyblue","skyblue","skyblue")
	for(i in 1:3){
		#i=1
		cex.line.note=0.8
		tempQI=a[,i];UB=ub[,i];tempQIabline=QIabline[i]
		b=c(tempQI,UB,tempQIabline,0);yQI=range(b,na.rm=T);yQI[2]=yQI[2]*1.1
		newY=(i-0.5)+(b-yQI[1])/(yQI[2]-yQI[1])
		tempQI=newY[x]
		newUB=newY[max(x)+c(1,2,3)]
		newQIabline=newY[max(x)+4]
		#abline(h=newUB,col=colquantile,lwd=1.5,lty=1)
		#abline(h=newQIabline,col="purple",lwd=2)

		#points(c(x,max(x)+1)-0.5,c(tempQI,tempQI[max(x)]),type="s",col="blue",pch="+")
		points(x,tempQI,col="blue",pch="+")
		polygonX=c(0.5,x,nind+0.5)
		polygonY=c(i-0.5,tempQI,i-0.5)
		polygon(polygonX,polygonY,col="yellow",pch="+",density=45,border="gray25")
		UBlab=round(UB*1000)/1000
		axis(4,newQIabline,labels=round(tempQIabline*100)/100,cex.axis=0.75,font.axis=3,col.axis="deeppink2",tick=F,line=-0.75,las=2)
		#leg.txt <- c(paste(UBlab[1],"-(95%)",sep=""),paste(UBlab[2],"-(97.5%)",sep=""),paste(UBlab[3],"-(99%)",sep=""))
		#legend(nind*0.9,yQI[2]*0.8, leg.txt, col=colquantile, cex=0.75)
		
		#text(nind*0.975,newUB[1]*0.99,labels=UBlab[1],cex=cex.line.note,font=2)
		#text(nind*0.175,newUB[2],labels=UBlab[2],cex=cex.line.note,font=2)
		#text(nind*0.975,newUB[3]*1.01,labels=UBlab[3],cex=cex.line.note,font=2)
		UBlab=round(UBlab*100)/100
		d=0.5
		QIat=range(tempQI)
		QIlab=range(a[,i]);QIlab=round(QIlab*100)/100
		axis(4,at=QIat,labels=QIlab,tick=F,line=-0.75,las=2,cex.axis=0.75)
		abline(h=newUB,col=colquantile,lwd=1.5,lty=1)
		abline(h=newQIabline,col="hotpink2",lwd=2.5)
		abline(h=QIat,col="gray25",lwd=1.75,lty=1)
		if(!is.na(newUB[1])) axis(4,as.numeric(newUB[1]*0.99),labels=UBlab[1],col.axis="darkblue",font.axis=4,tick=F,line=-3.25-d,las=2,cex.axis=0.75)
		if(!is.na(newUB[2])) axis(4,as.numeric(newUB[2]*1.00),labels=UBlab[2],col.axis="darkblue",font.axis=4,tick=F,line=-1.75-d,las=2,cex.axis=0.75)
		if(!is.na(newUB[3])) axis(4,as.numeric(newUB[3]*1.01),labels=UBlab[3],col.axis="darkblue",font.axis=4,tick=F,line=-3.25-d,las=2,cex.axis=0.75)

		#text(nind*0.925,newUB[1]*0.99,labels="(95%)",cex=cex.line.note,font=2)
		#text(nind*0.125,newUB[2],labels="(97.5%)",cex=cex.line.note,font=2)
		#text(nind*0.925,newUB[3]*1.01,labels="(99%)",cex=cex.line.note,font=2)
		
	#    QIat=range(tempQI)
	#    QIlab=range(a[,i]);QIlab=round(QIlab*100)/100
	#    axis(4,at=QIat,labels=QIlab,tick=F,line=-0.75,las=2,cex.axis=0.75)
	#    abline(h=newUB,col=colquantile,lwd=1.5,lty=1)
	#    abline(h=newQIabline,col="purple",lwd=2)
	#    abline(h=QIat,col="tan4",lwd=1.5,lty=1)

		#axis(2,at=c(i-0.5,i+0.5),labels=c("",""),tick=T,line=0.25)
	}
	if(index.2chip) abline(h=c(0.5,1.5,2.5),lwd=1.5,col="blue")
	box(col="blue",lwd=1)
	cexsize=1.75;ymove=0.5
	##-------------------------- X axis setting ----------------------------------
	nbroke=10;labsize=0.75
	ntemp=ceiling(nind/nbroke);
	#if(nind>=30)
	#{
	#  xind=c(1:ntemp)*nbroke
	#  if((nind-max(xind))>=ntemp/2){xind=c(xind,nind)}
	#  if((nind-max(xind))<ntemp/2){xind=c(xind[-length(xind)],nind)}
	#}
	#if(nind<30){xind=c(1:nind)}
	xind=c(1:ntemp)*nbroke
	axis(1,at=x,labels=ind_name,tick=F,line=-2.25,cex.axis=0.75,font=4,las=2,col.axis="blue")
	abline(v=xind+0.5,col="blue",lwd=1.5,lty=1)
	abline(v=c(0,nind)+0.5,col="black",lwd=2)
	abline(h=c(0,1,2,3)+0.5,col="blue",lwd=2)
	##-------------------------- Y axis setting ----------------------------------
	ymove=-0.25
	#axis(2,at=c(1:3),labels=c("Chip1","Chip2","Merge"),font.axis=4,tick=F,line=-0.25+ymove,col.axis="blue")
	axis(2,at=1:y.lim,labels=if(index.2chip){ c("Chip1","Chip2","Merge") } else unlist(strsplit(chip_name, " "))[1],
			font.axis=4,tick=F,line=-0.25+ymove,col.axis="blue")
	##--------------------------------------------------
	normal_chip_ratio=1-apply(z,2,mean)
	normal_chip_ratio=as.numeric(normal_chip_ratio)
	normal_chip_ratio=round(normal_chip_ratio*100)/100
	#for(i in 1:3){text(nind,i+0.45,labels=normal_chip_ratio[i],adj=c(1,1),cex=cexsize)}
	##--------------------------------------------------
	return()
}

##==============================================================================
##==============================================================================
##==============================================================================

##------------------------------------------------------------------------------
## subfunction for fillframe, it is for Mytkexamp
# Mytkexamp function from TeachingDemo then modified it by Hsin-Chi
##------------------------------------------------------------------------------
Mytkexamp_heatmap=function (FUN, param.list, vscale = 1.5, hscale = 1.5, wait = FALSE,plotloc = "top",Tkexample_title)  #plottype=1 heatmap; plottype=2 polygon 
{
    #FUN=HeatMap; param.list=HeatMaplist;vscale = 1.5; hscale = 2.75;Tkexample_title="QI HeatMap plot";plotloc = "top"; 
    #FUN=Polygon; param.list=Polygonlist;vscale = 1.5; hscale = 2.75;Tkexample_title="QI Polygon plot";plotloc = "top"
    if (!require("tkrplot")) {stop("The tkrplot package is needed")}
    require(tcltk2)
    ocl <- cl <- substitute(FUN)
    exargs <- as.list(quote(list()))
    replot <- function() {eval(cl)}
    tt <- tktoplevel()
    tkwm.title(tt, Tkexample_title)
    img <- tkrplot(tt, replot, vscale = vscale, hscale = hscale)
    tkpack(img, side = plotloc)
    hsc <- tclVar()
    tclvalue(hsc) <- hscale
    vsc <- tclVar()
    tclvalue(vsc) <- vscale
    ##--------------------------------------------------------------------------
    ##                                Main area
    ##--------------------------------------------------------------------------
    fillframe <- function(frame, lst, pkdir, prfx)
    {
        #frame=tt;lst=param.list;prfx="tkv";pkdir=plotloc
        for (i in seq_along(lst))
        {
            #i=1
            vname <- paste(prfx, ".", i, sep = "")
            el <- lst[[i]]
            eln <- names(lst)[i]
            if (is.list(el[[1]]))
            {
                fr <- tkframe(frame, relief = "ridge", borderwidth = 3)
                tkpack(fr, side = pkdir)
                if (length(eln) && nchar(eln))
                {
                  tkpack(tklabel(fr, text = eln), side = "top",anchor = "nw")
                }
                Recall(fr, el, ifelse(pkdir == "top", "left","top"), vname)
                next
            }
            ##--------------------------- slider -------------------------------
            if (tolower(el[[1]]) == "slider")
            {
                tkpack(fr <- tkframe(frame), side = pkdir)
                tkpack(tklabel(fr, text = eln), side = "left",
                  anchor = "s", pady = 4)
                tmp <- tclVar()
                tclvalue(tmp) <- if ("init" %in% names(el)) {
                  el$init
                }
                else if ("from" %in% names(el)) {
                  el$from
                }
                else {
                  1
                }
                alist <- list(fr, variable = tmp, orient = "horizontal",
                  command = function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)),
                    vscale = as.numeric(tclvalue(vsc))))
                el2 <- el[-1]
                el2$init <- NULL
                alist <- c(alist, el2)
                tkpack(do.call("tkscale", alist), side = pkdir)
                tmpcl <- as.list(cl)
                tmpl <- list(substitute(as.numeric(tclvalue(VNAME)),
                  list(VNAME = as.character(tmp))))
                names(tmpl) <- eln
                cl <<- as.call(c(tmpcl, tmpl))
                exargs <<- c(exargs, tmpl)
                next
            }
            ##--------------------------- checkbox ------------------------------
            if (tolower(el[[1]]) == "checkbox")
            {
                tmp <- tclVar()
                tclvalue(tmp) <- if ("init" %in% names(el)) {el$init}
                else {"F"}                                            ##"T"           "F"
                alist <- list(fr, variable = tmp, text = eln,onvalue = "T", offvalue = "F", command = function(...) tkrreplot(img,hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc))))
                el2 <- el[-1]
                el2$init <- NULL
                tmpvars <- if ("values" %in% names(el)) {el$values}
                else {""}
                el2$values <- NULL
                alist <- c(alist, el2)
                tkpack(do.call("tkcheckbutton", alist), side = pkdir)
                tmpcl <- as.list(cl)
                tmpl <- list(substitute(as.logical(tclvalue(VNAME)),list(VNAME = as.character(tmp))))
                names(tmpl) <- eln
                cl <<- as.call(c(tmpcl, tmpl))
                exargs <<- c(exargs, tmpl)
                next
            }
            ###==================================================================
            ##--------------------------- radiobuttons -------------------------
            if (tolower(el[[1]]) == "radiobuttons")
            {
                tkpack(fr <- tkframe(frame, relief = "groove",
                  borderwidth = 3), side = pkdir)
                tkpack(tklabel(fr, text = eln), side = "top",anchor = "nw")#,font=fontHeading
                tmp <- tclVar()
                tclvalue(tmp) <- if ("init" %in% names(el)) {
                  el$init
                }
                else {
                  el$values[1]
                }
                el2 <- el[-1]
                tmp.vals <- el2$values
                el2$values <- NULL
                el2$init <- NULL
                alist <- list(fr, variable = tmp, command = function() tkrreplot(img,
                  hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc))))
                pkdir2 <- ifelse(pkdir == "top", "left", "top")
                for (v in tmp.vals) {
                  tkpack(do.call("tkradiobutton", c(alist, value = v,
                    text = v)), side = pkdir2)
                }
                tmpcl <- as.list(cl)
                tmpl <- list(substitute(tclvalue(VNAME), list(VNAME = as.character(tmp))))
                names(tmpl) <- eln
                cl <<- as.call(c(tmpcl, tmpl))
                exargs <<- c(exargs, tmpl)
                next
            }
        }
       
    }
    
    ##--------------------------------------------------------------------------
    tkpack(tfr <- tkframe(tt), side = "bottom", fill = "x")

    tkpack(tkbutton(tfr, text = "Refresh",width=20, command = function() {
        tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc)))
    }), side = "left", anchor = "s")
    tkpack(tklabel(tfr, text = "      ",width=50), side = "left")
    
    ##---------modified here for figure print!!-----------------
    tkpack(tkbutton(tfr, text = "Save file",width=20,command = function() {
        tmp <- c(as.list(ocl), eval(as.call(exargs)))
        #cat(deparse(as.call(tmp)), "\n")
        figure_info=(eval(as.call(exargs)))
        #print("figure_info");print(figure_info)
        if(Tkexample_title=="QI HeatMap plot"){HeatMap_output(figure_info,output_popu);cat("HeatMap figure print out.\n")}
        if(Tkexample_title=="QI Polygon plot"){Polygon_output(figure_info,output_popu);cat("Polygon figure print out.\n")}
        #print("------------")
        flush.console()
    }), side = "left", anchor = "s")
    ##===========================================================
    tkpack(tkbutton(tfr, text = "Exit",width=20, command = function() tkdestroy(tt)), side = "right", anchor = "s")
    tkpack(tfr <- tkframe(tt), side = "bottom", fill = "x")
    tkpack(tklabel(tfr, text = "Hscale: "), side = "left")
    tkpack(tkentry(tfr, textvariable = hsc, width = 6), side = "left")
    tkpack(tklabel(tfr, text = "      Vscale: "), side = "left")
    tkpack(tkentry(tfr, textvariable = vsc, width = 6), side = "left")
    
    ##------------------------------------------------------------ 
   # range_for_figure_info=eval(as.call(exargs))
#    print(range_for_figure_info)
#    range_for_figure_info=unlist(range_for_figure_info) 
#    if(range_for_figure_info[1]==QI_name[1]){temp_QI_Range=QI_range[1,]} 
#    if(range_for_figure_info[1]==QI_name[2]){temp_QI_Range=QI_range[2,]} 
#    if(range_for_figure_info[1]==QI_name[3]){temp_QI_Range=QI_range[3,]} 
#    if(range_for_figure_info[1]==QI_name[4]){temp_QI_Range=QI_range[4,]} 
#                
#      param.list[[1]]$Color_Scale$Min_QI$from=temp_QI_Range[1] 
#      param.list[[1]]$Color_Scale$Min_QI$to=temp_QI_Range[3] 
#      param.list[[1]]$Color_Scale$Min_QI$init=temp_QI_Range[2]
#                   
#      param.list[[1]]$Color_Scale$Max_QI$from=temp_QI_Range[3] 
#      param.list[[1]]$Color_Scale$Max_QI$to=temp_QI_Range[5] 
#      param.list[[1]]$Color_Scale$Max_QI$init=temp_QI_Range[4]
     
              fillframe(tt, param.list, plotloc, "tkv")   ##<------------------------
              
    tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc)))

    ##==========================================================================
    ##==========================================================================
    if (wait) {tkwait.window(tt);return(eval(as.call(exargs)))}
    else {return(invisible(NULL))}
}
##------------------Mytkexamp for polygon--------------------------------------- 
Mytkexamp_polygon=function (FUN, param.list, vscale = 1.5, hscale = 1.5, wait = FALSE,plotloc = "top",Tkexample_title)  #plottype=1 heatmap; plottype=2 polygon 
{
    #FUN=HeatMap; param.list=HeatMaplist;vscale = 1.5; hscale = 2.75;Tkexample_title="QI HeatMap plot";plotloc = "top"; 
    #FUN=Polygon; param.list=Polygonlist;vscale = 1.5; hscale = 2.75;Tkexample_title="QI Polygon plot";plotloc = "top"
    if (!require("tkrplot")) {stop("The tkrplot package is needed")}
    require(tcltk2)
    ocl <- cl <- substitute(FUN)
    exargs <- as.list(quote(list()))
    replot <- function() {eval(cl)}
    tt <- tktoplevel()
    tkwm.title(tt, Tkexample_title)
    img <- tkrplot(tt, replot, vscale = vscale, hscale = hscale)
    tkpack(img, side = plotloc)
    hsc <- tclVar()
    tclvalue(hsc) <- hscale
    vsc <- tclVar()
    tclvalue(vsc) <- vscale
    ##--------------------------------------------------------------------------
    ##                                Main area
    ##--------------------------------------------------------------------------
    fillframe <- function(frame, lst, pkdir, prfx)
    {
        #frame=tt;lst=param.list;prfx="tkv";pkdir=plotloc
        for (i in seq_along(lst))
        {
            #i=1
            vname <- paste(prfx, ".", i, sep = "")
            el <- lst[[i]]
            eln <- names(lst)[i]
            if (is.list(el[[1]]))
            {
                fr <- tkframe(frame, relief = "ridge", borderwidth = 3)
                tkpack(fr, side = pkdir)
                if (length(eln) && nchar(eln))
                {
                  tkpack(tklabel(fr, text = eln), side = "top",anchor = "nw")
                }
                Recall(fr, el, ifelse(pkdir == "top", "left","top"), vname)
                next
            }
            ##--------------------------- slider -------------------------------
            if (tolower(el[[1]]) == "slider")
            {
                tkpack(fr <- tkframe(frame), side = pkdir)
                tkpack(tklabel(fr, text = eln), side = "left",
                  anchor = "s", pady = 4)
                tmp <- tclVar()
                tclvalue(tmp) <- if ("init" %in% names(el)) {
                  el$init
                }
                else if ("from" %in% names(el)) {
                  el$from
                }
                else {
                  1
                }
                alist <- list(fr, variable = tmp, orient = "horizontal",
                  command = function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)),
                    vscale = as.numeric(tclvalue(vsc))))
                el2 <- el[-1]
                el2$init <- NULL
                alist <- c(alist, el2)
                tkpack(do.call("tkscale", alist), side = pkdir)
                tmpcl <- as.list(cl)
                tmpl <- list(substitute(as.numeric(tclvalue(VNAME)),
                  list(VNAME = as.character(tmp))))
                names(tmpl) <- eln
                cl <<- as.call(c(tmpcl, tmpl))
                exargs <<- c(exargs, tmpl)
                next
            }
            ##--------------------------- checkbox ------------------------------
            if (tolower(el[[1]]) == "checkbox")
            {
                tmp <- tclVar()
                tclvalue(tmp) <- if ("init" %in% names(el)) {el$init}
                else {"F"}                                            ##"T"           "F"
                alist <- list(fr, variable = tmp, text = eln,onvalue = "T", offvalue = "F", command = function(...) tkrreplot(img,hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc))))
                el2 <- el[-1]
                el2$init <- NULL
                tmpvars <- if ("values" %in% names(el)) {el$values}
                else {""}
                el2$values <- NULL
                alist <- c(alist, el2)
                tkpack(do.call("tkcheckbutton", alist), side = pkdir)
                tmpcl <- as.list(cl)
                tmpl <- list(substitute(as.logical(tclvalue(VNAME)),list(VNAME = as.character(tmp))))
                names(tmpl) <- eln
                cl <<- as.call(c(tmpcl, tmpl))
                exargs <<- c(exargs, tmpl)
                next
            }
            ###==================================================================
            ##--------------------------- radiobuttons -------------------------
            if (tolower(el[[1]]) == "radiobuttons")
            {
                tkpack(fr <- tkframe(frame, relief = "groove",
                  borderwidth = 3), side = pkdir)
                tkpack(tklabel(fr, text = eln), side = "top",anchor = "nw")#,font=fontHeading
                tmp <- tclVar()
                tclvalue(tmp) <- if ("init" %in% names(el)) {
                  el$init
                }
                else {
                  el$values[1]
                }
                el2 <- el[-1]
                tmp.vals <- el2$values
                el2$values <- NULL
                el2$init <- NULL
                alist <- list(fr, variable = tmp, command = function() tkrreplot(img,
                  hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc))))
                pkdir2 <- ifelse(pkdir == "top", "left", "top")
                for (v in tmp.vals) {
                  tkpack(do.call("tkradiobutton", c(alist, value = v,
                    text = v)), side = pkdir2)
                }
                tmpcl <- as.list(cl)
                tmpl <- list(substitute(tclvalue(VNAME), list(VNAME = as.character(tmp))))
                names(tmpl) <- eln
                cl <<- as.call(c(tmpcl, tmpl))
                exargs <<- c(exargs, tmpl)
                next
            }
        }
        #figure_info=(eval(as.call(exargs)))
        #print(figure_info)
        #cat("---------------------------------- \n")
    }
    ##--------------------------------------------------------------------------
    tkpack(tfr <- tkframe(tt), side = "bottom", fill = "x")

    tkpack(tkbutton(tfr, text = "Refresh",width=20, command = function() {
        tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc)))
    }), side = "left", anchor = "s")
    tkpack(tklabel(tfr, text = "      ",width=50), side = "left")
    ##---------modified here for figure print!!-----------------
    tkpack(tkbutton(tfr, text = "Save file",width=20,command = function() {
        tmp <- c(as.list(ocl), eval(as.call(exargs)))
        #cat(deparse(as.call(tmp)), "\n")
        figure_info=(eval(as.call(exargs)))
        #print("figure_info");print(figure_info)
        if(Tkexample_title=="QI HeatMap plot"){HeatMap_output(figure_info,output_popu);cat("HeatMap figure print out.\n")}
        if(Tkexample_title=="QI Polygon plot"){Polygon_output(figure_info,output_popu);cat("Polygon figure print out.\n")}
        #print("------------")
        flush.console()
    }), side = "left", anchor = "s")
    ##===========================================================
    tkpack(tkbutton(tfr, text = "Exit",width=20, command = function() tkdestroy(tt)), side = "right", anchor = "s")
    tkpack(tfr <- tkframe(tt), side = "bottom", fill = "x")
    tkpack(tklabel(tfr, text = "Hscale: "), side = "left")
    tkpack(tkentry(tfr, textvariable = hsc, width = 6), side = "left")
    tkpack(tklabel(tfr, text = "      Vscale: "), side = "left")
    tkpack(tkentry(tfr, textvariable = vsc, width = 6), side = "left")
              fillframe(tt, param.list, plotloc, "tkv")   ##<------------------------
    tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc)))

    ##==========================================================================
    ##==========================================================================
    if (wait) {tkwait.window(tt);return(eval(as.call(exargs)))}
    else {return(invisible(NULL))}
}

##------------------------------------------------------------------------------
##                           HeatMap figure
##------------------------------------------------------------------------------
HeatMap_output=function(A,output_popu)
{
	#A=figure_info
	QIstat=as.character(A$QI_Selection)
	ChipSelection=as.character(A$Sort_By_Array)
	MinQI=A$Min_QI;MaxQI=A$Max_QI
	Ncolor=A$N_color
	##-------------------------------------------------------------
	cols <-terrain.colors(Ncolor)
	#if(Reset_option==T){MinQI=Zmin;MaxQI=Zmax}
	#print(QIstat);print(QI_name)
	if(QIstat==QI_name[1]){QI=MQI1}
	if(QIstat==QI_name[2]){QI=MQI2}
	if(QIstat==QI_name[3]){QI=WQI1}
	if(QIstat==QI_name[4]){QI=WQI2}
	# iRank<- switch(ChipSelect,2="Merge",3="Chip2",4="Chip1")
	if(ChipSelection==chip_name[1]){iRank=4}
	if(ChipSelection==chip_name[2]){iRank=3}
	if(ChipSelection==chip_name[3]){iRank=2}
	z=QI[order(QI[,iRank],decreasing=T),]
	indname=z[,1]
	z=z[,-1];z=as.matrix(z);z=matrix(as.numeric(z),nrow(z),ncol(z))
	nind=nrow(z)
	y.lim = if(index.2chip){
		3
	} else 1
	y=c(1:3)
	x=c(1:nind)
	tempz=z
	tempz[tempz<MinQI]=MinQI;tempz[tempz>MaxQI]=MaxQI
	QIstat=strsplit(QIstat,"  ");QIstat=unlist(QIstat);QIstat=QIstat[1]
	HeatMap_name=paste("HeatMap-",QIstat,"-(",round(MinQI*100)/100,",",round(MaxQI*100)/100,",",Ncolor,")",sep="")
	HeatMap_figure=paste(output_popu,"/",HeatMap_name,".png",sep="")
	png(HeatMap_figure,width = 1440, height = 800,pointsize=14)
	par(mai=c(0.75,0.35,0.25,0.35))
	##==========================  HeatMap figure ============================
	##==========================  HeatMap ===================================
	tempz=z
	tempz[tempz<MinQI]=MinQI;tempz[tempz>MaxQI]=MaxQI
	dspace=0.5
	#image(x,y,tempz, col=cols, ylim=c(0.5,3.5+dspace),xlim=c(0.5,(nind+0.5)),zlim=c(MinQI,MaxQI),axes=F,ylab="",xlab="")
	image(x,y,tempz, col=cols, ylim=c(0.5, y.lim + 0.5 + dspace),xlim=c(0.5,(nind+0.5)),zlim=c(MinQI,MaxQI),axes=F,ylab="",xlab="")
	##-------------------------------------------
	if(nind<61)
	{
		Xlist=(rep(x,3));Ylist=c(rep(1,nind),rep(2,nind),rep(3,nind))
		a=as.numeric(z);a=round(a*100)/100
		text(Xlist,Ylist,labels=a,srt = 60,col="blue",cex=0.75)
	}
	##-------------------------------------------
	block=8;n=round(nind/block);Xlab=c(1:n)*block
	Xlab[length(Xlab)]=nind;
	x=c(0,x);y=c(0,y)
	abline(v=c(x+0.5),lty=2,col="skyblue");
	abline(h=c(y+0.5),lty=1,col="blue");
	#axis(1,at=Xlab,labels=Xlab,tick=F,font=4,line=-0.5)
	#axis(1,at=median(x),labels="Sample",tick=F,font=4,line=1.5)
	abline(v=c(Xlab+0.5),lty=1,col="blue");
	#axis(2,at=c(1,2,3),labels=c("Chip1","Chip2","Merge"),tick=F,font=4,line=-0.5,col.axis="blue")
	axis(2,at=1:y.lim,labels=if(index.2chip){ c("Chip1","Chip2","Merge") } else unlist(strsplit(ChipSelection, " "))[1],
			tick=F,font=4,line=-0.5,col.axis="blue")
	x=x[-1]
	axis(1,at=x,labels=indname,col.axis="blue",tick=F,las=2,line=-2.25,cex.axis=0.75,font=4)
	#rect(0.5,3.5,nind+0.5,3.5+dspace,col="aliceblue")      #gray75
	rect(0.5, y.lim + 0.5,nind+0.5, y.lim + 0.5 + dspace,col="aliceblue")      #gray75
	xleft=0.5+nind*0.25;xright=0.5+nind*0.75
	#ybottom=3.5+dspace*0.375;ytop=3.5+dspace*0.75
	ybottom= y.lim + 0.5 + dspace*0.375;ytop = y.lim + 0.5 + dspace*0.75
	xspace=(xright-xleft)/Ncolor
	xblock=xleft+(0:Ncolor)*xspace
	for(i in 1:Ncolor){rect(xblock[i],ybottom,xblock[i+1],ytop,col=cols[Ncolor+1-i],border=F)}
	Z=range(tempz);Zset=Z[1]+c(0:6)*diff(Z)/6;Zset=round(Zset*1000)/1000;Zset=Zset[7:1]
	U=c(xleft,xright);Uset=U[1]+c(0:6)*diff(U)/6;Uset=round(Uset*1000)/1000;
	#text(Uset,rep(3.5+dspace*0.3,length(Uset)),labels=Zset,col="blue",font=4)
	text(Uset,rep(y.lim + 0.5 + dspace*0.3,length(Uset)),labels=Zset,col="blue",font=4)
	box(col="blue")
	dev.off()
}
##------------------------------------------------------------------------------
      Polygon_output=function(A,output_popu)
      {
        #A=figure_info
                QIstat=as.character(A$QI_Selection)
                ChipSelection=as.character(A$Sort_By_Array)
                MergeQIline=A$Merge_QI; #print(MergeQIline)
                Chip2QIline=A$Chip2_QI; #print(Chip2QIline)
                Chip1QIline=A$Chip1_QI; #print(Chip1QIline)
                #ChipSelection=strsplit(ChipSelection," ");ChipSelection=unlist(ChipSelection);ChipSelection=ChipSelection[1]
                title_QIstat=strsplit(QIstat,"  ");title_QIstat=unlist(title_QIstat);title_QIstat=title_QIstat[1]
                Polygon_name=paste("Polygon-",title_QIstat,"-(",round(MergeQIline*100)/100,",",round(Chip2QIline*100)/100,",",round(Chip1QIline*100)/100,")",sep="")
               Polygon_figure=paste(output_popu,"/",Polygon_name,".png",sep="")
               png(Polygon_figure,width = 1440, height = 800,pointsize=14)
                Polygon(QIstat,ChipSelection,Chip1QIline,Chip2QIline,MergeQIline)
               dev.off()
      }
      