printMeta <- function(plotname, pageno, jobname) {
    
    gp=gpar(fontsize=11, fontfamily="Helvetica", col="black", fontface="bold")
    gpfill=gpar(fill="gray90", lwd=0, lty=0)
    grid.rect(vp=viewport(layout.pos.row=1, layout.pos.col=1), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=1, layout.pos.col=2), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=1, layout.pos.col=3), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=1, layout.pos.col=4), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=1, layout.pos.col=5), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=1, layout.pos.col=6), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=5, layout.pos.col=1), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=5, layout.pos.col=2), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=5, layout.pos.col=3), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=5, layout.pos.col=4), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=5, layout.pos.col=5), gp=gpfill)
    grid.rect(vp=viewport(layout.pos.row=5, layout.pos.col=6), gp=gpfill)
    
    grid.text(plotname, vp=viewport(layout.pos.row=1, layout.pos.col=1), 
              just=c("left", "center"), 
              gp=gpar(fontface="bold", fontsize=14, fontfamily="Helvetica", col="darkBlue"))
    grid.text("Normalyzer Report", vp=viewport(layout.pos.row=1, layout.pos.col=6),
              just=c("right", "center"), gp=gp)
    grid.text(paste("Page ", pageno,sep=""), vp=viewport(layout.pos.row=5, layout.pos.col=6),
              just=c("right", "center"), gp=gp)
    grid.text(paste("Project: ", jobname, sep=""), 
              vp=viewport(layout.pos.row=5, layout.pos.col=1), just=c("left", "center"), gp=gp)
}








