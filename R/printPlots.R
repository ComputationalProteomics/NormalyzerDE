#' Generate PDF grid page filling it with provided list of plots
#' ! Needs refactoring to reduce redundancy in code
#' 
#' @param plotlist List of target plots to display.
#' @param plotname List of names corresponding to the provided plot list.
#' @param pageno Current page number.
#' @param jobname Name of ongoing job.
printPlots <- function(plotlist, plotname, pageno, jobname) {
    
    print("In print plots!")
    
    gp <- grid::gpar(fontsize=11, fontfamily="Helvetica", col="black", fontface="bold")
    gpfill <- grid::gpar(fill="gray90", lwd=0, lty=0)
    
    grid::grid.rect(vp=grid::viewport(layout.pos.row=1, layout.pos.col=1), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=1, layout.pos.col=2), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=1, layout.pos.col=3), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=1, layout.pos.col=4), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=1, layout.pos.col=5), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=1, layout.pos.col=6), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=5, layout.pos.col=1), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=5, layout.pos.col=2), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=5, layout.pos.col=3), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=5, layout.pos.col=4), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=5, layout.pos.col=5), gp=gpfill)
    grid::grid.rect(vp=grid::viewport(layout.pos.row=5, layout.pos.col=6), gp=gpfill)
    
    grid::grid.text(plotname, 
              vp=grid::viewport(layout.pos.row=1, layout.pos.col=1),
              just=c("left","center"), gp=grid::gpar(fontface="bold", fontsize=14, 
                                               fontfamily="Helvetica", col="darkBlue"))
    grid::grid.text("Normalyzer Report", 
              vp=grid::viewport(layout.pos.row=1, layout.pos.col=6),
              just=c("right", "center"), gp=gp)
    grid::grid.text(paste("Page ", pageno,sep=""), 
              vp=grid::viewport(layout.pos.row=5, layout.pos.col=6), just=c("right", "center"), gp=gp)
    grid::grid.text(paste("Project: ", jobname, sep=""),
              vp=grid::viewport(layout.pos.row=5, layout.pos.col=1), just=c("left", "center"), gp=gp)
    
    print(plotlist[[1]], vp=grid::viewport(layout.pos.row=2, layout.pos.col=2))
    print(plotlist[[2]], vp=grid::viewport(layout.pos.row=2, layout.pos.col=3))
    print(plotlist[[3]], vp=grid::viewport(layout.pos.row=2, layout.pos.col=4))
    print(plotlist[[4]], vp=grid::viewport(layout.pos.row=2, layout.pos.col=5))
    print(plotlist[[5]], vp=grid::viewport(layout.pos.row=3, layout.pos.col=2))
    
    if (length(plotlist) > 5) {
        
        print(plotlist[[6]], vp=grid::viewport(layout.pos.row=3, layout.pos.col=3))
        print(plotlist[[7]], vp=grid::viewport(layout.pos.row=3, layout.pos.col=4))
        print(plotlist[[8]], vp=grid::viewport(layout.pos.row=3, layout.pos.col=5))
        print(plotlist[[9]], vp=grid::viewport(layout.pos.row=4, layout.pos.col=2))
        print(plotlist[[10]], vp=grid::viewport(layout.pos.row=4, layout.pos.col=3))
        print(plotlist[[11]], vp=grid::viewport(layout.pos.row=4, layout.pos.col=4))
        
        if (length(plotlist) == 12) {
            print(plotlist[[12]], vp=grid::viewport(layout.pos.row=4, layout.pos.col=5))  
        }
    }
}