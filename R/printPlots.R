#' Generate PDF grid page filling it with provided list of plots
#' ! Needs refactoring to reduce redundancy in code
#' 
#' @param plotlist List of target plots to display.
#' @param plotname List of names corresponding to the provided plot list.
#' @param pageno Current page number.
#' @param jobname Name of ongoing job.
#' @return None
printPlots <- function(plotlist, plotname, pageno, jobname) {
    
    print("In print plots!")
    
    gp <- grid::gpar(fontsize=11, fontfamily="Helvetica", col="black", fontface="bold")
    gpfill <- grid::gpar(fill="gray90", lwd=0, lty=0)
    
    for (pos in 1:12) {
        if (pos <= 6) row <- 1
        else row <- 5

        col <- (pos - 1) %% 6 + 1
        # print(paste("row", row, "col", col))
        grid::grid.rect(vp=grid::viewport(layout.pos.row=row, layout.pos.col=col), gp=gpfill)
    }
    
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
    
    for (i in 1:length(plotlist)) {
        row <- (i - 1) %/% 4 + 2
        col <- (i - 1) %% 4 + 2
        
        print(plotlist[[i]], vp=grid::viewport(layout.pos.row=row, layout.pos.col=col))
    }
}