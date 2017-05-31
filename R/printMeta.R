#' Print meta information for Normalyzer plot page
#' ! Needs refactoring to reduce redundancy in code
#' ! Needs double check of functionality
#' 
#' @param plotname Name of current plot.
#' @param pageno Current page number.
#' @param jobname Name of ongoing job.
#' @return None
printMeta <- function(plotname, pageno, jobname, currentLayout) {

    ncol <- currentLayout$ncol
    nrow <- currentLayout$nrow
        
    gp <- grid::gpar(fontsize=11, fontfamily="Helvetica", 
                     col="black", fontface="bold")
    gpfill <- grid::gpar(fill="gray90", lwd=0, lty=0)
    
    for (pos in 1:(ncol-1)*(nrow-1)) {
        if (pos <= ncol) row <- 1
        else row <- nrow
        
        col <- (pos - 1) %% ncol + 1
        grid::grid.rect(vp=grid::viewport(layout.pos.row=row, 
                                          layout.pos.col=col), gp=gpfill)
    }
    
    grid::grid.text(plotname, 
                    vp=grid::viewport(layout.pos.row=1, layout.pos.col=1), 
                    just=c("left", "center"), 
                    gp=grid::gpar(fontface="bold", fontsize=14, fontfamily="Helvetica", col="darkBlue"))
    grid::grid.text("Normalyzer Report", 
                    vp=grid::viewport(layout.pos.row=1, layout.pos.col=ncol),
                    just=c("right", "center"), gp=gp)
    grid::grid.text(paste("Page ", pageno,sep=""), 
                    vp=grid::viewport(layout.pos.row=nrow, layout.pos.col=ncol),
                    just=c("right", "center"), gp=gp)
    grid::grid.text(paste("Project: ", jobname, sep=""), 
                    vp=grid::viewport(layout.pos.row=nrow, layout.pos.col=1), just=c("left", "center"), gp=gp)
}








