#' Generate PDF grid page filling it with provided list of plots
#' 
#' @param plotlist List of target plots to display.
#' @param plotname List of names corresponding to the provided plot list.
#' @param pageno Current page number.
#' @param jobname Name of ongoing job.
#' @param currentLayout Custom viewport layout.
#' @return None
#' @keywords internal
printPlots <- function(plotlist, plotname, pageno, jobname, currentLayout) {
    
    ncol <- currentLayout$ncol
    nrow <- currentLayout$nrow
    
    gp <- grid::gpar(fontsize=11, fontfamily="Helvetica", col="black", fontface="bold")
    gpfill <- grid::gpar(fill="gray90", lwd=0, lty=0)
    
    drawRectangles(nrow, ncol, gpfill, gp)
    
    grid::grid.text(plotname, 
                    vp=grid::viewport(layout.pos.row=1, layout.pos.col=1),
                    just=c("left","center"), gp=grid::gpar(fontface="bold", fontsize=14, 
                                                           fontfamily="Helvetica", col="darkBlue"))
    grid::grid.text("Normalyzer Report", 
                    vp=grid::viewport(layout.pos.row=1, layout.pos.col=ncol),
                    just=c("right", "center"), gp=gp)
    grid::grid.text(paste("Page ", pageno, sep=""), 
                    vp=grid::viewport(layout.pos.row=nrow, layout.pos.col=ncol), just=c("right", "center"), gp=gp)
    grid::grid.text(paste("Project: ", jobname, sep=""),
                    vp=grid::viewport(layout.pos.row=nrow, layout.pos.col=1), just=c("left", "center"), gp=gp)
    
    gridRows <- nrow - 2
    gridCols <- ncol - 2
    gridSize <- gridRows * gridCols

    row <- 2  # Start value
    col <- 2  # Start value

    posCounter <- 0
    for (i in seq_len(length(plotlist))) {
        
        if (gridSize == 1 || (i != 1 && i %% gridSize == 1)) {
            print(paste("New page, row:", row, "col:", col))
            grid::grid.newpage()
            grid::pushViewport(grid::viewport(layout=currentLayout))
            
            grid::grid.rect(vp=grid::viewport(layout.pos.row=row, layout.pos.col=col), gp=gpfill)
            posCounter <- 0
        }
        
        posCounter <- posCounter + 1
        
        row <- (posCounter - 1) %/% gridCols + 2
        col <- (posCounter - 1) %% gridCols + 2
        
        print(plotlist[[i]], vp=grid::viewport(layout.pos.row=row, layout.pos.col=col))
    }
}

drawRectangles <- function(nrow, ncol, gpfill, gp) {
    for (pos in seq_len((ncol - 1) * (nrow - 1))) {
        if (pos <= ncol) row <- 1
        else row <- nrow
        
        col <- (pos - 1) %% ncol + 1
        grid::grid.rect(vp=grid::viewport(layout.pos.row=row, layout.pos.col=col), gp=gpfill)
    }
}