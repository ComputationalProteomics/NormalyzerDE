#' Generate PDF grid page filling it with provided list of plots
#' 
#' @param plotlist List of target plots to display.
#' @param plotname List of names corresponding to the provided plot list.
#' @param pageno Current page number.
#' @param jobname Name of ongoing job.
#' @return None
printPlots <- function(plotlist, plotname, pageno, jobname, currentLayout) {
    
    print("In print plots!")
    
    ncol <- currentLayout$ncol
    nrow <- currentLayout$nrow
    
    gp <- grid::gpar(fontsize=11, fontfamily="Helvetica", col="black", fontface="bold")
    gpfill <- grid::gpar(fill="gray90", lwd=0, lty=0)
    
    draw_rectangles(nrow, ncol, gpfill, gp)
    
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
    
    grid_rows <- nrow - 2
    grid_cols <- ncol - 2
    grid_size <- grid_rows * grid_cols

    row <- 2  # Start value
    col <- 2  # Start value
        
    pos_counter <- 0
    for (i in 1:length(plotlist)) {
        
        # browser()
        
        # if (i != 1 && i %% grid_size == 1) {
        if (grid_size == 1 || (i != 1 && i %% grid_size == 1)) {
            print(paste("New page, row:", row, "col:", col))
            grid::grid.newpage()
            grid::pushViewport(grid::viewport(layout=currentLayout))
            
            grid::grid.rect(vp=grid::viewport(layout.pos.row=row, layout.pos.col=col), gp=gpfill)
            pos_counter <- 0
        }
        
        pos_counter <- pos_counter + 1
        
        row <- (pos_counter-1) %/% grid_cols + 2
        col <- (pos_counter-1) %% grid_cols + 2

        # print(paste("i", i, "pos_counter", pos_counter, "row", row, "col", col, "grid cols", grid_cols, "grid rows", grid_rows))
        
        print(plotlist[[i]], vp=grid::viewport(layout.pos.row=row, layout.pos.col=col))
    }
}

draw_rectangles <- function(nrow, ncol, gpfill, gp) {
    for (pos in 1:(ncol-1)*(nrow-1)) {
        if (pos <= ncol) row <- 1
        else row <- nrow
        
        col <- (pos - 1) %% ncol + 1
        grid::grid.rect(vp=grid::viewport(layout.pos.row=row, layout.pos.col=col), gp=gpfill)
    }
}