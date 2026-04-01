#' @title Plotting coefficients
#'
#' @description This function is to plot the maximum likelihood estimates and their respective confidence intervals.
#' @param mles A vector containing the maximum likelihood estimates to be plotted.
#' @param upper A vector containing the upper limits of the confidence intervals.
#' @param lower A vector containing the lower limits of the confidence intervals.
#' @param species_labels A vector containing the labels of each species.
#' @param xlab A label for the x axis.
#' @param cex.axis The size of the x axis.
#' @param y_spa Space to be added to the minimum and maximum values. This is to improve visualization.
#' @param rect Should assign different backgrounds to predators and non-predators? Default is to FALSE.
#' @param rect_lim Limit of the background.
#' @param ... Other graphical parameters.

#' @export

#'



My_coefplot <- function (mles, upper, lower, species_labels = NULL, xlab = NULL, cex.axis = 1,y_spa = 0, rect = F,rect_lim = 5, xlim = NULL, xlimit = 15, species_font = 1, axis_line = 0, ...)
{
  
  col.seq <- rep("grey70", length(mles))
  col.seq[which(lower < 0 & upper < 0)] <- "black"
  col.seq[which(lower > 0 & upper > 0)] <- "black"
  
  #col.seq[which(lower < 0 & upper < 0 & mles < -10)] <- "indianred1"
  #col.seq[which(lower > 0 & upper > 0 & mles > 10)] <- "steelblue1"
  #col.seq[which(upper < 0.0000001 & lower > -0.0000001 & mles > -0.0000001 & mles < 0.0000001)] <- "grey"
  
  
  
  lwd.seq <- rep(1, length(mles))
  lwd.seq[which(lower < 0 & upper < 0)] <- 2
  lwd.seq[which(lower > 0 & upper > 0)] <- 2
  
  #lwd.seq[which(mles < -10 | mles > 10)] <- 1
  
  cex.seq <- rep(0.8, length(mles))
  cex.seq[which(lower < 0 & upper < 0)] <- 1
  cex.seq[which(lower > 0 & upper > 0)] <- 1
  
  #cex.seq[which(mles < -10 | mles > 10)] <- 1
  
  
  At.y <- rev(1:length(mles))
  ylim <- c(min(At.y-y_spa), max(At.y+y_spa))
  
  if(is.null(xlim)){
    xlim <- c(min(lower), max(upper))
    
    if(xlim[1] < -xlimit){
      lower2 <- lower
      lower2[lower2 < -xlimit] <- NA
      xlim[1] <- min(c(mles,lower2), na.rm = TRUE)*1.1
    }
    if(xlim[2] > xlimit){
      upper2 <- upper
      upper2[upper2 > xlimit] <- NA
      xlim[2] <- max(c(mles,upper2), na.rm = TRUE)*1.1    
    }
  }
  
  
  plot(x = NULL, y = NULL, yaxt = "n", xaxt = "n", ylim = ylim,
       ylab = "", xlab = xlab, xlim = xlim, ...)
  
  if(isTRUE(rect)){
    rect(xleft = min(lower)-1, ybottom = 0, xright =  max(upper)+1, ytop = max(At.y)-rect_lim+0.5, density = NULL, border = "transparent", col = rgb(col2rgb("forestgreen", alpha = FALSE)[1],
                                                                                                                                                     col2rgb("forestgreen", alpha = FALSE)[2],
                                                                                                                                                     col2rgb("forestgreen", alpha = FALSE)[3],
                                                                                                                                                     alpha = 40, maxColorValue = 255))
    rect(xleft = min(lower)-1, ybottom = max(At.y)-rect_lim+0.5, xright =  max(upper)+1, ytop = max(At.y)+1, density = NULL, border = "transparent", col = rgb(col2rgb("gold3", alpha = FALSE)[1],
                                                                                                                                                               col2rgb("gold3", alpha = FALSE)[2],
                                                                                                                                                               col2rgb("gold3", alpha = FALSE)[3],
                                                                                                                                                               alpha = 40, maxColorValue = 255))
  }
  
  points(x = mles, y = At.y, col = col.seq, pch = 16, cex = cex.seq)
  #segments(x0 = lower,
  #         x1 = upper,
  #         y1 = At.y, y0 = At.y, col = col.seq, lwd =lwd.seq)
  
  
  arrows(y1 = At.y, y0 = At.y, x1 = upper, x0 = lower,
         code = 3, angle = 90, length = 0.025,col = col.seq, lwd =lwd.seq)
  
  
  abline(v = 0, lty = 3)
  axis(2, at = At.y, labels = FALSE, tcl = -0.25)
  for(sp in 1:length(sp_font)){
    axis(2, at = At.y[sp], labels = species_labels[sp], tick = FALSE, las = 1, cex.axis = cex.axis, font = species_font[sp], line = axis_line, gap.axis = -10)
  }
  
  axis(1, labels = FALSE, tcl = -0.25)
  axis(1, cex.axis = cex.axis, line = axis_line, gap.axis = -10, tick = FALSE)
  
}
