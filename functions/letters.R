letters <- function(x, y, label, cex = 2, font = 2){
  oldmar <- par()
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA, xlim = c(0,100), ylim = c(0,100), xlab = "", ylab = "", xaxs = "i", yaxs = "i", bty = "n", xaxt = "n", yaxt = "n")
  text(label, x = x, y = y, cex = cex, font = font)
  par(new = TRUE, mar = oldmar$mar)
}
