letters <- function(x, y, label, cex = 2, font = 2, adj = c(0,1)){
  oldmar <- par()
  par(new = TRUE, mar = c(0,0,0,0))
  plot(NA, xlim = c(0,100), ylim = c(0,100), xlab = "", ylab = "", xaxs = "i", yaxs = "i", bty = "n", xaxt = "n", yaxt = "n")
  text(label, x = x, y = y, cex = cex, font = font, adj = adj)
  par(new = oldmar$new, mar = oldmar$mar)
}
