text_contour <- function(x, y, range_y, range_x, thick = 0.002, ...){
  
  range_y <- par("usr")[3:4]
  range_x <- par("usr")[1:3]

  pos_y <- range_y*thick
  pos_x <- range_x*thick
  
  text(x = x+pos_x, y = y, ...)
  text(x = x-pos_x, y = y, ...)
  text(x = x, y = y+pos_y, ...)
  text(x = x, y = y-pos_y, ...)
  
  text(x = x+pos_x, y = y-pos_y, ...)
  text(x = x-pos_x, y = y+pos_y, ...)
  text(x = x+pos_x, y = y+pos_y, ...)
  text(x = x-pos_x, y = y-pos_y, ...)
  
  text(x = x+(pos_x), y = y-(pos_y/2), ...)
  text(x = x-(pos_x), y = y+(pos_y/2), ...)
  text(x = x+(pos_x), y = y+(pos_y/2), ...)
  text(x = x-(pos_x), y = y-(pos_y/2), ...)
  
  text(x = x+(pos_x/2), y = y-(pos_y), ...)
  text(x = x-(pos_x/2), y = y+(pos_y), ...)
  text(x = x+(pos_x/2), y = y+(pos_y), ...)
  text(x = x-(pos_x/2), y = y-(pos_y), ...)
}
