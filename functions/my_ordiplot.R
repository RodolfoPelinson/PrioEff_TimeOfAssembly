

my_ordiplot <- function(axis = NULL, gllvm_mod = NULL, com = NULL, groups, covariable = NULL,
                        family = "negative.binomial", num.lv = 2, method = "EVA", seed = 3,
                        col = NULL, scale_sites = 1, scale_species = 1, species = TRUE,
                        sites = TRUE, kind = "sd", nmds = FALSE, range = TRUE, axis_name = c("Axis"),
                        zero_lines = TRUE,xlim = NULL, ylim = NULL, mar = c(4,5,0.5,0.5),cex = 1, bty = "o",
                        yline = 2,xline = 2, transparent_pt = NULL, cex.lab = 1, cex.axis = 1, ...){
  
levels <- levels(groups)

if(is.null(col)){
  col <- rep("black", length(levels))
}

if(is.null(axis)){
  if(is.null(gllvm_mod)){
    if(isTRUE(nmds)){
      
      scaled_lvs <- list()
      
      NMDS <- metaMDS(com, distance = "bray", k = 2)
      NMDS$points
      NMDS$species
      scaled_lvs$species <- NMDS$species
      scaled_lvs$sites <- NMDS$points
      
    }else{
      if(is.null(covariable)){
        gllvm_plot <- gllvm(y = com, family = family, num.lv = num.lv, method = method, control.start  = list(n.init = 1000, ...), seed = seed, col = col)
      }else{
        gllvm_plot <- gllvm(y = com, X = data.frame(covariable), family = family, num.lv = num.lv, method = method, control.start  = list(n.init = 1000, ...), seed = seed, col = col)
        
      }
      
      scaled_lvs <- get_scaled_lvs(gllvm_plot)
    }
    
  }else{
    scaled_lvs <- get_scaled_lvs(gllvm_mod)
    
  }
}else{
  scaled_lvs <- list()
  
  scaled_lvs$species <- axis$species
  scaled_lvs$sites <- axis$sites
}








scaled_lvs$sites <- scaled_lvs$sites * scale_sites
scaled_lvs$species <- scaled_lvs$species * scale_species

if(isTRUE(species)){
  xmin <- min(c(scaled_lvs$sites[,1], scaled_lvs$species[,1]))*1.1
  xmax <- max(c(scaled_lvs$sites[,1], scaled_lvs$species[,1]))*1.1
  ymin <- min(c(scaled_lvs$sites[,2], scaled_lvs$species[,2]))*1.1 - 0.005
  ymax <- max(c(scaled_lvs$sites[,2], scaled_lvs$species[,2]))*1.1 + 0.005
}else{
  xmin <- min(c(scaled_lvs$sites[,1]))*1.1
  xmax <- max(c(scaled_lvs$sites[,1]))*1.1
  ymin <- min(c(scaled_lvs$sites[,2]))*1.1 - 0.005
  ymax <- max(c(scaled_lvs$sites[,2]))*1.1 + 0.005
}





par(mar = mar, cex = cex, bty = bty)


if(isTRUE(nmds)){
  ylab <- "NMDS 2"
  xlab <- "NMDS 1"
}else{
  ylab <- "LV 2"
  xlab <- "LV 1"
}

if(is.null(axis)==FALSE){
  ylab <- paste(axis_name,"2")
  xlab <- paste(axis_name,"1")
}

if(is.null(xlim)){
  xlim <- c(xmin, xmax)
}

if(is.null(ylim)){
  ylim <- c(ymin, ymax)
}

par(cex.axis = cex.axis)
plot(NA,xlim = xlim, ylim = ylim, ylab = "", xlab = "", axes = F)
axis(1 , gap.axis = -10)
axis(2 , gap.axis = -10)
title(ylab = ylab, line = yline, cex.lab = cex.lab)
title(xlab = xlab, line = xline, cex.lab = cex.lab)

if(isTRUE(zero_lines)){
  abline(h=0,v=0, lty =2)
} 

if(isTRUE(range)){
  ordiellipse(scaled_lvs$sites,
              groups = groups,
              draw = "polygon", border = col,
              col = "transparent", kind = "ehull", lty = 2,alpha = 0)
}


ordiellipse(scaled_lvs$sites,
            groups = groups,
            draw = "polygon", border = col,
            col = col, kind = kind, alpha = 70)

if(is.null(transparent_pt) == FALSE){
  col.pt <- transparent(col, trans.val = transparent_pt)
}else{
  col.pt <- col
}




if(isTRUE(sites)){
  for(i in 1:length(levels)){
    points(x = scaled_lvs$sites[groups == levels[i],], col = col.pt[i], pch = 16, cex = 1.5)
  }
}



species_labels <- gsub("_", "", rownames(scaled_lvs$species))

if(isTRUE(species)){
  for(i in 1:nrow(scaled_lvs$species)){
    points(x = scaled_lvs$species[i,1], y = scaled_lvs$species[i,2], cex = 4, pch = 21, bg = "white")
    text(x = scaled_lvs$species[i,1], y = scaled_lvs$species[i,2], labels = substr(species_labels[i], 1, 3), cex = 0.9)
  }
}
box()


#points(centroid_controle[1], centroid_controle[2], bg = col_controle, pch = 23, cex = 1.5)
#points(centroid_atrasado[1], centroid_atrasado[2], bg = col_atrasado, pch = 24, cex = 1.5)
#points(centroid_fechado[1], centroid_fechado[2], bg = col_fechado, pch = 22, cex = 1.5)
#points(centroid_fechado_dia[1], centroid_fechado_dia[2], bg = col_fechado_dia, pch = 21, cex = 1.5)
#points(centroid_fechado_noite[1], centroid_fechado_noite[2], bg = col_fechado_noite, pch = 21, cex = 1.5)

#legend(y = ymax, x = xmax,
#       legend = c("Delayed" ,"Control", "Closed", "Night exc.", "Day exc."),
#       pch = 16, col = c(col_atrasado, col_controle, col_fechado, col_fechado_noite, col_fechado_dia),
#       xjust = 1, cex = 1, bty = "o")

#legend(y = ymax, x = xmin,
#       legend = c("Delayed ≠ Night exc",
#                  "Delayed ≠ Day exc.",
#                  "Closed ≠ Night exc"),
#       xjust = 0, cex = 0.75, bty = "n")
}



