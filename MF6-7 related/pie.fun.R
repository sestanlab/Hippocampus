library(scatterpie)
library(ggforce)
library(tidyr)
geom_scatterpie_new <- function(mapping = NULL, data, cols, legend_name = "species", scale.expression = TRUE, ...) {
      names(mapping)[match(c("x", "y"), names(mapping))] <- c("x0", "y0")
      yvar <- get_aes_var(mapping, "y0")
      xvar <- get_aes_var(mapping, "x0")
      
      
      col_use <- names(cols)
      df <- data[rowSums(data[, col_use]) > 0, ] %>%
               gather(., "type", "value", !!enquo(col_use)) %>% 
               group_by(!!sym(xvar), !!sym(yvar)) %>%
               mutate(scaleexp = value/mean(value)) %>%
               mutate(scaleexp = MinMax(scaleexp, min = 0.5, max = 2)) %>% 
               ungroup()


      nbks <- 1000
      df$border_color <- cols[as.character(df$type)]      
      ## set the border & fill color:
      if (scale.expression){
         df$fill_color <- mapply(FUN = function(color, value) {
                           return(colorRampPalette(colors = c("grey95", color))(nbks)[value])
                  }, color = df$border_color, value = as.numeric(cut(df$scaleexp, nbks)))
         df$border_color[df$scaleexp < 1.5] <- NA
      } else{
         df$fill_color <- df$border_color
         df$border_color <- NA
      }

      
      ## Factorize the type column
      df$type <- factor(df$type, levels = col_use)
      names(df)[which(names(df) == "type")] = legend_name
      df$r0 <- 0
      return(geom_arc_bar(mapping, data = df, stat = "pie", inherit.aes = FALSE, ...))
}
environment(geom_scatterpie_new) <- environment(geom_scatterpie)



PlotScatterPieInteraction <- function(pie.data, group.col = "cluster", feature.col = "gene", r.col = "radius", split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), group.order = NULL, feature.order = NULL, rsf = 2, scale.expression = TRUE, bg.col = "lightgrey"){ ## rsf:radius-scale-factor
   ## set the colors 
   all.sp <- c("mett7b", "coexp", "other", "noexp")
   sp.colors <- setNames(c("#d73027", "#ffffbf", "#1a9850", "#e5e5e5"), c(all.sp)) %>%   
            .[split.order]


   ## Set the order of x & y Axis
   if (is.null(group.order)){
      group.order <- levels(as.factor(pie.data[, group.col]))
   }
   if (is.null(feature.order)){
      feature.order <- levels(as.factor(pie.data[, feature.col]))
   }


   ## pie data
   pie.plot <- pie.data %>%
               mutate(!!group.col := as.character(!!sym(group.col))) %>%
               mutate(!!feature.col := as.character(!!sym(feature.col))) %>%
               mutate(!!group.col := as.numeric(factor(!!sym(group.col), levels = group.order))) %>%
               mutate(!!feature.col := as.numeric(factor(!!sym(feature.col), levels = rev(feature.order))))

   ## legend data
   lg.ra <- pie.plot[, r.col]
   lg.ra[lg.ra == 0] <- 0.00001


   ## Plots
   p <- ggplot(data = pie.plot) + 
      geom_scatterpie_new(mapping = aes_string(x = group.col, y = feature.col, r0 = "r0", r = r.col, amount = "value", fill = "fill_color", color = "border_color"), data = pie.plot, cols = sp.colors, scale.expression = scale.expression) + 
      theme_cowplot() + 
      coord_fixed() + 
      scale_x_continuous(breaks = 1:length(group.order), labels = setNames(group.order, as.character(1:length(group.order)))) +
      scale_y_continuous(breaks = 1:length(feature.order), labels = setNames(rev(feature.order), as.character(1:length(feature.order)))) +
      xlim(0, length(group.order )+ 1) +
      scale_fill_identity() + 
      scale_color_identity() + 
      RotatedAxis() + 
      theme(panel.grid.major = element_line(colour=bg.col, size = 0.05), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_text(size = 2))#+
      ##geom_scatterpie_legend(lg.ra, x=length(group.order)+2, y=ceiling(length(feature.order)/4), n = 4, labeller = function(x) rsf * x)
   return(p)
}

 
PlotScatterPie <- function(pie.data, group.col = "cluster", feature.col = "gene", r.col = "radius", pie_scale = 0.5, split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), group.order = NULL, feature.order = NULL, rsf = 2, cluster.coord){ ## rsf:radius-scale-factor
   ## set the colors 
   all.sp <- c("mett7b", "coexp", "other", "noexp")
   sp.colors <- setNames(c("#d73027", "#ffffbf", "#1a9850", "#e5e5e5"), c(all.sp)) %>%   
            .[split.order]


   ## Set the order of x & y Axis
   if (is.null(group.order)){
      group.order <- levels(as.factor(pie.data[, group.col]))
   }
   if (is.null(feature.order)){
      feature.order <- levels(as.factor(pie.data[, feature.col]))
   }


   pie.plot <- pie.data %>%
               mutate(!!group.col := as.character(!!sym(group.col))) %>%
               mutate(!!feature.col := as.character(!!sym(feature.col))) %>%
               ##mutate(!!group.col := as.numeric(factor(!!sym(group.col), levels = group.order))) %>%
               mutate(!!feature.col := as.numeric(factor(!!sym(feature.col), levels = rev(feature.order)))) %>%
               mutate(!!group.col := cluster.coord[!!sym(group.col)])


   lg.ra <- pie.plot[, r.col]
   lg.ra[lg.ra == 0] <- 0.00001

   p <- ggplot(data = pie.plot) + 
      geom_scatterpie(mapping = aes_string(x = group.col, y = feature.col, r = r.col), data = pie.plot, pie_scale = pie_scale, cols = split.order, color = NA) + 
      theme_cowplot() + 
      coord_fixed() + 
      ##scale_x_continuous(breaks = 1:length(group.order), labels = setNames(group.order, as.character(1:length(group.order)))) +
      scale_x_continuous(breaks = setNames(cluster.coord, NULL), labels = setNames(names(cluster.coord), as.character(setNames(cluster.coord, NULL)))) +
      scale_y_continuous(breaks = 1:length(feature.order), labels = setNames(rev(feature.order), as.character(1:length(feature.order)))) +
      scale_fill_manual(values = sp.colors, labels = c("METTL7B", "coexpression", "other genes", "no expression") %>% setNames(., c(c("mett7b", "coexp", "other", "noexp"))))+
      RotatedAxis() + 
      theme(panel.grid.major = element_line(colour="lightgrey", size = 0.2), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())##+
      ##geom_scatterpie_legend(lg.ra, x=length(group.order)+2, y=ceiling(length(feature.order)/4), n = 4, labeller = function(x) rsf * x)
   return(p)
}
##pie.plot[, r.col]





