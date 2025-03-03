# helper functions

publication_theme <- 
  function(axis_title_size = 14, axis_text_size = 12,
           legend_text_size = 12, legend_title_size = 14, 
           strip_text_size = 12, my_font = 'Helvetica',
           grid_colour = 'grey90', background_fill = 'white',
           background_colour = 'white', strip_colour = 'grey60') {
    theme(plot.background = element_rect(fill = background_fill,
                                         colour = background_colour),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = NA, colour = strip_colour,
                                          linewidth = 1),
          strip.text = element_text(colour = 'grey20', size = strip_text_size,
                                    family = my_font),
          axis.line = element_line(colour = 'grey60', linewidth = 0.8),
          axis.ticks = element_line(colour = 'grey60'),
          axis.title = element_text(colour = 'grey20', size = axis_title_size,
                                    family = my_font),
          axis.text = element_text(colour = 'grey20', size = axis_text_size,
                                   family = my_font),
          legend.title = element_text(colour = 'grey20', size = legend_title_size,
                                      family = my_font),
          legend.text = element_text(colour = 'grey20', size = legend_text_size,
                                     family = my_font),
          legend.background = element_blank(),
          legend.key = element_blank(),
          plot.title = element_text(colour = 'grey20', size = legend_title_size,
                                    family = my_font))
  }