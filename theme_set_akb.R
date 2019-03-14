require(ggplot2)

# @ R.T. Nakano
theme_RTN <- theme_bw() +
    theme(text = element_text(size = 30, hjust = 0.3, vjust = 0.3),
          panel.border=element_blank(),
          panel.grid=element_blank(),
          legend.key=element_blank(),
          legend.text.align=0,
          legend.position="top",
          strip.text=element_text(face="bold"),
          axis.line=element_line(),
          axis.line.x=element_line(),
          axis.line.y=element_line(),
          panel.background=element_rect(fill="transparent", colour=NA),
          plot.background=element_rect(fill="transparent", colour=NA),
          strip.background=element_rect(fill="transparent", colour=NA),
          strip.placement="outside")

# Colour Palette
col.pal <- c(c_dark_grey,  "indianred4", c_green, c_cudo_skyblue, c_dark_red, c_very_dark_green, 
             c_red,  c_blue,  "darkseagreen4","dodgerblue1", 
             c_dark_brown, "forestgreen", "dodgerblue4","goldenrod4", "skyblue3",
             "deeppink3", "chocolate3", "limegreen", "orangered3", "brown3", "darkorchid3", "khaki3", "sienna3", 
             "slateblue3", "palegreen3", "salmon3", "olivedrab3", "plum3", "thistle3", "royalblue3"
)
