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