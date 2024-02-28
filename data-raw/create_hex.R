## code to prepare `create_hex` sticker goes here

library(ggplot2)
library(dittoSeq)
library(hexSticker)

basepairs<-c("A","T","G","C")
colors<-dittoColors()[1:4]
df<-data.frame(x=rep(1:40,each=40),
               y=1:40,
               base=sample(basepairs,1600,replace=TRUE),
               color=sample(colors,1600,replace=TRUE))

p<-ggplot(df,aes(x=x,y=y)) + 
  geom_text(label = df$base, size = 6, colour = df$color) + 
  theme(panel.background = element_rect(fill = 'black', color = 'black'),
                                        panel.grid.major = element_blank(), 
                                        panel.grid.minor = element_blank())
print(p)

s <- sticker(p,
             package="demuxSNP", p_size=25, s_x=1, s_y=1, s_width=3, s_height=3, p_y=1,
             filename="inst/figures/baseplot.png")
library(cropcircles)
img_cropped <- hex_crop(
  images = 'inst/figures/baseplot.png',
  border_colour = "black",
  border_size = 24
)

use_logo(img_cropped, geometry = "240x278", retina = TRUE)
#usethis::use_data(create_hex, overwrite = TRUE)
