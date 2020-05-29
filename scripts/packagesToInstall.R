install.packages("survminer")
install.packages("survival")
install.packages("ggfortify")
install.packages("ggplot2")
install.packages("gdtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("psichomics")

install.packages("knitr")
install.packages("kableExtra")
install.packages("extrafont")
install.packages("showtext")
library(extrafont)
X11.options()
X11(type = "cairo")
library(ggplot2)
require(extrafont)
extrafont::font_import()  ## first time could take a while
extrafont::loadfonts()
write.table(extrafont::fonttable(), file="fonts.txt", sep="\t")
font_import(pattern="Helvetica")
extrafont::choose_font("FreeSans")
ggplot_font <- extrafont::choose_font(c("Helvetica", "Arial"))  ## 

f = data.frame(one = c(1:5), two = c(6:10))
options(bitmapType="cairo")
ggplot(f, aes(one,two)) + geom_point()
+
  ggtitle("This is a default font")+
  theme(plot.title = element_text(size=30, face="bold", vjust=1))

dev.off()
sessionInfo

ggplot(nmmaps, aes(x = date, y = temp)) + geom_point(color="red")+
  ggtitle("This is a default font")+
  theme(plot.title = element_text(size=30, face="bold", vjust=1))
getwd()

