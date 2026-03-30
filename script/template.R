# install.packages("tidyverse")
library(tidyverse)
read.xlsx("filename.xlsx")
X2sample_CoJo$n<-as.factor(X2sample_CoJo$n)
ggplot(data = X2sample_CoJo, aes(x = n , y = m ,fill = group))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_hline(yintercept = 3, color = "black")+
  labs(x="number of sample", y="-log(P value)")+
  scale_fill_manual(values = c( "#FFC900", "#FF9000", "#22D222", "#005800", "#D879FF", "#7100A4"))

# jpgs/pngs to pdf
library(magick)
img_dir <- "/Users/annabethlu/jpgs"
files <- list.files(img_dir, pattern = "\\.(jpe?g|png)$", full.names = TRUE, ignore.case = TRUE)

# (optional) sort in a natural way if filenames include numbers like 1,2,10
# install.packages("gtools")
# files <- gtools::mixedsort(files)

imgs <- image_read(files) |> image_convert(format = "pdf", quality = 60)
image_write(imgs, path = "/Users/annabethlu/jpgs/output.pdf", format = "pdf")

# pdf to png
pdf_dir <- 'figs/fig2.pdf'
pdf_convert(pdf_dir, format = "png", dpi = 300)

# png to pdf (high resolution)
img <- image_read("input.png")
img_a4 <- image_resize(img, "5000")
image_write(img_a4, path = "output.pdf", format = "pdf", density = 300)

# compress pics.pdf
input <- image_read_pdf("/Users/annabethlu/jpgs/1_handout.pdf")
out <- image_convert(input, format = "pdf", quality = 60)
image_write(image_join(out), path = "/Users/annabethlu/jpgs/smaller_file.pdf")

# merge pdfs
library(qpdf)
pdf_a <- '/Users/annabethlu/jpgs/1_handout.pdf' 
pdf_b <- '/Users/annabethlu/jpgs/output.pdf'
pdf_combine(input = c(pdf_a, pdf_b), output = "merged.pdf")


