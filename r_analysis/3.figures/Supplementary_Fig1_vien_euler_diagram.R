library(eulerr)
library(RColorBrewer)
brewer.pal(3, "Dark2")

wilkinson <- euler(c("A" = 5, "B" = 7, "C" = 8,
                     "A&B" = 9, "A&C" = 1, "B&C" = 1,
                     "A&B&C" = 3))

plot(wilkinson, quantities = list(type = "counts",font = 3),
     fills = brewer.pal(3, "Dark2"), alpha = 0.4,
     labels = c("Shortwave", "NIR+SWIR", "VIS"))