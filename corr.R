library(data.table)
library(corrplot)

data <- fread("langevin_corr.dat")

res <- cor(t(data))

corrplot(res, method="number")

