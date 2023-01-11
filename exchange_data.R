
path = "exchange_rate.tsv"

dati_exchange <- read.table(file = path, sep = '\t', header = TRUE)

dati_exchange <- dati_exchange[,-1]

dati_exchange <- t(dati_exchange)

rownames(dati_exchange) <- NULL

colnames(dati_exchange) <- c('EU/CHF','EU/GBP','EU/USD')

dati_exchange <- as.data.frame(dati_exchange)

# Convert data in numeric
dati_exchange <- as.data.frame(sapply(dati_exchange , as.numeric))

# Delate NA objects

dati_exchange <- na.omit(dati_exchange)

# Set in LOG
dati_exchange[,c(1,2,3)] <- -log(dati_exchange[,c(1,2,3)])


