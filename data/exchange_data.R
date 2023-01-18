path = ".\\data\\exchange_rate.tsv"

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
dati_exchange[,c(1,2,3)] <- log(dati_exchange[,c(1,2,3)])
dati_exchange = dati_exchange*100

# Divided data
dati_CHF = dati_exchange[,1]
dati_GBP = dati_exchange[,2]
dati_USD = dati_exchange[,3]

# Plot the data
plot(dati_CHF,type = "l", lwd = 2, col = "darkred", xlab="Time", ylab="Exchange rate",
     ylim = c(min(dati_exchange),max(dati_exchange)))
grid()
lines(dati_GBP, lwd = 2, col = "darkblue")
lines(dati_USD, lwd = 2, col = "darkgreen")
legend(0,0, legend = c("EU/CHF", "EU/GBP", "EU/USD"), lwd = 3, 
       cex=0.42, col = c("darkred","darkblue","darkgreen"))

dati_exchange = as.data.frame(cbind(dati_CHF, dati_USD))



