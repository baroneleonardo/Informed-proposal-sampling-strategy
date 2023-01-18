
path = ".\\data\\dpc-covid19-ita-regioni.csv"

dati_covid_regioni <- read.csv(file = path)

# 

nord      <- c("Emilia-Romagna","Friuli Venezia Giulia", "Liguria", "Lombardia", "P.A. Bolzano", "P.A. Trento", "Piemonte","Valle d'Aosta","Veneto")

centro    <- c("Lazio","Toscana","Umbria","Marche")

sud_isole <- c("Abruzzo","Basilicata","Calabria","Campania","Molise","Puglia","Sardegna","Sicilia")

for (i in 1:nrow(dati_covid_regioni)){
  
  if(dati_covid_regioni[i,"denominazione_regione"] %in% nord){
    dati_covid_regioni[i,"zona"] = "Nord"}
  
  if(dati_covid_regioni[i,"denominazione_regione"] %in% centro){
    dati_covid_regioni[i,"zona"] = "Centro"}
  
  if(dati_covid_regioni[i,"denominazione_regione"] %in% sud_isole){
    dati_covid_regioni[i,"zona"] = "Sud e Isole"}
  
}

date <- names(table(dati_covid_regioni$data))


casi_tot <- as.data.frame(matrix(NA, nrow = length(date), ncol = 4))

colnames(casi_tot) = c("Date","Nord","Centro","Sud e Isole")

casi_tot$Date <- date

tot_casi_nord       <- as.numeric()
tot_casi_centro     <- as.numeric()
tot_casi_sud_isole  <- as.numeric()

for(data in date){
  
  dati_data <- dati_covid_regioni[which(dati_covid_regioni$data == data),]
  
  row_nord      = which(dati_covid_regioni[which(dati_covid_regioni$data == data), "zona"] == "Nord")
  
  row_centro    = which(dati_covid_regioni[which(dati_covid_regioni$data == data), "zona"] == "Centro")
    
  row_sud_isole = which(dati_covid_regioni[which(dati_covid_regioni$data == data), "zona"] == "Sud e Isole")
  
  casi_tot[which(casi_tot$Date == data), "Nord"]        = sum(dati_data[row_nord,"nuovi_positivi"])
  
  casi_tot[which(casi_tot$Date == data), "Centro"]      = sum(dati_data[row_centro,"nuovi_positivi"])
  
  casi_tot[which(casi_tot$Date == data), "Sud e Isole"] = sum(dati_data[row_sud_isole,"nuovi_positivi"])
  
}
                          

plot(casi_tot$Nord,type = "l")
plot(casi_tot$Centro,type = "l")
plot(casi_tot$`Sud e Isole`,type = "l")

# delete days with 0 new cases in South and Islands
remove <- which(casi_tot$`Sud e Isole` == 0) 
casi_tot = casi_tot[-remove,]
#

# there's a day with -222 new cases in South and Islands, we remove it
remove <- which(casi_tot$`Sud e Isole` < 0) 
casi_tot = casi_tot[-remove,]
#

log_nord   = log(casi_tot$Nord/casi_tot$`Sud e Isole`)
log_centro = log(casi_tot$Centro/casi_tot$`Sud e Isole`)

plot(log_nord,type = "l", lwd = 2, col = "darkred", ylim = c(-1.5,4))
lines(log_centro, lwd = 2)

dati_finali_log <- as.data.frame(cbind(log_nord, log_centro))

