library(readxl)

# ejemplo de los intervalos de confianza
#x=c(45,35,37,39.5,43.2,47.1,38,39.1,36.4,40)

#t.test(x,mu=40)

TIPSenfermos <- read_excel("C:/Users/SIMON/Google Drive/Doctorado/Semillero/HERVs in Cancer/nuevo/TIPsXpaciente.xlsx", sheet="Enfermos")
TIPSsanos <- read_excel("C:/Users/SIMON/Google Drive/Doctorado/Semillero/HERVs in Cancer/nuevo/TIPsXpaciente.xlsx", sheet="Sanos")

windows()
boxplot(TIPSenfermos$Paciente~TIPSenfermos$TIPs, main= "boxplot",horizontal=TRUE)

windows()
boxplot(TIPSsanos$Paciente~TIPSsanos$TIPs, main= "boxplot",horizontal=TRUE)

#t.test(presionSanguinea$Pacientes~presionSanguinea$grupo,alternative="two.sided")

#prueba de normalidad
pvalorEnfermos=shapiro.test(TIPSenfermos$TIPs) #p-valor = 0.1155

#prueba de normalidad
pvalorSanos=shapiro.test(TIPSsanos$TIPs) #p-valor = 0.6172

#analisis de asociaciones de cada uno de los TIPs en los cromosomas
TIPs <- read_excel("C:/Users/SIMON/Google Drive/Doctorado/Semillero/HERVs in Cancer/nuevo/analisisEnfermosVSsanos.xlsx", sheet="Sheet1")
dataframe = data.frame(TIPs)
asociaciones = data.frame()
names(asociaciones)<-c("Insertion","Case", "Control", "P-valor")
for (i in 1:6929){
  casoSi = dataframe[i,2]
  casoNo = 30-dataframe[i,2]
  controlSi = dataframe[i,3]
  controlNo = 30-dataframe[i,3]
  data = matrix(c(casoSi,controlSi,casoNo,controlNo), nrow=2, byrow=TRUE)
  chicuadrado = chisq.test(data)
  medidaAsoc = sqrt(as.double(chicuadrado$statistic)/(as.double(chicuadrado$statistic) + 60))
  pval = chicuadrado$p.value
  if(pval < 0.05 ){
    message(paste0(c(dataframe[i,1], dataframe[i,2], dataframe[i,3], pval, medidaAsoc), " "))
    de<-data.frame(dataframe[i,1], dataframe[i,2], dataframe[i,3], pval, medidaAsoc)
    names(de)<-c("Insertion","Case", "Control", "P-valor", "Asociacion")
    asociaciones <- rbind(asociaciones, de)
  }
}

write.table(asociaciones, file = "C:/Users/SIMON/Google Drive/Doctorado/Semillero/HERVs in Cancer/nuevo/TIPsConAsociaciones.csv", sep = ",", row.names = FALSE)

  