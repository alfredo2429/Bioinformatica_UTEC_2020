# COLORES 
colores1 <- c(rgb(219/255,84/255,97/255,alpha = 1),
              rgb(255/255,200/255,87/255,alpha = 1),
              rgb(87/255,204/255,153/255,alpha = 1),
              rgb(56/255,163/255,165/255,alpha = 1),
              rgb(98/255,76/255,171/255,alpha = 1)
)

colores2 <- c(rgb(219/255,84/255,97/255,alpha = .5),
              rgb(255/255,200/255,87/255,alpha = .5),
              rgb(87/255,204/255,153/255,alpha = .5),
              rgb(56/255,163/255,165/255,alpha = .5),
              rgb(98/255,76/255,171/255,alpha = .5)
)

# Importación de una secuencia desde una secuencia de una base de datos remota
library(seqinr)
library(ape)

# definimos la ruta
ruta <- paste0(getwd(), .Platform$file.sep, "data", .Platform$file.sep, "spike.fasta")
print(ruta)

# creamos objeto que contiene la secuencia
spike <- seqinr::read.fasta(file = paste0(getwd(),.Platform$file.sep, "data",.Platform$file.sep,"spike.txt"), 
                            seqtype = "AA")

# imprimimos en pantalla la secuencia
seqinr::getSequence(spike)

# generamos una tabla de contingencia para ver la frecuencia absoluta de aminoácidos
p1 <- table(seqinr::getSequence(spike))
print(p1)

# Generamos un gráfico de barras
barplot(p1, las = 1, main = "Frecuencia absoluta de AA", ylab = "frecuencia absoluta", xlab = "aminoacido", col = colores1)

# Definimos la longitud
long <- length(unlist(seqinr::getSequence(spike)))
long

# gráfico de barras ordenado
barplot(sort(p1/long, decreasing = T), las = 1, main = "Frecuencia absoluta de AA", ylab = "frecuencia absoluta", xlab = "aminoacido", col = colores1)

# Para el análisis de la estadistica descriptiva de la secuencia
# definimos la secuencia
secuencia <- seqinr::getSequence(spike)

# hacemos el análisis
seqinr::AAstat(unlist(secuencia))



