setwd("~/Documentos/Conesa_Lab/VSCODE/practica/datos")

# Librerías

library(ggplot2)
library(ggpp)
library(ggpointdensity)
library(viridis)


seq <- "ont"
plataforma <- "bambu"
ruta <-  "~/Documentos/Conesa_Lab/VSCODE/practica/datos/"
outdir <- "~/Documentos/Conesa_Lab/VSCODE/practica/graficos/"
modo <- "rq"
combi <- paste(seq, plataforma, sep="_")
combi2 <- paste(modo, combi, sep="_")

# Lectura de datos
df_class <- read.table(paste0(ruta,paste(modo, combi, "rescued_classification", sep ="_"), ".txt"), sep="\t", header=TRUE)

df_counts <- read.table(paste0(ruta, paste("counts", combi, sep="_"), ".tsv"), sep="\t", header=TRUE)

# Generar data frame
# Media de los datos y transformación pseudologarítimica
df <- data.frame(tr_id = df_counts[,"superPBID"], 
                          K100 = apply(df_counts[, 2:6], 1, mean),
                          B100 = apply(df_counts[, 3:7], 1, mean),
                          B20 =  apply(df_counts[, 8:12], 1, mean) ,
                          B80 =  apply(df_counts[, 13:17], 1, mean) 
                          )
# Calcular mezclas teóricas para comparar con las experimentales
df["B20_ex"] <- df[, "B100"] * 0.2 + df[, "K100"] * 0.8
df["B80_ex"] <- df[, "B100"] * 0.8 + df[, "K100"] * 0.2       

# Unir los dos data frames
library(dplyr)

df <- df %>% 
        left_join(
                df_class %>% select("isoform", "associated_transcript", 
                                           "structural_category"),
                by = c("tr_id" ="isoform" )
        
)
# Filtramos para quitar noveles y que no sean FSM
df_filtered <- df[which(
        df[, "associated_transcript"] != "novel" &
        df[, "structural_category"] == "full-splice_match"), ]

# Transformación logarítmica base 10 con pseudocont 0.01

df_filtered[,2:7] <- apply(df_filtered[,2:7], 2, function(x) log10(x + 0.01))

# modelo lineal B20 K80

lin1 <- lm( B20_ex ~ B20, data= df_filtered)
lin_res1 <- summary(lin)
r1 <- lin_res$r.squared
nombre1 <- paste(paste(combi2, "B20K80")) 

# modelo lineal B80 K20

lin1 <- lm( B20_ex ~ B20, data= df_filtered)
lin_res1 <- summary(lin)
r1 <- lin_res$r.squared
nombre2 <- paste(paste(combi2, "B80K20")) 

# hacer el plot con función para etiquetar exponentes

label_10_pow <- function(x) {
        parse(text = paste0("10^", x))
}

pl <- ggplot(df_filtered, aes(x=B20_ex, y=B20)) +
        #geom_point(aes(color=B20_ex)) +
        # scale_color_gradientn(colors=topo.colors(6)) +
        geom_pointdensity(size= 0.8) +
        #scale_color_viridis_c(name="Point Density")+
        scale_color_viridis_c(
                name = "Point Density",
                guide = guide_colorbar(
                        barheight = unit(0.8, "npc"),            
                        barwidth = unit(1.2, "lines"),         
                        title.position = "right",               
                        title.theme = element_text(angle = -90, hjust = 0.5)
                )
        ) +
        geom_smooth(method="lm", formula= y~ x, se=T, level= 0.95, color="red") +
        geom_text_npc(
                data = data.frame(npcx = 0.05, npcy = 0.95, label = paste0("r = ", round(r1, 3))),
                aes(npcx = npcx, npcy = npcy, label = label),
                inherit.aes = FALSE,
                size = 4.5
        ) +
        scale_x_continuous(labels = label_10_pow) +
        scale_y_continuous(labels = label_10_pow) +
        labs(x= "Expected expression", 
             y= "Observed expression", 
             title= nombre1) +
        theme_minimal() +
        theme(plot.title=element_text(hjust=0.5))

ggsave(filename=nombre1, path=outdir, create.dir=T, device="png")


