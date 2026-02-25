# Formal data analysis
library(tidyverse)
library(terra)

month.names<-list.files("./species/mig_ENM/monthly/", full.names = F) |> str_remove("_ENMm")
month.folder<-list.files("./species/mig_ENM/monthly/", full.names = T)

season.names<-list.files("./species/mig_ENM/seasonaly/", full.names = F) |> str_remove("_ENMs")
season.folder<-list.files("./species/mig_ENM/seasonaly/", full.names = T)

month.list <- list()
season.list <- list()

# Apendice S1 -------------------------------------------------------------

# Cuantos registros de cada especie hay por mes? 
library(sf)
shp.list <- list.files("./species/mig_shapes/Selected/", pattern="*.shp", full.names = T)

rcs.list <- list()
for(x in 94:length(shp.list)){
  shp <- 
    read_sf(shp.list[x]) 
  
  rcs.list[[x]] <-  
    shp |>
    tibble() |>
    dplyr::summarise(n = n(), .by = c(month, species)) |>
    mutate(spp = basename(shp.list[x]) |> str_remove(".shp"))
  
  print(paste(basename(shp.list[x]) |> str_remove(".shp"), x))
}

# df.rcs <- 
do.call(rbind, rcs.list) |> 
  mutate(month =  case_when(
    month == 1 ~ "Ene",
    month == 2 ~ "Feb",
    month == 3 ~ "Mar",
    month == 4 ~ "Abr",
    month == 5 ~ "May",
    month == 6 ~ "Jun",
    month == 7 ~ "Jul",
    month == 8 ~ "Ago",
    month == 9 ~ "Sep",
    month == 10 ~ "Oct",
    month == 11 ~ "Nov",
    month == 12 ~ "Dic")) |> 
  pivot_wider(names_from = month, values_from = n) |> 
  select(!spp) |> 
  write.table("./species/mig_ENM/recordsxspp.txt",sep="\t", dec = ".", row.names=F)


# Apendice S2 y S3 --------------------------------------------------------

for(x in 1:length(month.names)){
  month.list[[x]] <- read.table(list.files(month.folder[x], full.names = T, pattern="ENM_mods_"), sep = "\t", dec=".", header=T) |> 
    mutate(type="month")
  season.list[[x]] <- read.table(list.files(season.folder[x], full.names = T, pattern="ENM_mods_"), sep = "\t", dec=".", header=T) |> 
    mutate(type="season")
}

 # enm.df <- 
  do.call(rbind, month.list) |> 
    tibble() |> 
    rename(ovlp = "month", OR = "om_rate_train", AUC_ratio="env_bg_paucratio") |> 
    select(species, ovlp, OR, AUC_ratio, type) |> 
    mutate(ovlp =  case_when(
      ovlp == 1 ~ "Ene",
      ovlp == 2 ~ "Feb",
      ovlp == 3 ~ "Mar",
      ovlp == 4 ~ "Abr",
      ovlp == 5 ~ "May",
      ovlp == 6 ~ "Jun",
      ovlp == 7 ~ "Jul",
      ovlp == 8 ~ "Ago",
      ovlp == 9 ~ "Sep",
      ovlp == 10 ~ "Oct",
      ovlp == 11 ~ "Nov",
      ovlp == 12 ~ "Dic")) |> 
    pivot_wider(names_from = ovlp, values_from = c(OR, AUC_ratio), names_vary="slowest")  |> 
    write.table("./species/mig_ENM/ENM_tenm_res_M.txt",sep="\t", dec = ".", row.names=F)
  
  
  do.call(rbind, season.list) |>
    tibble() |> 
    rename(ovlp = "season", OR = "om_rate_train", AUC_ratio="env_bg_paucratio") |> 
    select(species, ovlp, OR, AUC_ratio, type) |> 
    mutate(ovlp =  case_when(
      ovlp == 1 ~ "Inv",
      ovlp == 2 ~ "Pri",
      ovlp == 3 ~ "Ver",
      ovlp == 4 ~ "Oto")) |> 
    # filter(OR<0.06)
    pivot_wider(names_from = ovlp, values_from = c(OR, AUC_ratio), names_vary="slowest") |> 
    write.table("./species/mig_ENM/ENM_tenm_res_S.txt",sep="\t", dec = ".", row.names=F)
    
# Apendice S4 y S5 -------------------------------------------------------------

  data.ovrlp.m <- read.table(
    "./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt",
    sep = "\t",
    dec = ".",
    header = T)
  
  data.ovrlp.s <- read.table(
    "./species/mig_ENM/overlaps_tables/data.ovrlp_S_df.txt",
    sep = "\t",
    dec = ".",
    header = T)

# data.ovrlp.s |> 
data.ovrlp.m |>
  tibble() |> 
  select(spp, overlap, obs.D, p.D) |> 
  rename("Especie" = spp, "Superposicion" = overlap,"D Schoenner" = obs.D, "p-value" = p.D) |> 
  # write.table("./species/mig_ENM/overlaps_tables/D_Broenniman_S.txt",sep="\t", dec = ".", row.names=F)
  write.table("./species/mig_ENM/overlaps_tables/D_Broenniman_M.txt", sep="\t", dec = ".", row.names=F)


# Apendice S6 -------------------------------------------------------------
# Rangos del PCA
pca.paths.m <- list.files("./species/mig_ENM/overlaps_tables/PCA_m/", full.names=T)
pca.paths.s <- list.files("./species/mig_ENM/overlaps_tables/PCA_s/", full.names=T)

pca.m <- list()
pca.s <- list()

for (x in 1:length(pca.paths.m)) {
  pca.m[[x]] <- 
    read_tsv(pca.paths.m[x]) |> 
    filter(str_starts(elip,"abex_"))
  pca.s[[x]] <- 
    read_tsv(pca.paths.s[x]) |> 
    filter(str_starts(elip,"abex_"))
  
  print(paste(basename(pca.paths.s[x] |> str_remove("_s_PCA.txt")), x, "ready"))
}

pca.df.m <- 
  do.call(rbind, pca.m)

pca.df.s <- 
  do.call(rbind, pca.s)

amp.m<- pca.df.m |> 
  summarise(min.a1=min(Axis1), 
            max.a1=max(Axis1), 
            min.a2=min(Axis2), 
            max.a2=max(Axis2), 
            .by=elip) |> 
  extract(elip, 
          into = c("spp", "month"),
          regex = "(.*)_([0-9]+$)") |> 
  mutate(spp = spp |> str_remove_all("abex_"),
         cod = case_when(
           month == 1 ~ "In", 
           month == 2 ~ "In", 
           month == 3 ~ "In", 
           month == 4 ~ "T1", 
           month == 5 ~ "T1", 
           month == 6 ~ "T1", 
           month == 7 ~ "Ve", 
           month == 8 ~ "Ve", 
           month == 9 ~ "Ve", 
           month == 10 ~ "T2", 
           month == 11 ~ "T2", 
           month == 12 ~ "T2")) |> 
  summarise(min.a1.tot=min(min.a1),
            max.a1.tot=max(max.a1),
            min.a2.tot=min(min.a2),
            max.a2.tot=max(max.a2),
            .by=c(spp, cod)) |>
  # Calculo de la amplitud para cada temporada por especie
  mutate(
    # Primero detectamos los extremos de cada fila
    min_extremo = pmin(min.a1.tot, min.a2.tot),
    max_extremo = pmax(max.a1.tot, max.a2.tot),
    # min_extremo = pmin(min.a1, min.a2),
    # max_extremo = pmax(max.a1, max.a2),
    # Calculamos la diferencia
    amplitud = max_extremo - min_extremo) |> 
  select(spp, cod, min_extremo, max_extremo, amplitud)

amp.m |> 
  select(spp, cod, amplitud) |> 
  pivot_wider(values_from = amplitud, names_from = cod) |> 
  write.table("./outputs/tablas/ampl_m_total.txt", sep="\t", dec = ".", row.names=F)

