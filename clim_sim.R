# environmental coincidence degree - ecode
# climatic similarities analysis - clisa

# Analisis de similaridad de nicho existente con MOP
# La idea es hacer MOPs cruzados entre las temporadas (o meses) de las especies
library(tidyverse)
library(terra)

spp.path <- list.files("./Phd/borrar_phd/spp/", full.names = T)[-2]
spp.names <-basename(spp.path) |> str_remove("_ENMs")

comb.mop <- data.frame(a = rep(c(1:4),each=3), b= c(2,3,4,1,3,4,1,2,4,1,2,3))

r.path <- list.files("D:/Mauro/Phd/MOP/mean_season_5k", full.names = T)

for(y in 2:length(spp.path)) {
	# x <- 1
	env_temp <- new.env()
	load(list.files(paste0(spp.path[y]),pattern = ".Rdata$", full.names = T),envir =  env_temp)
	
	for(x in 1:nrow(comb.mop)) {
		# Cargar puntos
		pts.a <- env_temp$abbg.list[[comb.mop[x, 1]]]$temporal_df
		pts.b <- env_temp$abbg.list[[comb.mop[x, 2]]]$temporal_df
		# vectorizar
		pts.a.vect <- vect(pts.a, geom = c("lon", "lat"), crs = "EPSG:4326")
		pts.b.vect <- vect(pts.b, geom = c("lon", "lat"), crs = "EPSG:4326")
		
		# Definición LCC recomendada para cubrir América
		lambert_america <- 
			"+proj=lcc +lat_1=-10 +lat_2=65 +lat_0=0 +lon_0=-100 +datum=WGS84 +units=m +no_defs"
		
		pts.a.vec2 <- project(pts.a.vect, lambert_america)  # CRS UTM adecuado para México central
		pts.b.vec2 <- project(pts.b.vect, lambert_america)
		
		# Crear buffer de 250,000 metros (250 km - 50 celdas de 5km)
		buffer_250km.a <- buffer(pts.a.vec2, width = 250000) |>
			aggregate() |>
			project('EPSG:4326')
		buffer_250km.b <- buffer(pts.b.vec2, width = 250000) |>
			aggregate() |>
			project('EPSG:4326')
		
		# Ahora cargar los raster
		r.a <- rast(list.files(r.path[which(basename(r.path) == comb.mop[x, 1])], full.names = T))
		r.b <- rast(list.files(r.path[which(basename(r.path) == comb.mop[x, 2])], full.names = T))
		
		# Cortar los raster a las áreas de calibracion
		a.r.crop <- crop(r.a, buffer_250km.a, mask = T)
		b.r.crop <- crop(r.b, buffer_250km.b, mask = T)
		
		# MOP
		future::plan("future::multisession", workers = 10)
		mop_seas <- smop::mop(
			M_calibra = a.r.crop,
			G_transfer = b.r.crop,
			comp_each = 10000,
			percent = 50,
			standardize_vars = TRUE,
			normalized = TRUE)
		future::plan("future::sequential")
		
		writeRaster(mop_seas,
								paste0("./Phd/borrar_phd/mop_cross/mop_cross_", spp.names[y], "_",
											 comb.mop[x, 1], "_", comb.mop[x, 2], ".tif"))
		
		print(paste0(spp.names[y], " ", comb.mop[x, 1], " in ", comb.mop[x, 2], " done"))
		gc()
	}
}


## Ahora a extraer la información. Los puntos deben ser los del mes a comparar, el segundo
spp.path

mop.path <- list.files("./Phd/borrar_phd/mop_cross/", full.names = T)
spp.names <- 	mop.path |> 
	basename() |> 
	str_remove("mop_cross_") |> 
	str_remove(".tif") |> 
	substr(1,9) |> 
	unique()

spp.ext.mop <- list()
# x <- 1
for (x in 1:length(spp.names)){
	mop.months <-
		mop.path |>
		basename() |>
		str_remove("mop_cross_") |>
		str_remove(".tif") |>
		keep(~ str_detect(.x, spp.names[x])) |>
		str_sub(-3)
	
	mop.path2 <- 
		mop.path |>
		keep(~ str_detect(.x, spp.names[x]))
	
	env_temp <- new.env()
	load(list.files(spp.path[x], pattern = ".Rdata$", full.names = T), envir =  env_temp)
	
	ext.mop <- list()
	# Abrir mop
	# y <- 1
	for (y in 1:length(mop.path2)){
		mop.r <- rast(mop.path2[y])
		
		# Abrir datos de la especie
		occs <- 
			env_temp$abex.list[[as.numeric(mop.months[y] |> str_sub(-1))]]$temporal_df |> 
			select(2,3) |> 
			mutate(spp=spp.names[x],
						 mop.cros=mop.months[y]) |> 
			relocate(spp) |> 
			vect(geom=c("lon", "lat"), crs="EPSG:4326")
		
		ext.mop[[y]] <- 
			terra::extract(mop.r, occs, ID=F, bind=T) |> 
			as.data.frame()
	}
	
	spp.ext.mop[[x]] <- do.call(rbind, ext.mop)
	
}

# Ahora voy a intentar hacer un remuestreo aleatorio (n=100) en el área de calibración (MOP).
# el resultado ira en una lista


spp.ext.mop.df <- do.call(rbind, spp.ext.mop)
spp.ext.mop.df$mop.cros |> unique()

spp.ext.mop.df |> 
	summarise(mop.mean=mean(MOP), mop.sd=sd(MOP), .by=c(mop.cros, spp)) |> 
	mutate(ID=rep(1:12, 5),
				 lower = mop.mean - 1.96 * mop.sd,
				 upper = mop.mean + 1.96 * mop.sd) |> 
	ggplot(aes(x = ID, y = mop.mean, color=spp)) +
	geom_line(lwd=0.8) +
	geom_point() +
	# geom_hline(yintercept = me, color=c("#753BBDFF", "#00A9E0FF", "#006271FF", "#954E4CFF", "#010101FF"), linetype="dashed") +
	# geom_hline(yintercept = med, color="red", linetype="dashed") +
	# geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linetype = "dashed") +
	scale_x_continuous(breaks = 1:12, labels = spp.ext.mop.df$mop.cros |> unique()) +
	scale_color_manual(values=c("#753BBDFF", "#00A9E0FF", "#006271FF", "#954E4CFF", "#010101FF")) +
	labs(x="Crossing MOP", y="MOP (mean)", color= "Species") +
	theme_minimal()

spp.ext.mop.df |> 
	# summarise(mop.mean=mean(MOP), mop.sd=sd(MOP), .by=c(mop.cros, spp)) |> 
	summarise(mop.sd=sd(MOP), .by=spp) |> 
	ggplot(aes(x=spp, y=sort(mop.sd))) +
	geom_bar(stat='identity')

spp.ext.mop.df |> 
	ggplot(aes(x=mop.cros, y=MOP, fill=spp)) +
	ggstream::geom_stream()

spp.ext.mop.df |> 
	summarise(n=n(), .by=c(spp, mop.cros))

ggstream::blockbusters |> 
	summarise(n=n(), .by=c(genre, year))


spp.ext.mop.df$MOP |> range()

