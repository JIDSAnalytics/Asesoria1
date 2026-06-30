################################################################################
## JIDS ANALYTICS — Consultoría en Análisis Estadístico de la Sismicidad
## PRE-INFORME DE ASESORÍA 1
## Encargo: Alaska-Aleutianas oriental versus sur de Chile
##          (sismicidad en márgenes de alta latitud)
##
##
## Fuente de datos: catálogo USGS (FDSN), archivos crudos Alaska.csv y Chile.csv.
## base procesada final: base_procesada.csv   (en el script es "datos")



## Carpetas de salida ----------------------------------------------------------

dir.create("salidas",        showWarnings = FALSE)
dir.create("salidas/figuras", showWarnings = FALSE)
dir.create("salidas/tablas",  showWarnings = FALSE)

################################################################################


#### 0. PAQUETES ---------------------------------------------------------------

library(tidyverse)        # dplyr, ggplot2, tidyr, readr, etc.
library(lubridate)        # manejo de fechas/horas UTC
library(sf)               # objetos espaciales y mapas georreferenciados
library(rnaturalearth)    # mapa base mundial de costas/países
library(scales)           # formato de ejes
library(janitor)          # nombres de columnas limpios (opcional)
library(patchwork)        # combinación de varios gráficos en una sola figura.
library(ineq)             # cálculo de medidas de desigualdad, como el coeficiente de Gini
library(moments)          # Asimetría y curtosis

## TEMA GRÁFICO COMÚN ----------------------------------------

tema_jids <- theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.caption  = element_text(size = 8, color = "grey40"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
col_zonas <- c("A — Alaska-Aleutianas oriental" = "#1f78b4",
               "B — Sur de Chile"                = "#e31a1c")
## -----------------------------------------------------------


#### 1. PARÁMETROS DEL ENCARGO (trazabilidad de la descarga) -------------------

## Período de referencia y umbral de magnitud.
PERIODO_INI <- as.Date("2000-01-01")
PERIODO_FIN <- as.Date("2025-12-31")
N_ANIOS     <- 26          # años completos del período
N_MESES     <- 312         # 26 * 12 meses
M_MIN       <- 5.0         # umbral de magnitud para comparación internacional

## Cajas de coordenadas asignadas por la contraparte (N, S, O, E) 
zonas_box <- tibble::tribble(
  ~zona,                          ~norte, ~sur,  ~oeste, ~este,
  "A — Alaska-Aleutianas oriental",  62,    51,   -170,   -130,
  "B — Sur de Chile",               -38,   -56,    -77,    -68
)

#### 2. LECTURA DE LAS BASES CRUDAS --------------------------------------------
# Bases originales cuentan ambas con 22 columnas 

ruta_alaska <- "Alaska.csv"   # nombre de archivo en la ruta
ruta_chile  <- "Chile.csv"

raw_alaska <- readr::read_csv(ruta_alaska, show_col_types = FALSE) |>
  mutate(zona = "A — Alaska-Aleutianas oriental")

raw_chile  <- readr::read_csv(ruta_chile,  show_col_types = FALSE) |>
  mutate(zona = "B — Sur de Chile")

## Unión de ambas bases
datos_raw <- bind_rows(raw_alaska, raw_chile)



#### 3. CONTROL DE CALIDAD — diagnóstico sobre la base CRUDA -------------------

diagnostico_calidad <- datos_raw |>
  group_by(zona) |>
  summarise(
    registros_crudos      = n(),
    ids_duplicados        = sum(duplicated(id)),
    faltantes_tiempo      = sum(is.na(time)),
    faltantes_magnitud    = sum(is.na(mag)),
    faltantes_profundidad = sum(is.na(depth)),
    faltantes_coordenadas = sum(is.na(latitude) | is.na(longitude)),
    eventos_no_tectonicos = sum(type != "earthquake", na.rm = TRUE),
    fuera_de_periodo      = sum(as.Date(time) < PERIODO_INI |
                                as.Date(time) > PERIODO_FIN, na.rm = TRUE),
    bajo_umbral_M         = sum(mag < M_MIN, na.rm = TRUE),
    estado_no_reviewed    = sum(status != "reviewed", na.rm = TRUE),
    .groups = "drop"
  )

diagnostico_calidad
readr::write_csv(diagnostico_calidad, "salidas/tablas/diagnostico_calidad.csv")


#### 4. DEPURACIÓN Y VARIABLES DERIVADAS → BASE PROCESADA ----------------------

## Decisiones: se conservan solo eventos tectónicos
## Se documenta el evento de tipo "landslide" detectado en la zona A, que se excluye.

datos <- datos_raw |>
  mutate(
    fecha_hora = ymd_hms(time, tz = "UTC"),
    fecha      = as_date(fecha_hora)
  ) |>
  filter(
    type == "earthquake",                        # excluir no tectónicos (landslide)
    fecha >= PERIODO_INI, fecha <= PERIODO_FIN,  # fecha entre los límites
    mag   >= M_MIN                               # magnitud sobre el umbral  
  ) |>
  # Variables derivadas
  mutate(
    anio = year(fecha_hora),
    mes  = floor_date(fecha_hora, "month"),
    prof_categoria = case_when(           # clasificación operativa 
      depth <  70             ~ "Superficial (0–70)",
      depth >= 70 & depth <= 300 ~ "Intermedio (70–300)",
      depth >  300            ~ "Profundo (>300)",
      TRUE                    ~ NA_character_
    ) |> factor(levels = c("Superficial (0–70)",
                            "Intermedio (70–300)",
                            "Profundo (>300)")),

  ) |>
  arrange(zona, fecha_hora)  # ordena según zona y luego cronológicamente

## Exportar base procesada (reproducibilidad) (primer checkpoint)
readr::write_csv(datos, "salidas/base_procesada1.csv")




#### 5. HOMOGENEIZACIÓN DE MAGNITUDES a magnitud de momento --------------------
## Conversión de Scordilis (2006). Se crea la columna `Mw` con la magnitud
## homogeneizada (la variable original `mag` NO se modifica).

datos <- datos |>
  mutate(
    Mw = case_when(
      magType == "mb" ~ 0.85 * mag + 1.03,                         # mb  -> Mw
      magType == "ms" & mag >= 3.0 & mag <= 6.1 ~ 0.67*mag + 2.07, # Ms -> Mw
      magType == "ms" & mag >= 6.2 & mag <= 8.2 ~ 0.99*mag + 0.08, # Ms -> Mw
      magType %in% c("mw","mww","mwb","mwc","mwr") ~ mag,          # ya es Mw
      magType == "ml" ~ mag,                                       # ml conservada
      TRUE ~ mag                                                   # resto
    ),
    fuerte_60 = Mw >= 6.0,               # indicadores de eventos fuertes
    fuerte_65 = Mw >= 6.5,
    fuerte_70 = Mw >= 7.0
  )

## Exportar base procesada y con homogenización (segundo checkpoint) 
readr::write_csv(datos, "salidas/base_procesada2.csv")


## Boxplot de magnitud por tipo (justifica la homogeneización) 
fig_box_magtype <- ggplot(datos, aes(x = magType, y = mag)) +
  geom_boxplot(fill = "#1f78b4", color = "#16293d", alpha = .7,
               outlier.size = .8) +
  labs(title = "Distribución de la magnitud por tipo (magType)",
       subtitle = "Para analizar previo a homogeneización a Mw",
       x = "Tipo de magnitud (magType)", y = "Magnitud original (mag)",
       caption = "Fuente: catálogo USGS.") +
  tema_jids
fig_box_magtype
ggsave("salidas/figuras/fig_box_magtype.png", fig_box_magtype,
       width = 10, height = 5, dpi = 300)




#### * COMPONENTES GEOGRÁFICAS ----------------------------

## Mapa base mundial (costas/países) en proyección geográfica WGS84.
mundo <- rnaturalearth::ne_countries(scale = "medium",
                                     returnclass = "sf") # objeto tipo simple feature   

## Eventos como objeto espacial, separados por zona. Agrega "geometry" 
datos_sf <- st_as_sf(datos, coords = c("longitude", "latitude"),  
                     crs = 4326, remove = FALSE)   
# crs (coordinate reference system) = 4326. Proyección geográfica WGS84. 
# Sistema geodésico (estandar GPS)

alaska_sf <- filter(datos_sf, zona == "A — Alaska-Aleutianas oriental")
chile_sf  <- filter(datos_sf, zona == "B — Sur de Chile")


## Convertir las cajas de coordenadas a polígonos sf 
caja_a_poligono <- function(norte, sur, oeste, este) {
  st_polygon(list(matrix(c(
    oeste, sur,  este, sur,  este, norte,  oeste, norte,  oeste, sur
  ), ncol = 2, byrow = TRUE)))
}

rect_sf <- zonas_box |>
  rowwise() |>
  mutate(geometry = st_sfc(caja_a_poligono(norte, sur, oeste, este),
                           crs = 4326)) |>
  ungroup() |>
  st_as_sf()


# Límites de placas tectónicas (archivo externo PB2002) 
# Requiere PB2002_boundaries.shp (+ .shx, .dbf, .prj) en el directorio de trabajo.
placas <- st_read("PB2002_boundaries.shp")
st_crs(placas) <- 4326  # Mismo sistema de coordenadas de referencia

# Chile: Sudamericana (SA), Nazca (NZ), Antártica (AN).
# Alaska: Pacífico (PA), Norteamérica (NA).
placas_chile <- placas[placas$PlateA %in% c("SA","NZ","AN") |
                         placas$PlateB %in% c("SA","NZ","AN"), ]
placas_alaska <- placas[placas$PlateA %in% c("PA","NA") |
                          placas$PlateB %in% c("PA","NA"), ]


####-----------------------------------------------------------

#### * AGREGAMOS VARIABLE DERIVADA: Distancia a la fosa ------------------------

## Filtrar límites de subducción desde la capa de placas (PB2002)
fosa <- placas |>
  filter(Type == "subduction")

fosa_chile <- fosa |>
  filter(PlateA %in% c("SA","NZ","AN") |
           PlateB %in% c("SA","NZ","AN"))

fosa_alaska <- fosa |>
  filter(PlateA %in% c("PA","NA") |
           PlateB %in% c("PA","NA"))

## Proyectar a sistemas métricos (distancias en metros)
chile_m       <- st_transform(chile_sf,    32719)
alaska_m      <- st_transform(alaska_sf,   3338)
fosa_chile_m  <- st_transform(fosa_chile,  32719)
fosa_alaska_m <- st_transform(fosa_alaska, 3338)

## Distancia mínima de cada sismo a la fosa (km), por zona
alaska_sf$dist_fosa_km <- apply(st_distance(alaska_m, fosa_alaska_m), 1, min) / 1000
chile_sf$dist_fosa_km  <- apply(st_distance(chile_m,  fosa_chile_m),  1, min) / 1000

## Incorporar la columna a `datos` mediante un join por `id`.
## (alaska_sf y chile_sf provienen de datos_sf, que conserva `id`.)
dist_fosa_lookup <- bind_rows(
  st_drop_geometry(alaska_sf) |> select(id, dist_fosa_km),
  st_drop_geometry(chile_sf)  |> select(id, dist_fosa_km)
)

datos <- datos |>
  left_join(dist_fosa_lookup, by = "id")

## Base procesada ahora con la distancia incluida (base procesada final)
readr::write_csv(datos, "salidas/base_procesada.csv")




#### * Longitud de fosa por zona (km) para normalización por exposición --------
# Límites de las ventanas (los mismos de la descarga USGS):
#   Alaska: lat 51 a 62 ; Chile: lat -56 a -38

# Caja (bbox) de cada zona en coordenadas geográficas (EPSG:4326)
bbox_alaska <- st_as_sfc(st_bbox(c(xmin = -170, xmax = -130,
                                   ymin = 51,  ymax = 62),
                                 crs = 4326))
bbox_chile  <- st_as_sfc(st_bbox(c(xmin = -77,  xmax = -68,
                                   ymin = -56, ymax = -38),
                                 crs = 4326))

# Recortar la fosa de cada zona a su ventana y proyectar a métrico
fosa_alaska_clip <- st_intersection(st_transform(fosa_alaska, 4326), bbox_alaska) |>
  st_transform(3338)
fosa_chile_clip  <- st_intersection(st_transform(fosa_chile, 4326), bbox_chile) |>
  st_transform(32719)

# Longitud total de fosa (km) por zona
long_fosa_alaska_km <- as.numeric(sum(st_length(fosa_alaska_clip))) / 1000
long_fosa_chile_km  <- as.numeric(sum(st_length(fosa_chile_clip)))  / 1000



#### 6. TABLAS DE ESTADÍSTICOS DESCRIPTIVOS ------------------------------------

## Frecuencias de categorías de magnitud (Con datos originales) -----

# table(raw_alaska$magType)
# table(raw_chile$magType)

tabla_magtype <- datos |>
  count(zona, magType, name = "frecuencia") |>
  tidyr::pivot_wider(names_from = zona,
                     values_from = frecuencia,
                     values_fill = 0) |>
  arrange(magType)

print(tabla_magtype)
write_csv(tabla_magtype, "salidas/tablas/tabla_magtype.csv")


## Tabla descriptiva magType ----
# ¿sesga el tipo de magnitud? Esto sigue siendo con mag

tabla_magtype_resumen <- datos |>
  group_by(zona, magType) |>
  summarise(n = n(), mag_media = mean(mag), mag_min = min(mag),
            mag_max = max(mag), .groups = "drop")
print(tabla_magtype_resumen)
write_csv(tabla_magtype_resumen, "salidas/tablas/tabla_magtype_resumen.csv")



## Conteos y tasas de ocurrencia (Ya con datos nuevos) -----
# (tasa anual = n/26; tasa mensual = n/312). 

tabla_conteo_tasas <- datos |>
  group_by(zona) |>
  summarise(
    n_eventos    = n(),
    tasa_anual   = n() / N_ANIOS,
    tasa_mensual = n() / N_MESES,
    .groups = "drop"
  )
print(tabla_conteo_tasas)
write_csv(tabla_conteo_tasas, "salidas/tablas/tabla_conteos_tasas.csv")



## Conteos y tasas de ocurrencia (con normalización por longitud de fosa) ----

long_fosa <- tabla_conteo_tasas |>
  distinct(zona) |>
  mutate(long_fosa_km = if_else(grepl("Alaska", zona),
                                long_fosa_alaska_km,
                                long_fosa_chile_km))
print(long_fosa)


tabla_normalizada <- tabla_conteo_tasas |>
  left_join(long_fosa, by = "zona") |>
  mutate(
    # Tasa por exposición: eventos por año por 100 km de fosa
    tasa_fosa_100km = tasa_anual / long_fosa_km * 100
  )
print(tabla_normalizada)
write_csv(tabla_normalizada, "salidas/tablas/tabla_tasas_normalizadas.csv")

## Razones entre zonas (A respecto de B) ----
razon_bruta <- with(tabla_normalizada,
                    tasa_anual[grepl("^A", zona)] / tasa_anual[grepl("^B", zona)])

razon_fosa  <- with(tabla_normalizada,
                    tasa_fosa_100km[grepl("^A", zona)] / tasa_fosa_100km[grepl("^B", zona)])

cat("Razón de tasas A/B — bruta:", round(razon_bruta, 2),
    "| por longitud de fosa:", round(razon_fosa, 2), "\n")

## Distribución de magnitud -----

tabla_dist_Mw <- datos |>
  group_by(zona) |>
  summarise(
    minimo   = min(Mw),
    media    = mean(Mw),
    mediana  = median(Mw),
    maximo   = max(Mw),
    desv_est = sd(Mw),
    q25 = quantile(Mw, .25),
    q75 = quantile(Mw, .75),
    q90 = quantile(Mw, .90),
    q95 = quantile(Mw, .95),
    asimetria  = skewness(Mw),
    curtosis   = kurtosis(Mw),
    .groups = "drop"
  )
print(tabla_dist_Mw)
write_csv(tabla_dist_Mw, "salidas/tablas/tabla_dist_Mw.csv")


fig_dist_mag <- ggplot(datos, aes(Mw, fill = zona, color = zona)) +
  geom_density(alpha = .25, linewidth = .7) +
  scale_fill_manual(values = col_zonas) +
  scale_color_manual(values = col_zonas) +
  labs(title = "Distribución de magnitud por zona",
       x = "Magnitud (Mw)", y = "Densidad", fill = NULL, color = NULL,
       caption = "Fuente: catálogo USGS. Eventos M ≥ 5,0, 2000–2025.") +
  tema_jids
fig_dist_mag
ggsave("salidas/figuras/fig_dist_mag.png", fig_dist_mag,
       width = 9, height = 5, dpi = 300)

## Distribución de profundidad y composición por categoría -----

tabla_dist_depth <- datos |>
  group_by(zona) |>
  summarise(
    min_prof = min(depth), media_prof = mean(depth),
    mediana_prof = median(depth), max_prof = max(depth),
    asimetria  = skewness(depth), curtosis = kurtosis(depth),
    .groups = "drop"
  )

## Distribución de profundidad (densidad por zona) ────
fig_dist_depth <- ggplot(datos, aes(depth, fill = zona, color = zona)) +
  geom_density(alpha = .25, linewidth = .7) +
  scale_fill_manual(values = col_zonas) +
  scale_color_manual(values = col_zonas) +
  labs(title = "Distribución de profundidad por zona",
       x = "Profundidad (Km)", y = "Densidad", fill = NULL, color = NULL,
       caption = "Fuente: catálogo USGS. Eventos M ≥ 5,0, 2000–2025.") +
  tema_jids
fig_dist_depth
ggsave("salidas/figuras/fig_dist_depth.png", fig_dist_depth,
       width = 9, height = 5, dpi = 300)


tabla_depth_cat <- datos |>
  count(zona, prof_categoria) |>
  group_by(zona) |>
  mutate(proporcion = n / sum(n)) |>
  ungroup()

print(tabla_dist_depth); print(tabla_depth_cat)
write_csv(tabla_dist_depth,      "salidas/tablas/tabla_dist_depth.csv")
write_csv(tabla_depth_cat, "salidas/tablas/tabla_depth_cat.csv")


## Distribución y tabla de distancias a la fosa ----

tabla_fosa <- datos |>
  group_by(zona) |>
  summarise(n = n(),
            distancia_media   = mean(dist_fosa_km),
            distancia_mediana = median(dist_fosa_km),
            q25 = quantile(dist_fosa_km, .25),
            q75 = quantile(dist_fosa_km, .75),
            max_dist = max(dist_fosa_km),
            .groups = "drop")
print(tabla_fosa)
write_csv(tabla_fosa, "salidas/tablas/tabla_distancia_fosa.csv")



fig_dist_fosa <- ggplot(datos, aes(dist_fosa_km, fill = zona, color = zona)) +
  geom_density(alpha = .25, linewidth = .8) +
  scale_fill_manual(values = col_zonas) +
  scale_color_manual(values = col_zonas) +
  labs(title = "Distribución de la distancia a la zona de subducción",
       x = "Distancia a la fosa (km)", y = "Densidad", fill = NULL, color = NULL,
       caption = "Fuente: catálogo USGS. Eventos M ≥ 5,0, 2000–2025.") +
  tema_jids
fig_dist_fosa
ggsave("salidas/figuras/fig_distancia_fosa.png", fig_dist_fosa,
       width = 9, height = 5, dpi = 300)


## Eventos fuertes y extremos (n y proporción) -----

tabla_eventos_fuertes <- datos |>
  group_by(zona) |>
  summarise(
    total = n(),
    n_M60 = sum(fuerte_60), p_M60 = mean(fuerte_60),
    n_M65 = sum(fuerte_65), p_M65 = mean(fuerte_65),
    n_M70 = sum(fuerte_70), p_M70 = mean(fuerte_70),
    .groups = "drop"
  )
print(tabla_eventos_fuertes) 
write_csv(tabla_eventos_fuertes, "salidas/tablas/tabla_eventos_fuertes.csv")


## Temporalidad ----

conteos_mensuales <- datos |>
  count(zona, mes, name = "n") |>
  group_by(zona) |>
  complete(mes = seq(floor_date(PERIODO_INI, "month"),
                     floor_date(PERIODO_FIN, "month"), by = "month"),
           fill = list(n = 0)) |>
  arrange(zona,mes) |>
  ungroup()


## ACF descriptiva de conteos mensuales ----

acf_df <- conteos_mensuales |>
  group_by(zona) |>
  group_modify(~ {
    n_obs <- sum(!is.na(.x$n)) 
    ac <- acf(.x$n, plot = FALSE, na.action = na.pass, lag.max = 36)
    
    tibble(
      lag = as.numeric(ac$lag), 
      acf = as.numeric(ac$acf),
      ci = 1.96 / sqrt(n_obs) # Cada zona tiene su propio valor aquí
    )
  }) |>
  ungroup()

write_csv(acf_df, "salidas/tablas/acf_conteos_mensuales.csv")

# Gráfico con bandas segmentadas por zona 
fig_acf <- ggplot(acf_df, aes(x = lag, y = acf, color = zona)) +
  geom_hline(yintercept = 0, color = "grey40") +
  
  # SOLUCIÓN: geom_segment obliga a ggplot a respetar el 'ci' de cada panel
  geom_segment(aes(x = min(lag), xend = max(lag), y = ci, yend = ci), 
               linetype = "dashed", color = "firebrick", alpha = 0.6) +
  geom_segment(aes(x = min(lag), xend = max(lag), y = -ci, yend = -ci), 
               linetype = "dashed", color = "firebrick", alpha = 0.6) +
  
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.7) +
  facet_wrap(~ zona) +
  scale_color_manual(values = col_zonas, guide = "none") + 
  labs(
    title = "ACF descriptiva de los conteos mensuales",
    subtitle = "Líneas discontinuas rojas indican el límite de significancia (95%) propio de cada zona",
    x = "Rezago mensual",
    y = "Autocorrelación",
    caption = "Fuente: catálogo USGS. Período 2000–2025."
  ) +
  tema_jids

fig_acf
ggsave("salidas/figuras/fig_acfs.png", fig_acf,
       width = 9, height = 5, dpi = 300)


## Tiempos entre eventos ----

tiempos_entre <- datos |>
  group_by(zona) |>
  mutate(dt_dias = as.numeric(difftime(fecha_hora, lag(fecha_hora),
                                       units = "days"))) |>
  summarise(
    dt_medio   = mean(dt_dias, na.rm = TRUE),
    dt_mediana = median(dt_dias, na.rm = TRUE),
    .groups = "drop"
  )

tiempos_entre_eventos <- datos |>
  group_by(zona) |>
  arrange(fecha_hora, .by_group = TRUE) |>
  mutate(dt_dias = as.numeric(difftime(fecha_hora, lag(fecha_hora), units = "days"))) |>
  ungroup() |>
  filter(is.finite(dt_dias), dt_dias >= 0)

resumen_dt <- tiempos_entre_eventos |>
  group_by(zona) |>
  summarise(
    dt_media = mean(dt_dias, na.rm = TRUE),
    dt_mediana = median(dt_dias, na.rm = TRUE),
    dt_max = max(dt_dias, na.rm = TRUE),
    dt_q25 = quantile(dt_dias, 0.25, na.rm = TRUE),
    dt_q75 = quantile(dt_dias, 0.75, na.rm = TRUE),
    dt_q90 = quantile(dt_dias, 0.90, na.rm = TRUE),
    dt_q95 = quantile(dt_dias, 0.95, na.rm = TRUE),
    asimetria  = skewness(dt_dias), 
    curtosis = kurtosis(dt_dias),
    .groups = "drop"
  )
resumen_dt

write_csv(resumen_dt, "salidas/tablas/resumen_tiempos_entre_eventos.csv")


## Se estima la tasa de cada zona (1/media) para definir la exponencial teórica.
tasas_exp <- tiempos_entre_eventos |>
  group_by(zona) |>
  summarise(rate = 1 / mean(dt_dias, na.rm = TRUE), .groups = "drop")

## QQ-plot frente a la exponencial, una recta de referencia por zona.
fig_qq_dt <- ggplot(tiempos_entre_eventos, aes(sample = dt_dias, color = zona)) +
  stat_qq(distribution = qexp, alpha = .5, size = 1) +
  stat_qq_line(distribution = qexp, color = "grey20", linewidth = .7) +
  scale_color_manual(values = col_zonas) +
  facet_wrap(~ zona, scales = "free") +
  labs(title = "QQ-plot del tiempo entre eventos frente a la Exponencial",
       x = "Cuantiles teóricos (Exponencial)",
       y = "Cuantiles muestrales (días entre eventos)",
       color = NULL,
       caption = "Fuente: catálogo USGS. Recta = ajuste exponencial teórico.") +
  tema_jids + theme(legend.position = "none")
fig_qq_dt
ggsave("salidas/figuras/fig_qq_tiempos.png", fig_qq_dt,
       width = 10, height = 5, dpi = 300)



## Magnitud máxima anual -----

# Magnitud máxima por año y zona (producto central del encargo).
mag_max_anual <- datos |>
  group_by(zona, anio) |>
  summarise(mag_max = max(Mw), .groups = "drop")
print(mag_max_anual) 
write_csv(mag_max_anual,"salidas/tablas/mag_max_anual.csv")

# Tabla resumen
tabla_mag_max <- mag_max_anual |> 
  group_by(zona) |>
  summarise(
    mag_max_anual_media   = mean(mag_max),
    mag_max_anual_mediana = median(mag_max),
    mag_max_anual_sd      = sd(mag_max),          # variabilidad año a año
    mag_max_anual_min     = min(mag_max),         # el "peor" año (más tranquilo)
    mag_max_periodo       = max(mag_max),         # máximo absoluto del período
    anio_mag_max          = anio[which.max(mag_max)],  # año del sismo mayor
    n_anios               = n(),                  # años con registro
    .groups = "drop"
  )

print(tabla_mag_max)
write_csv(tabla_mag_max,    "salidas/tablas/tabla_mag_max.csv")




#### 7. ANÁLISIS BIVARIADO -----

## Matriz de correlación ---- 

vars_cont <- c("latitude", "longitude", "depth", "Mw", "dist_fosa_km")

# Función que devuelve la matriz de correlación en formato largo (tidy),
# etiquetando método y zona, para poder apilar ambas zonas y exportar.
matriz_cor_tidy <- function(df, metodo) {
  m <- cor(df[, vars_cont], method = metodo, use = "complete.obs")
  as.data.frame(as.table(m)) |>
    rename(var1 = Var1, var2 = Var2, correlacion = Freq) |>
    mutate(metodo = metodo)
}

tabla_correlaciones <- datos |>
  group_split(zona) |>
  lapply(function(df) {
    z <- unique(df$zona)
    bind_rows(matriz_cor_tidy(df, "pearson"),
              matriz_cor_tidy(df, "spearman")) |>
      mutate(zona = z)
  }) |>
  bind_rows() |>
  relocate(zona, metodo)

print(tabla_correlaciones)
write_csv(tabla_correlaciones, "salidas/tablas/tabla_correlaciones.csv")

## Figura: corrplot (Spearman)  
# Se usa Spearman por la no normalidad.
fig_corrplot <- tabla_correlaciones |>
  filter(metodo == "spearman") |>
  ggplot(aes(var1, var2, fill = correlacion)) +
  geom_tile(color = "white", linewidth = .6) +
  geom_text(aes(label = sprintf("%.2f", correlacion)), size = 3.2) +
  facet_wrap(~ zona) +
  scale_fill_gradient2(name = "Spearman ρ",
                       low = "#e31a1c", mid = "white", high = "#1f78b4",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Matriz de correlación de Spearman entre variables continuas",
       subtitle = "Asociación monótona por zona (lat, lon, profundidad, magnitud, distancia a fosa)",
       x = NULL, y = NULL,
       caption = "Fuente: catálogo USGS. Medida descriptiva de asociación, no causal.") +
  tema_jids +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
fig_corrplot
ggsave("salidas/figuras/fig_corrplot.png", fig_corrplot,
       width = 11, height = 5.5, dpi = 300)



## Profundidad y Magnitud categóricas ----
mag_depth_cat <- datos %>%
  filter(!is.na(depth), !is.na(Mw)) %>%
  mutate(
    cat_profundidad = case_when(depth <= 70 ~ "Superficial (0-70 km)", depth > 70 & depth <= 300 ~ "Intermedio (70-300 km)", TRUE ~ "Profundo (>300 km)"),
    cat_profundidad = factor(cat_profundidad, levels = c("Superficial (0-70 km)", "Intermedio (70-300 km)", "Profundo (>300 km)")),
    rango_magnitud = case_when(Mw >= 5.0 & Mw < 5.5 ~ "5.0 ≤ Mw < 5.5 (Moderado)", Mw >= 5.5 & Mw < 6.5 ~ "5.5 ≤ Mw < 6.5 (Fuerte)", TRUE ~ "Mw ≥ 6.5 (Extremo)"),
    rango_magnitud = factor(rango_magnitud, levels = c("5.0 ≤ Mw < 5.5 (Moderado)", "5.5 ≤ Mw < 6.5 (Fuerte)", "Mw ≥ 6.5 (Extremo)"))
  )

ggplot(mag_depth_cat, aes(x = cat_profundidad, fill = rango_magnitud)) +
  geom_bar(position = "fill", width = 0.55) + 
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~zona) +
  # CAMBIO AQUÍ: Usamos la paleta "Blues" para tonos celestes/azules
  scale_fill_brewer(palette = "Blues", name = "Severidad (Magnitud)") +
  labs(
    title = "Composición Estructural de Magnitudes por Estrato de Profundidad",
    subtitle = "Análisis del porcentaje acumulado para verificar la severidad según el nivel tectónico",
    x = "Categoría de Profundidad Operativa", y = "Proporción Relativa (100%)"
  ) + 
  theme_minimal() + 
  theme(
    plot.title = element_text(face="bold", size=13), 
    plot.subtitle = element_text(size=10, face="italic"), 
    legend.position = "bottom", 
    strip.text = element_text(face="bold")
  )

cramer_v <- function(df) {
  # droplevels() quita los niveles de factor que no tienen eventos en este df
  prof <- droplevels(df$cat_profundidad)
  mag  <- droplevels(df$rango_magnitud)
  tabla <- table(prof, mag)
  chi2  <- suppressWarnings(chisq.test(tabla)$statistic)
  n     <- sum(tabla)
  k     <- min(nrow(tabla), ncol(tabla)) - 1   # gl mínimo
  as.numeric(sqrt(chi2 / (n * k)))
}

tabla_cramer <- mag_depth_cat |>
  group_split(zona) |>
  lapply(function(df) {
    tibble(zona = unique(df$zona), V_Cramer = cramer_v(df))
  }) |>
  bind_rows()

print(tabla_cramer)
write_csv(tabla_cramer, "salidas/tablas/tabla_cramer.csv")


## Profundidad magnitud numéricas -----

## Hexbin Profundidad ↔ Magnitud 
#  "¿la magnitud cambia con la profundidad del hipocentro?".

fig_depth_mag <- ggplot(datos, aes(depth, Mw)) +
  geom_hex(bins = 30) +
  geom_vline(xintercept = 70, linetype = "dashed", color = "grey40") +
  facet_wrap(~ zona, scales = "free_x") +
  scale_fill_viridis_c(name = "N° eventos", option = "C") +
  labs(title = "Relación profundidad–magnitud por zona",
       subtitle = "Densidad hexagonal. Línea discontinua: límite superficial/intermedio (70 km)",
       x = "Profundidad (km)", y = "Magnitud (Mw)",
       caption = "Fuente: catálogo USGS. Exploratorio.") +
  tema_jids
fig_depth_mag

ggsave("salidas/figuras/fig_depth_mag.png", fig_depth_mag,
       width = 11, height = 5, dpi = 300)


## Correlación depth–mag por categoría de profundidad
tabla_depthmag_porcat <- datos |>
  group_by(zona, prof_categoria) |>
  summarise(
    n          = n(),
    rho_spearman = if (n() >= 3) cor(depth, mag, method = "spearman") else NA_real_,
    mag_media  = mean(mag),
    p_M65      = mean(fuerte_65),     # proporción de eventos fuertes en la categoría
    .groups = "drop"
  )
print(tabla_depthmag_porcat)
write_csv(tabla_depthmag_porcat, "salidas/tablas/tabla_depthmag_porcat.csv")



## Latitud/Longitud vs magnitud ----

## Perfil de magnitud media por banda LATITUDINAL (1°)
perfil_lat_mag <- datos |>
  mutate(banda_lat = floor(latitude)) |>
  group_by(zona, banda_lat) |>
  summarise(mag_media = mean(Mw), mag_max = max(Mw), n = n(),
            .groups = "drop")

fig_lat_mag <- ggplot(perfil_lat_mag, aes(banda_lat, mag_media, color = zona)) +
  geom_line(linewidth = .8) +
  geom_point(aes(size = n), alpha = .7) +
  facet_wrap(~ zona, scales = "free_x") +
  scale_color_manual(values = col_zonas, guide = "none") +
  scale_size_continuous(name = "N° eventos") +
  labs(title = "Perfil latitudinal de la magnitud media",
       subtitle = "Magnitud media por banda de 1° de latitud",
       x = "Latitud (°)", y = "Magnitud media (M)",
       caption = "Fuente: catálogo USGS. Exploratorio.") +
  tema_jids
fig_lat_mag
ggsave("salidas/figuras/fig_lat_mag.png", fig_lat_mag,
       width = 11, height = 5, dpi = 300)

## Perfil de magnitud media por banda LONGITUDINAL (1°)
perfil_lon_mag <- datos |>
  mutate(banda_lon = floor(longitude)) |>
  group_by(zona, banda_lon) |>
  summarise(mag_media = mean(Mw), mag_max = max(Mw), n = n(),
            .groups = "drop")

fig_lon_mag <- ggplot(perfil_lon_mag, aes(banda_lon, mag_media, color = zona)) +
  geom_line(linewidth = .8) +
  geom_point(aes(size = n), alpha = .7) +
  facet_wrap(~ zona, scales = "free_x") +
  scale_color_manual(values = col_zonas, guide = "none") +
  scale_size_continuous(name = "N° eventos") +
  labs(title = "Perfil longitudinal de la magnitud media",
       subtitle = "Magnitud media por banda de 1° de longitud",
       x = "Longitud (°)", y = "Magnitud media (M)",
       caption = "Fuente: catálogo USGS. Exploratorio.") +
  tema_jids
fig_lon_mag
ggsave("salidas/figuras/fig_lon_mag.png", fig_lon_mag,
       width = 11, height = 5, dpi = 300)



## Magnitud vs. distancia a la fosa -----
fig_mag_dist <- ggplot(datos, aes(dist_fosa_km, mag, color = zona)) +
  geom_point(alpha = .15) +
  # geom_smooth() eliminado para mostrar solo la nube de puntos pura
  scale_color_manual(values = col_zonas) +
  labs(title = "Magnitud en función de la distancia a la fosa",
       x = "Distancia a la fosa (km)", y = "Magnitud", color = NULL,
       caption = "Fuente: catálogo USGS. Gráfico de dispersión puro.") + # Breve ajuste al caption
  tema_jids
fig_mag_dist
ggsave("salidas/figuras/fig_magnitud_distancia_fosa.png", fig_mag_dist,
       width = 9, height = 5, dpi = 300)

## Profundidad vs. distancia a la fosa -----
fig_depth_dist <- ggplot(datos, aes(dist_fosa_km, depth, color = zona)) +
  geom_point(alpha = .15) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = col_zonas) +
  labs(title = "Profundidad en función de la distancia a la fosa",
       x = "Distancia a la fosa (km)", y = "Profundidad (km)", color = NULL,
       caption = "Fuente: catálogo USGS. Suavizado loess descriptivo.") +
  tema_jids
fig_depth_dist
ggsave("salidas/figuras/fig_profundidad_distancia_fosa.png", fig_depth_dist,
       width = 9, height = 5, dpi = 300)


## Distancia a la fosa por categoría de magnitud -----
datos_umbrales <- datos |>
  select(zona, dist_fosa_km, fuerte_60, fuerte_65, fuerte_70) |>
  tidyr::pivot_longer(
    cols      = c(fuerte_60, fuerte_65, fuerte_70),
    names_to  = "umbral",
    values_to = "supera"
  ) |>
  filter(supera) |>                      # conserva solo los que superan cada umbral
  mutate(umbral = recode(umbral,
                         fuerte_60 = "Mw ≥ 6,0",
                         fuerte_65 = "Mw ≥ 6,5",
                         fuerte_70 = "Mw ≥ 7,0") |>
           factor(levels = c("Mw ≥ 6,0", "Mw ≥ 6,5", "Mw ≥ 7,0")))

fig_fuertes_fosa <- ggplot(datos_umbrales, aes(umbral, dist_fosa_km, fill = zona)) +
  geom_boxplot(alpha = .8) +
  scale_fill_manual(values = col_zonas) +
  labs(title = "Distancia a la fosa según umbral de magnitud",
       subtitle = "Categorías acumulativas: cada evento aparece en todos los umbrales que supera",
       x = NULL, y = "Distancia a la fosa (km)", fill = NULL,
       caption = "Fuente: catálogo USGS. Eventos M ≥ 5,0, 2000–2025.") +
  tema_jids
fig_fuertes_fosa
ggsave("salidas/figuras/fig_fuertes_distancia_fosa.png", fig_fuertes_fosa,
       width = 9, height = 5, dpi = 300)



## Distancia a la fosa por categoría de profundidad -----
fig_profcat_fosa <- ggplot(datos, aes(prof_categoria, dist_fosa_km, fill = zona)) +
  geom_boxplot(alpha = .8) +
  scale_fill_manual(values = col_zonas) +
  labs(title = "Distancia a la fosa por categoría de profundidad",
       x = NULL, y = "Distancia a la fosa (km)", fill = NULL) +
  tema_jids
fig_profcat_fosa
ggsave("salidas/figuras/fig_profcat_distancia_fosa.png", fig_profcat_fosa,
       width = 9, height = 5, dpi = 300)

## Tabla: distancia a la fosa de eventos fuertes -----
## NOTA: median/mean sobre subconjuntos vacíos devuelven NaN/NA. Si una zona no
## tiene eventos M ≥ 7,0, esas celdas saldrán NaN; es correcto, indica ausencia.
tabla_fuertes_fosa <- datos |>
  group_by(zona) |>
  summarise(mediana_total = median(dist_fosa_km),
            mediana_M65 = median(dist_fosa_km[Mw >= 6.5]),
            mediana_M70 = median(dist_fosa_km[Mw >= 7.0]),
            media_total = mean(dist_fosa_km),
            mediana_M65 = median(dist_fosa_km[Mw >= 6.5]),
            mediana_M70 = median(dist_fosa_km[Mw >= 7.0]),
            .groups = "drop")
print(tabla_fuertes_fosa)
write_csv(tabla_fuertes_fosa, "salidas/tablas/tabla_fuertes_distancia_fosa.csv")


## Correlación de Spearman entre distancia y profundidad por zona -----
datos |>
  group_by(zona) |>
  summarise(rho = cor(dist_fosa_km, depth, method = "spearman"), .groups = "drop")




#### 8. ANÁLISIS MULTIVARIADO -------------------------------------------------------


## Regresión lineal múltiple exploratoria ----
## Tres modelos por zona. Devuelve R² y R² ajustado de cada uno.
## NOTA: medida descriptiva del poder explicativo, no inferencia formal.

ajustar_r2 <- function(df, formula_txt) {
  m <- lm(as.formula(formula_txt), data = df)
  s <- summary(m)
  tibble(
    modelo      = formula_txt,
    R2          = s$r.squared,
    R2_ajustado = s$adj.r.squared
  )
}

modelos <- c(
  "Mw ~ latitude + longitude + depth",
  "Mw ~ dist_fosa_km",
  "depth ~ latitude + longitude + Mw",
  "depth ~ dist_fosa_km"
)

tabla_r2 <- datos |>
  group_split(zona) |>
  lapply(function(df) {
    z <- unique(df$zona)
    lapply(modelos, function(f) ajustar_r2(df, f)) |>
      bind_rows() |>
      mutate(zona = z)
  }) |>
  bind_rows() |>
  relocate(zona)

print(tabla_r2)
write_csv(tabla_r2, "salidas/tablas/tabla_r2_regresiones.csv")




## Latitud ↔ Profundidad, Magnitud  y  Longitud ↔ Profundidad, Magnitud ----
##   EL PERFIL DE SUBDUCCIÓN 
##         Al graficar la profundidad del hipocentro contra la coordenada
##         que CRUZA el margen, debe aparecer el plano de Wadati–Benioff:
##         los eventos se hacen más profundos al alejarse de la fosa.


## Profundidad, Magnitud vs LONGITUD (clave para Chile)
fig_lon_depth <- ggplot(datos, aes(longitude, depth, color = Mw)) +
  geom_point(alpha = .6, size = 1.4) +
  # geom_smooth() eliminado para mostrar solo los puntos
  facet_wrap(~ zona, scales = "free_x") +
  scale_y_reverse() +
  scale_color_viridis_c(name = "Magnitud (M)", option = "B") +
  labs(title = "Perfil de profundidad vs longitud (corte perpendicular al margen)",
       subtitle = "Estructura de subducción: la profundidad crece al alejarse de la fosa",
       x = "Longitud (°)", y = "Profundidad (km) — eje invertido",
       caption = "Fuente: catálogo USGS. Exploratorio.") +
  tema_jids
fig_lon_depth
ggsave("salidas/figuras/fig_lon_depth.png", fig_lon_depth,
       width = 11, height = 5.5, dpi = 300)

## Profundidad, Magnitud vs LATITUD (clave para el arco de Alaska)
fig_lat_depth <- ggplot(datos, aes(latitude, depth, color = Mw)) +
  geom_point(alpha = .6, size = 1.4) +
  # geom_smooth() eliminado para mostrar solo los puntos
  facet_wrap(~ zona, scales = "free_x") +
  scale_y_reverse() +
  scale_color_viridis_c(name = "Magnitud (M)", option = "B") +
  labs(title = "Perfil de profundidad vs latitud",
       subtitle = "Complemento del perfil de subducción (relevante en el arco oblicuo de Alaska)",
       x = "Latitud (°)", y = "Profundidad (km) — eje invertido",
       caption = "Fuente: catálogo USGS. Exploratorio.") +
  tema_jids
fig_lat_depth
ggsave("salidas/figuras/fig_lat_depth.png", fig_lat_depth,
       width = 11, height = 5.5, dpi = 300)




#### 9. FIGURAS PRELIMINARES DESCRIPTIVAS --------------------------------------

## Serie temporal de conteos mensuales -----

fig_conteos_mensuales <- ggplot(conteos_mensuales, aes(mes, n, color = zona)) +
  geom_line(alpha = .9, linewidth = .4) +
  scale_color_manual(values = col_zonas) +
  labs(title = "Conteo mensual de sismos por zona (M ≥ 5,0)",
       x = "Mes", y = "N° de eventos", color = NULL,
       caption = "Fuente: catálogo USGS. Período 2000–2025.") +
  tema_jids
fig_conteos_mensuales
ggsave("salidas/figuras/fig_conteos_mensuales.png", fig_conteos_mensuales,
       width = 9, height = 5, dpi = 300)

## Serie temporal de conteos anuales -----
conteos_anuales <- datos |> count(zona, anio)
fig_conteos_anuales <- ggplot(conteos_anuales, aes(anio, n, color = zona)) +
  geom_line(linewidth = .8) + geom_point(size = 1.6) +
  scale_color_manual(values = col_zonas) +
  scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  labs(title = "Conteo anual de sismos por zona (M ≥ 5,0)",
       x = "Año", y = "N° de eventos", color = NULL,
       caption = "Fuente: catálogo USGS. Período 2000–2025.") +
  tema_jids
fig_conteos_anuales
ggsave("salidas/figuras/fig_conteos_anuales.png", fig_conteos_anuales,
       width = 9, height = 5, dpi = 300)

## ST Magnitud máxima anual -----
fig_mag_max_anual <- ggplot(mag_max_anual, aes(anio, mag_max, color = zona)) +
  geom_line(linewidth = .8) + geom_point(size = 1.8) +
  scale_color_manual(values = col_zonas) +
  scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  labs(title = "Magnitud máxima anual por zona",
       x = "Año", y = "Magnitud máxima (Mw)", color = NULL,
       caption = "Fuente: catálogo USGS. Un valor por año y zona.") +
  tema_jids
fig_mag_max_anual
ggsave("salidas/figuras/fig_mag_max_anual.png", fig_mag_max_anual,
       width = 9, height = 5, dpi = 300)



## Distribución de Tiempo entre eventos  ----
# Con curva exponencial visual con media empírica por zona
exp_grid <- resumen_dt |>
  group_by(zona) |>
  group_modify(~ {
    media <- .x$dt_media[1]
    xmax <- quantile(tiempos_entre_eventos$dt_dias[tiempos_entre_eventos$zona == .y$zona], 0.98, na.rm = TRUE)
    tibble(dt_dias = seq(0, xmax, length.out = 200), dens_exp = dexp(dt_dias, rate = 1 / media))
  }) |>
  ungroup()

fig_dt <- ggplot(tiempos_entre_eventos, aes(x = dt_dias, fill = zona, color = zona)) +
  geom_histogram(aes(y = after_stat(density)), bins = 45, alpha = 0.25, position = "identity") +
  geom_line(data = exp_grid, aes(x = dt_dias, y = dens_exp, color = zona), linewidth = 1.1) +
  facet_wrap(~ zona, scales = "free") +
  scale_fill_manual(values = col_zonas) +
  scale_color_manual(values = col_zonas) +
  coord_cartesian(xlim = c(0, quantile(tiempos_entre_eventos$dt_dias, 0.98, na.rm = TRUE))) +
  labs(
    title = "Distribución de tiempos entre eventos consecutivos",
    subtitle = "Histograma y curva exponencial visual con media empírica",
    x = "Días entre eventos",
    y = "Densidad",
    fill = "Zona",
    color = "Zona",
    caption = "Fuente: catálogo USGS. Eventos M ≥ 5,0, 2000–2025."
  ) +
  tema_jids
print(fig_dt)
ggsave("salidas/figuras/fig_dt.png", fig_dt,
       width = 9, height = 5, dpi = 300)

## Gráfico Curva de Lorenz: concentración temporal de la sismicidad ----
## Ajustes previos -----------------------------------------------------

# 1. Conteos mensuales por zona, rellenando con 0 los meses sin sismos 
sismos_ventanas <- datos |>
  count(zona, mes, name = "conteo_sismos") |>
  group_by(zona) |>
  complete(mes = seq(floor_date(PERIODO_INI, "month"),
                     floor_date(PERIODO_FIN, "month"), by = "month"),
           fill = list(conteo_sismos = 0)) |>
  ungroup()

# 2. Función que devuelve coordenadas de Lorenz + Gini
obtener_lorenz <- function(df) {
  conteos <- sort(df$conteo_sismos)
  curva   <- ineq::Lc(conteos)
  data.frame(p_tiempo = curva$p,
             p_sismos = curva$L,
             gini     = ineq::ineq(conteos, type = "Gini"))
}

# 3. Curvas por zona, con el Gini en la etiqueta de la leyenda 
lorenz_alaska <- sismos_ventanas |>
  filter(zona == "A — Alaska-Aleutianas oriental") |>
  obtener_lorenz()
gini_alaska <- unique(lorenz_alaska$gini)
lorenz_alaska <- lorenz_alaska |>
  mutate(etiqueta = paste0("A — Alaska-Aleutianas (Gini = ",
                           round(gini_alaska, 3), ")"),
         zona = "A — Alaska-Aleutianas oriental")

lorenz_chile <- sismos_ventanas |>
  filter(zona == "B — Sur de Chile") |>
  obtener_lorenz()
gini_chile <- unique(lorenz_chile$gini)
lorenz_chile <- lorenz_chile |>
  mutate(etiqueta = paste0("B — Sur de Chile (Gini = ",
                           round(gini_chile, 3), ")"),
         zona = "B — Sur de Chile")

datos_lorenz <- bind_rows(lorenz_alaska, lorenz_chile)

## Mapear cada etiqueta-con-Gini al color de su zona 
etq <- distinct(datos_lorenz, zona, etiqueta)
col_lorenz <- setNames(col_zonas[etq$zona], etq$etiqueta)

## ---------------------------------------------------------------------
## Gráfico -----
fig_lorenz <- ggplot(datos_lorenz, aes(p_tiempo, p_sismos, color = etiqueta)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40", linewidth = .8) +
  geom_line(linewidth = 1.3) +
  geom_point(data = data.frame(x = c(0, 1), y = c(0, 1)),
             aes(x, y), inherit.aes = FALSE, size = 2.5) +
  scale_color_manual(values = col_lorenz, name = "Zona (índice de Gini)") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  coord_fixed(ratio = 1) +
  labs(title = "Curva de Lorenz de la concentración sísmica",
       subtitle = "Eventos M ≥ 5,0 agrupados por mes (2000–2025)",
       x = "Proporción acumulada de meses ordenados por actividad",
       y = "Proporción acumulada del total de sismos",
       caption = "Fuente: catálogo USGS. Mayor desviación de la diagonal = más concentración temporal.") +
  tema_jids + theme(legend.position = "bottom")
fig_lorenz
ggsave("salidas/figuras/fig_lorenz.png", fig_lorenz,
       width = 8, height = 8, dpi = 300)






#### 10. MAPAS GEORREFERENCIADOS ------------------------------------------------

## Mapa zonas con placas -----
rect_alaska_sf <- filter(rect_sf, zona == "A — Alaska-Aleutianas oriental")
rect_chile_sf  <- filter(rect_sf, zona == "B — Sur de Chile")


mapa_zonas_alaska <- ggplot() +
  geom_sf(data = mundo, fill = "gray95", color = "black") +
  geom_sf(data = placas_alaska, color = "black", linewidth = 1) +
  geom_sf(data = rect_alaska_sf, fill = NA, color = "blue", linewidth = 1.5) +
  annotate("text", x = -150, y = 54, label = "Placa del\nPacífico",
           fontface = "bold", size = 4) +
  annotate("text", x = -155, y = 60, label = "Placa\nNorteamericana",
           fontface = "bold", size = 4) +
  coord_sf(xlim = c(-172, -128), ylim = c(49, 64), expand = FALSE) +
  labs(x = "Longitud (°)", y = "Latitud (°)") +
  theme_bw() + ggtitle("Zona A — Alaska-Aleutianas")

mapa_zonas_chile <- ggplot() +
  geom_sf(data = mundo, fill = "gray95", color = "black") +
  geom_sf(data = placas_chile, color = "black", linewidth = 1) +
  geom_sf(data = rect_chile_sf, fill = NA, color = "red", linewidth = 1.5) +
  annotate("text", x = -78, y = -40, label = "Placa de\nNazca",
           fontface = "bold", size = 4) +
  annotate("text", x = -68, y = -40, label = "Placa\nSudamericana",
           fontface = "bold", size = 4) +
  annotate("text", x = -74, y = -55, label = "Placa\nAntártica",
           fontface = "bold", size = 4) +
  coord_sf(xlim = c(-81, -64), ylim = c(-58, -36), expand = FALSE) +
  labs(x = "Longitud (°)", y = "Latitud (°)") +
  theme_bw() + ggtitle("Zona B — Sur de Chile")

mapa_zonas_placas <- (mapa_zonas_alaska | mapa_zonas_chile) +
  plot_annotation(
    title = "Contexto tectónico y delimitación de las zonas de estudio",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)))
mapa_zonas_placas

ggsave("salidas/figuras/mapa_zonas_placas.png", mapa_zonas_placas,
       width = 13, height = 7, dpi = 300)

## Rangos comunes para que las leyendas coincidan entre zonas -------

rango_mag   <- range(datos$mag,   na.rm = TRUE)
rango_depth <- range(datos$depth, na.rm = TRUE)

##-------------------------------------------------------------------

## Mapa sismos en Chile con placas -----
mapa_chile <- ggplot() +
  geom_sf(data = mundo, fill = "gray95", color = "black") +
  geom_sf(data = placas_chile, color = "black", linewidth = 1) +
  geom_sf(data = chile_sf, aes(size = Mw, color = depth), alpha = 0.7) +
  annotate("text", x = -78, y = -40, label = "Placa de\nNazca",
           size = 4, fontface = "bold") +
  annotate("text", x = -68, y = -40, label = "Placa\nSudamericana",
           size = 4, fontface = "bold") +
  annotate("text", x = -77, y = -55, label = "Placa\nAntártica",
           size = 4, fontface = "bold") +
  coord_sf(xlim = c(-81, -64), ylim = c(-58, -36), expand = FALSE) +
  scale_size_continuous(name = "Magnitud (Mw)", limits = rango_mag) +
  scale_color_gradient(name = "Profundidad (km)",
                       low = "lightblue", high = "darkblue",
                       limits = rango_depth) +
  labs(x = "Longitud (°)", y = "Latitud (°)") +
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Zona B — Sur de Chile")

## Mapa sismos en Alaska con placas -----
mapa_alaska <- ggplot() +
  geom_sf(data = mundo, fill = "gray95", color = "black") +
  geom_sf(data = placas_alaska, color = "black", linewidth = 1) +
  geom_sf(data = alaska_sf, aes(size = Mw, color = depth), alpha = 0.7) +
  annotate("text", x = -145, y = 54, label = "Placa del\nPacífico",
           size = 4, fontface = "bold") +
  annotate("text", x = -150, y = 62, label = "Placa\nNorteamericana",
           size = 4, fontface = "bold") +
  coord_sf(xlim = c(-172, -128), ylim = c(49, 64), expand = FALSE) +
  scale_size_continuous(name = "Magnitud (Mw)", limits = rango_mag) +
  scale_color_gradient(name = "Profundidad (km)",
                       low = "lightblue", high = "darkblue",
                       limits = rango_depth) +
  labs(x = "Longitud (°)", y = "Latitud (°)") +
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Zona A — Alaska-Aleutianas")

## Unión de gráficos -----
mapa_sismos <- (mapa_alaska | mapa_chile) +
  plot_layout(guides = "collect") +
  scale_size_continuous(name = "Magnitud (Mw)", limits = rango_mag) + 
  plot_annotation(
    title = "Distribución espacial de los sismos y principales límites de placas tectónicas",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

mapa_sismos
ggsave("salidas/figuras/mapa_sismos.png", mapa_sismos,
       width = 13, height = 7, dpi = 300)







#### 11. CORROBORACIÓN DE Mc ------------------------------------------------------

## Se considera base con sismos con magnitud >=3.5
ruta_alaska2 <- "Alaska2.csv"   # nombre de archivo en la ruta
ruta_chile2  <- "Chile2.csv"

raw_alaska2 <- readr::read_csv(ruta_alaska2, show_col_types = FALSE) |>
  mutate(zona = "A — Alaska-Aleutianas oriental")

raw_chile2  <- readr::read_csv(ruta_chile2,  show_col_types = FALSE) |>
  mutate(zona = "B — Sur de Chile")

## Unión de ambas bases
datos_raw2 <- bind_rows(raw_alaska2, raw_chile2)

## Depuración de la base (mismas condiciones que la original, pero magnitud 3.5) -----

sismos2 <- datos_raw2 |>
  mutate(
    fecha_hora = ymd_hms(time, tz = "UTC"),
    fecha      = as_date(fecha_hora)
  ) |>
  filter(
    type == "earthquake",                        # excluir no tectónicos (landslide)
    fecha >= PERIODO_INI, fecha <= PERIODO_FIN,  # fecha entre los límites
    mag   >= 3.5                               # magnitud sobre el umbral  
  ) |>
  # Variables derivadas
  mutate(
    anio = year(fecha_hora),
    mes  = floor_date(fecha_hora, "month"),
    prof_categoria = case_when(           # clasificación operativa 
      depth <  70             ~ "Superficial (0–70)",
      depth >= 70 & depth <= 300 ~ "Intermedio (70–300)",
      depth >  300            ~ "Profundo (>300)",
      TRUE                    ~ NA_character_
    ) |> factor(levels = c("Superficial (0–70)",
                           "Intermedio (70–300)",
                           "Profundo (>300)")),
    #fuerte_60 = mag >= 6.0,               # indicadores de eventos fuertes
    #fuerte_65 = mag >= 6.5,
    #fuerte_70 = mag >= 7.0
  ) |>
  arrange(zona, fecha_hora)  # ordena según zona y luego cronológicamente


## Homogeneización -----
sismos2 <- sismos2 |>
  mutate(
    Mw = case_when(
      magType == "mb" ~ 0.85 * mag + 1.03,                         # mb  -> Mw
      magType == "ms" & mag >= 3.0 & mag <= 6.1 ~ 0.67*mag + 2.07, # Ms -> Mw
      magType == "ms" & mag >= 6.2 & mag <= 8.2 ~ 0.99*mag + 0.08, # Ms -> Mw
      magType %in% c("mw","mww","mwb","mwc","mwr") ~ mag,          # ya es Mw
      magType == "ml" ~ mag,                                       # ml conservada
      TRUE ~ mag                                                   # resto
    ),
    fuerte_60 = Mw >= 6.0,               # indicadores de eventos fuertes
    fuerte_65 = Mw >= 6.5,
    fuerte_70 = Mw >= 7.0
  )

## Exportar base procesada con umbral de 3.5 (reproducibilidad) 
readr::write_csv(sismos2, "salidas/base_procesada_35.csv")



## FMD (En conjunto) -----

fmd <- sismos2 %>%
  mutate(M = floor(Mw*10)/10) %>%  
  count(M)

ggplot(fmd, aes(M, n)) +
  geom_col(width = 0.09) +
  labs(
    title = "Distribución Frecuencia–Magnitud del catálogo preliminar",
    subtitle = "Eventos con magnitud igual o superior a 3.5",
    x = "Magnitud",
    y = "Frecuencia"
  ) +
  theme_minimal()

Mc_maxc <- fmd$M[which.max(fmd$n)]
Mc_maxc

fmd_acum <- sismos2 %>%
  mutate(M = floor(Mw*10)/10) %>%
  count(M) %>%
  arrange(desc(M)) %>%
  mutate(N = cumsum(n)) %>%
  arrange(M)


ggplot(fmd_acum,
       aes(M, log10(N))) +
  geom_point() +
  geom_line() +
  labs(
    title = "Distribución Frecuencia–Magnitud Acumulada",
    subtitle = "Representación de la ley de Gutenberg–Richter",
    x = "Magnitud",
    y = expression(log[10](N))
  ) +
  theme_minimal()


## Estimación de b (en cjto) -----

aki_b <- function(M, Mc){
  
  M <- M[M >= Mc]
  
  b <- log10(exp(1)) /
    (mean(M, na.rm = TRUE) - Mc)
  
  return(b)
}

b_est <- data.frame(
  Mc = seq(4.5, 5.2, 0.1),
  b = sapply(
    seq(4.5, 5.2, 0.1),
    function(x) aki_b(sismos2$Mw, x)
  )
)

ggplot(b_est, aes(Mc, b)) +
  geom_line() +
  geom_point(size = 2) +
  geom_vline(xintercept = 5,
             linetype = 2) +
  labs(
    title = "Estabilidad del parámetro b según la magnitud de completitud",
    x = expression(M[c]),
    y = "Parámetro b"
  ) +
  theme_minimal()




## FMD (Por zonas) -----

fmd_region <- sismos2 %>%
  mutate(M = floor(Mw*10)/10) %>%
  count(zona, M)

ggplot(fmd_region,
       aes(M, n, fill = zona)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = col_zonas, name = NULL) +
  labs(
    title = "Distribución Frecuencia–Magnitud por región",
    x = "Magnitud",
    y = "Frecuencia"
  ) +
  tema_jids

gr_zona <- sismos2 %>%
  mutate(M = floor(Mw*10)/10) %>%
  count(zona, M) %>%
  group_by(zona) %>%
  arrange(desc(M), .by_group = TRUE) %>%
  mutate(N = cumsum(n)) %>%
  arrange(M, .by_group = TRUE)

ggplot(gr_zona,
       aes(M, log10(N),
           color = zona)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = col_zonas, name = NULL) +
  labs(
    title = "Ley de Gutenberg–Richter por región",
    x = "Magnitud",
    y = expression(log[10](N))
  ) +
  tema_jids


