## USO DE UN PAQUETE DE BIOCONDUCTOR
# Paquete: SummarizedExperiment

# Uso de tablas para la información de un experimento: genes, cuentas, IDs, etc.
# Almacenaje en 3 familias de tablas: colData, rowRanges (info. de genes, cada
# renglon = gen, cada columna = muestra. Se obtiene un objeto GRanges), assay(s).

# Objeto GRanges:
# Contiene a seqnames donde se almacenan cromosomas (se almacena
# en un objeto Run Length Encoding (Rle) para almacenar la info. de manera
# eficiente con un objeto más pequeño). La info. de las coordenadas se almacena
# en un objeto Interval Ranges (IRanges) con posiciones de inicio y final así
# como los nombres. En información strand (Rle) se guarda la cadena (-, + o *),
# útil para saber cuál información se sobrelapa.
# Se puede agregar más información como score, GC, etc.

# Paquete: rtracklayer usado para leer archivos con info. de genes.

# ------------------------------------------------------------------------------
## EJERCICIO 1

## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment --> Pedir ayuda

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6

## Números al azar de cuentas
set.seed(20210223)  # Especificar semilla
## Números al azar siguiendo una distribución uniforme
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),  # Longitud de 100 pb c/u
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))

## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
# Objeto summarizedExperiment
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse


## Funciones para manejar nuestro objeto
## Número de genes y muestras
dim(rse)

## IDs de nuestros genes y muestras
dimnames(rse)

## Nombres de tablas de cuentas que tenemos (RPKM, CPM, counts, logcounts, etc)
assayNames(rse)

## El inicio de nuestra tabla de cuentas
head(assay(rse))

## Información de los genes en un objeto de Bioconductor
rowRanges(rse)

## Tabla con información de los genes
rowData(rse) # es idéntico a 'mcols(rowRanges(rse))'

## Tabla con información de las muestras
colData(rse)

# Explicar que sucede aquí
## Comando 1
rse[1:2, ]
# Explicación: Solamente se toman dos genes (renglones) y todas las columnas se
# conservan

## Comando 2
rse[, c("A", "D", "F")]
# Explicación: Solamente se conservan las columnas A, D y F pero se mantienen todos
# los genes.
# ------------------------------------------------------------------------------

# Paquete: iSEE
# Usado para explorar de manera interactiva a summarizedExperiment

## Explora el objeto rse de forma interactiva
library("iSEE")
iSEE::iSEE(rse)

## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
iSEE::iSEE(sce_layer)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

## DATOS DE RNA-SEQ
# ReCount, recount y recount3

## Load recount3 R package
library("recount3")

## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()

## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)
## Crea un objetio de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)

## Explora el objeto RSE
rse_gene_SRP009615

## Explora los proyectos disponibles de forma interactiva
proj_info_interactive <- interactiveDisplayBase::display(human_projects)
## Selecciona un solo renglón en la tabla y da click en "send".

## Aquí verificamos que solo seleccionaste un solo renglón.
stopifnot(nrow(proj_info_interactive) == 1)
## Crea el objeto RSE
rse_gene_interactive <- create_rse(proj_info_interactive)

## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)

## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]


