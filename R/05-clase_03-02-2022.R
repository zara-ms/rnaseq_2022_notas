## MODELOS ESTADÍSTICOS

# Regresiones lineales: por cada incremento de 1 en el eje x aumenta B1 el valor
# de y

## ?model.matrix
# Sintaxis de fórmula: Y ~ x1 + x2
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat

colnames(mat)

# lm es la función que nos ayuda a calcular una regresión lineal
summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))

## Paquete: ExploreModelMatrix.
# Interfaz gráfica que nos ayuda a ver los modelos de regresión lineal

library("ExploreModelMatrix")

## Datos de ejemplo
(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))

## Creemos las imágenes usando ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)

## Veamos las imágenes
cowplot::plot_grid(plotlist = vd$plotlist)

# De forma interactiva podemos correr:
## Usaremos shiny otra ves
# Nos genera una página con la info. util para refrescar la memoria o explicar
# lo que sucede
app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment
)
if (interactive()) shiny::runApp(app)
