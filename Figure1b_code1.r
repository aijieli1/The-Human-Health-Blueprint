library(data.table)
library(GenomicSEM)
library(tidyr)
library(dplyr)
library(semPlot)
library(lavaan)
library(ggsci)

semPlotModel_GSEM <- function(gsem.object, est.label = "STD_All"){ 
  object <- gsem.object$results
  object$free <- 0
  object$free[object$op != "~~"] <- seq_len(sum(object$op != "~~"))
  varNames  <- lavaanNames(object, type = "ov")
  factNames <- setdiff(lavaanNames(object, type = "lv"), varNames)
  if (is.null(object$label)) object$label <- rep("", nrow(object))
  object$est <- object[, est.label]
  if (is.null(object$group)) object$group <- ""
  semModel <- new("semPlotModel")
  semModel@Pars <- data.frame(
    label = object$label,
    lhs   = ifelse(object$op %in% c("~","~1"), object$rhs, object$lhs),
    edge  = "--",
    rhs   = ifelse(object$op %in% c("~","~1"), object$lhs, object$rhs),
    est   = object$est,
    std   = NA,
    group = object$group,
    fixed = object$free == 0,
    par   = object$free,
    stringsAsFactors = FALSE
  )
  semModel@Pars$edge[object$op == "~~"]  <- "<->"
  semModel@Pars$edge[object$op == "~*~"] <- "<->"
  semModel@Pars$edge[object$op == "~"]   <- "~>"
  semModel@Pars$edge[object$op == "=~"]  <- "->"
  semModel@Pars$edge[object$op == "~1"]  <- "int"
  semModel@Pars <- semModel@Pars[!object$op %in% c(":=","<",">","==","|"), ]
  semModel@Vars <- data.frame(
    name = c(varNames, factNames),
    manifest = c(varNames, factNames) %in% varNames,
    exogenous = NA,
    stringsAsFactors = FALSE
  )
  semModel@Thresholds <- semModel@Pars[grepl("\\|", semModel@Pars$edge), -(3:4)]
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  semModel@Computed <- FALSE
  semModel@Original <- list(object)
  semModel
}

aid_factor <- readRDS("//1.RDS")

aid_factor$results$STD_Genotype     <- formatC(aid_factor$results$STD_Genotype, digits = 2, format = "f")
aid_factor$results$STD_Genotype_SE  <- formatC(as.numeric(aid_factor$results$STD_Genotype_SE), digits = 2, format = "f")
aid_factor$results$STD_All          <- formatC(aid_factor$results$STD_All, digits = 2, format = "f")

npg_colors <- adjustcolor(pal_npg("nrc")(8), alpha.f = 0.6)
npg_colors <- c("#E64B3580","#4DBBD580","#00A08780","#3C548880","#F39B7F80","#8491B480","#91D1C280")

pdf("//1.pdf", width = 15, height = 12)
semPaths(
  semPlotModel_GSEM(aid_factor),
  as.expression = c("nodes","edges"),
  sizeMan = 6,
  sizeLat = 9,
  curve = 3,
  bg = "white",
  groups = "latents",
  intercepts = FALSE,
  borders = TRUE,
  label.norm = "O",
  what = "diagram",
  whatLabels = "par",
  residuals = TRUE,
  rotation = 2,        
  theme = "colorblind",
  node.width = 0.6,
  edge.label.cex = 0.8,
  color = list(lat = npg_colors),
  asize = 2
)
dev.off()

