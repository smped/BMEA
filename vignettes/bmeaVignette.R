## ---- echo=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(message = FALSE,warning = FALSE)

## ----packagesForGenerationOfVignetteOnly, echo=FALSE-----------------------
library(ggdag)
library(pander)
library(tidyverse)
library(data.tree)

## ----DAG, echo=FALSE, fig.height=6, fig.cap="Full BMEA model expressed as a DAG. All priors are as defined in this figure, with hyperparameters for background signal estimated from the observed data."----
dagify(PM ~ B + S,
       B ~ l,
       B ~ d,
       S ~ e,
       S ~ s,
       e ~ c,
       e ~ p,
       e ~ phi,
       p ~ sp,
       c ~ mu,
       c ~ sm,
       labels = c(PM = "PM[hijk]",
                  B = "B[hijk]",
                  S = "S[hijk]",
                  l = "lambda[il]",
                  d = "delta[il]",
                  e = "eta[hijk]",
                  s = "sigma[S]",
                  c = "c[hi]",
                  p = "p[k]",
                  phi = "phi[hj]",
                  sp = "sigma[p]",
                  mu = "mu[h]",
                  sm = "sigma[mu]")
       ) %>%
  tidy_dagitty() %>%
  mutate(x = case_when(
    name == "PM" ~ 0,
    name == "B" ~ -2,
    name == "S" ~ 2,
    name == "e" ~ -0.5,
    name == "c" ~ -2.5,
    name == "p" ~ -1,
    name == "phi" ~ 0.5,
    name == "mu" ~ -3.5,
    name == "sm" ~ -2.5,
    name == "sp" ~ -1,
    name == "s" ~ 2,
    name %in% c("d", "l") ~ -3.5
  ),
  xend = case_when(
    to == "PM" ~ 0,
    to == "B" ~ -2,
    to == "S" ~ 2,
    to == "e" ~ -0.5,
    to == "c" ~ -2.5,
    to == "p" ~ -1,
    to == "phi" ~ 0.5,
    to == "mu" ~ -3.5,
    to == "sm" ~ -2.5,
    to == "sp" ~ -1,
    to == "s" ~ 2,
    to %in% c("d", "l") ~ -3.5
  ),
  y = case_when(
    name == "PM" ~ 1,
    name %in% c("B", "S") ~ 2,
    name == "e" ~ 3,
    name %in% c("c", "p", "phi") ~ 4,
    name %in% c("mu", "sm", "sp", "s") ~ 5.5,
    name == "d" ~ 1.8,
    name == "l" ~ 2.7
  ),
  yend = case_when(
    to == "PM" ~ 1,
    to %in% c("B", "S") ~ 2,
    to == "e" ~ 3,
    to %in% c("c", "p", "phi") ~ 4,
    to %in% c("mu", "sm", "sp", "s") ~ 5.5,
    to == "d" ~ 1.8,
    to == "l" ~ 2.7
  ),
  shape = name %in% c("d", "l"),
  label = paste0("italic(", label, ")")) %>% 
  ggplot(aes(x, y, xend = xend, yend= yend)) +
  geom_dag_edges() +
  geom_dag_node(aes(shape = shape),
                colour = "white", internal_colour = "black") +
  geom_dag_text(aes(label = label), 
                parse = TRUE, 
                colour = "black",
                size = 5,
                family = "Times") +
  geom_rect(aes(xmin = xmin, xmax = xmax, 
                ymin = ymin, ymax = ymax),
            data = data_frame(ymin = 0.5,
                              ymax = 4.8,
                              xmin = -3,
                              xmax = 2.5), 
            fill = rgb(1,1,1, 0),
            colour = "black",
            inherit.aes = FALSE) +
  geom_text(aes(x, y, label = label),
            data = data_frame(
              label = c("italic(PM[hijk]) == italic(B[hijk] + S[hijk])",
                        "log(italic(S[hijk])) %~%~N(italic(eta[hijk],sigma[S]))",
                        "log(italic(B[hijk])) %~%~N(italic(lambda[il],delta[il]))",
                        "italic(eta[hijk] == c[hi] + p[k] + log(phi[hj]))",
                        "italic(phi[hj] %~%~U(0,1))",
                        "italic(p[k] %~%~N(0, sigma[p]))",
                        "italic(c[hi] %~%~N(mu[h], sigma[mu]))",
                        "italic(sigma[S] %prop% sigma[S]^-1)",
                        "italic(sigma[p] %~%~U(0, 10))",
                        "italic(sigma[mu] %~%~U(0, 5))",
                        "italic(mu[h] %~%~U(0, 2^16))")
              ) %>%
              mutate(x = 2.8,
                     y = seq(0.5, by = 0.52, length.out = nrow(.))),
            inherit.aes = FALSE, 
            family = "Times",
            parse = TRUE,
            size = 5,
            hjust = 0,
            colour = "black") +
  scale_shape_manual(values = c(21, 22)) +
  scale_x_continuous(limits = c(-4, 5)) +
  scale_y_continuous(limits = c(0.5, 7)) +
  guides(shape = FALSE) +
  theme_void()

## ----minimalPath, echo=FALSE, results='markup'-----------------------------
path <- c(
    "parentDirectory/annotationData/chipTypes/chipType/chipType.cdf",
    "parentDirectory/annotationData/chipTypes/chipType/chipType_bgProbes.bgp",
    "parentDirectory/annotationData/chipTypes/chipType/chipType.probe_tab",
    "parentDirectory/rawData/exptName/chipType/File1.CEL",
    "parentDirectory/rawData/exptName/chipType/File2.CEL"
)
data.tree::as.Node(data.frame(pathString = path))

## ----myPath, echo=FALSE, results='markup'----------------------------------
path <- c(
  "parentDirectory/annotationData/chipTypes/HuEx-1_0-st-v2/HuEx-1_0-st-v2.cdf",
  "parentDirectory/annotationData/chipTypes/HuEx-1_0-st-v2/HuEx-1_0-st-v2.r2.antigenomic.bgp",
  "parentDirectory/annotationData/chipTypes/HuEx-1_0-st-v2/HuEx-1_0-st-v2,Custom.cdf",
  "parentDirectory/annotationData/chipTypes/HuEx-1_0-st-v2/HuEx-1_0-st-v2,Custom_probe_tab",
  paste0("parentDirectory/rawData/myExpt/HuEx-1_0-st-v2/File", 1:8, ".CEL")
)
data.tree::as.Node(data.frame(pathString = path))

## --------------------------------------------------------------------------
library(BMEA)
library(limma)
library(snow)
library(magrittr)
library(tidyverse)

## ---- echo=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(eval = FALSE)

## --------------------------------------------------------------------------
#  chipType <- "HuEx-1_0-st-v2"
#  affyCdf <- AffymetrixCdfFile$byChipType(chipType)
#  myCdf <- AffymetrixCdfFile$byChipType(chipType, tags="Custom")

## --------------------------------------------------------------------------
#  exptName <- "myExpt"
#  cs <- AffymetrixCelSet$byName(exptName, cdf=affyCdf)

## ----qn--------------------------------------------------------------------
#  qn <- QuantileNormalization(cs)
#  csN <- process(qn, verbose=verbose)

## ----bgParam---------------------------------------------------------------
#  bgParam <- fitBackgroundParameters(csN, cdf=affyCdf, bgProbes="r2.antigenomic.bgp", method="MAT")

## --------------------------------------------------------------------------
#  bgBins <- defineMatBins(bgParam)

## --------------------------------------------------------------------------
#  c("lambda", "delta") %>%
#    lapply(function(x){
#      set_rownames(bgBins[[x]], paste0("Bin", 1:nrow(bgBins[[x]]))) %>%
#        as.data.frame(stringsAsFactors = FALSE) %>%
#        rownames_to_column("Bin") %>%
#        gather("Array", "Estimate", -Bin) %>%
#        mutate(Parameter = paste0("hat(", x, ")[il]"),
#               Bin = factor(Bin, levels = paste0("Bin", 1:nrow(.))))
#    }) %>%
#    bind_rows() %>%
#    as_tibble() %>%
#    ggplot(aes(Bin, Estimate)) +
#    geom_boxplot() +
#    facet_wrap(~Parameter, scales = "free", ncol = 1, labeller = label_parsed) +
#    theme_bw() +
#    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## --------------------------------------------------------------------------
#  setCdf(csN, myCdf)

## --------------------------------------------------------------------------
#  pmSequenceFile <- "HuEx-1_0-st-v2,Custom.probe_tab"
#  bgCelSet <- assignBgPriors(csN,
#                             seqFile=pmSequenceFile,
#                             bgBins=bgBins,
#                             bgParam=bgParam,
#                             overWrite=TRUE)

## ----contrasts-------------------------------------------------------------
#  n <- length(csN)
#  conditions <- rep(c("A", "B"), each = floor((n+1)/2))[seq_len(n)]
#  conditions <- factor(conditions, levels = c("A", "B"))
#  contrastMatrix <- makeContrasts(BVsA = B - A, levels = conditions)

## ----mcmcParam, eval=TRUE--------------------------------------------------
mcmcParam <- list(nChains = 3L, nIter = 12000L, nBurnin = 6000L, nThin = 6L)

## ----units-----------------------------------------------------------------
#  units <- data_frame(unitID = getUnitNames(myCdf), unit = seq_along(unitID))

## ----fitBmea.Snow----------------------------------------------------------
#  nCores <- 3
#  fitUnits <- units$unit[seq_len(nCores*20)]
#  cl <- makeCluster(nCores, type="SOCK")
#  bmeaFit <- fitBmea.Snow(celSet = csN,
#                          bgCelSet = bgCelSet,
#                          cl = cl,
#                          units = fitUnits,
#                          batchSize = 10,
#                          conditions = conditions,
#                          contMatrix = contrastMatrix,
#                          paramToWrite = c("c", "mu", "phi"),
#                          mcmcParam = mcmcParam)
#  stopCluster(cl)

## ----mergeNodes------------------------------------------------------------
#  mergeNodes(celSet=csN,
#             as.list(names(bmeaFit$units)),
#             paramToWrite=c("c", "mu","phi"))

## ----clearNodes------------------------------------------------------------
#  clearNodes(names(bmeaFit$units))

## ----csLogFC---------------------------------------------------------------
#  csLogFC <- AffymetrixCelSetList(csN, type="contrast", tags="logFC")

## ----bmeaLogFC-------------------------------------------------------------
#  bmeaLogFC <- extractBmeaArray(csLogFC)[, c("2.5%", "mean", "97.5%", "B"),"BVsA"] %>%
#    as.data.frame() %>%
#    rownames_to_column("unitID") %>%
#    as_data_frame() %>%
#    mutate(mean = mean / log(2),
#           `2.5%` = `2.5%` / log(2),
#           `97.5%` = `97.5%` / log(2)) %>%
#    filter(is.finite(mean)) %>% # Get rid of genes not fitted
#    arrange(desc(abs(B)))

## ----csMu------------------------------------------------------------------
#  csMu <- AffymetrixCelSetList(csN, type="model", tags="mu")
#  mu <- extractBmeaArray(csMu, units = fitUnits)[, "mean", c("A", "B")] %>%
#    as.data.frame() %>%
#    rownames_to_column("unitID") %>%
#    as_tibble() %>%
#    gather("Array", "Mu", -unitID)

## ----muTopQ----------------------------------------------------------------
#  muTopQ <- mu %>%
#    filter(Mu > 0) %>%
#    group_by(unitID) %>%
#    summarise(Mu_mean = mean(Mu)) %>%
#    filter(Mu_mean > quantile(Mu_mean, probs = 0.75)) %>%
#    left_join(units)

## ----logFCLowQ-------------------------------------------------------------
#  logFCLowQ <- filter(bmeaLogFC,
#                      abs(mean) < quantile(abs(mean), probs = 0.25))

## ----candUnits-------------------------------------------------------------
#  candidateUnits <- intersect(logFCLowQ$unitID, muTopQ$unitID)

## ----phiLogFC--------------------------------------------------------------
#  phiLogFC <- extractBmeaArray(csPhiLogFC,
#                               units = candidateUnits,
#                               firstOnly = FALSE) %>%
#    magrittr::extract(, c("2.5%", "50%", "97.5%", "B"),"BVsA") %>%
#    as.data.frame() %>%
#    rownames_to_column("groupID") %>%
#    as_data_frame() %>%
#    rename(median = `50%`) %>%
#    mutate(median = median / log(2),
#           `2.5%` = `2.5%` / log(2),
#           `97.5%` = `97.5%` / log(2))

## --------------------------------------------------------------------------
#  phiLogFC %>%
#    separate(groupID, into = c("unitID", "groupID"), sep = "\\.") %>%
#    arrange(desc(abs(B)))

## --------------------------------------------------------------------------
#  fitUgc <- getUnitGroupCellMap(myCdf, units = fitUnits, retNames = TRUE)
#  fitPM <- getIntensities(csN, indices = fitUgc$cell) %>%
#    set_colnames(celNames) %>%
#    as.data.frame() %>%
#    split(f = fitUgc$unit) %>%
#    parallel::mclapply(as.matrix, mc.cores = 3)
#  fitLambda <- getIntensities(bgCelSet$lambda, indices = fitUgc$cell) %>%
#    set_colnames(celNames) %>%
#    log() %>%
#    as.data.frame() %>%
#    split(f = fitUgc$unit) %>%
#    parallel::mclapply(as.matrix, mc.cores = 3)
#  fitDelta <- getIntensities(bgCelSet$delta, indices = fitUgc$cell) %>%
#    set_colnames(celNames) %>%
#    log() %>%
#    as.data.frame() %>%
#    split(f = fitUgc$unit) %>%
#    parallel::mclapply(as.matrix, mc.cores = 3)
#  fitZ <- names(fitPM) %>%
#    parallel::mclapply(function(x){
#      zScore(fitPM[[x]], fitLambda[[x]], fitDelta[[x]],
#             exons = droplevels(filter(fitUgc, unit ==x))$group)
#    }, mc.cores = 3) %>%
#    set_names(names(fitPM))
#  exonZ <- fitZ %>%
#    lapply(function(x){x$exon}) %>%
#    unlist %>%
#    as.data.frame() %>%
#    set_names("Z") %>%
#    rownames_to_column("unitID") %>%
#    as_data_frame() %>%
#    mutate(groupID = gsub(".+\\.(ENSG[0-9_]+)", "\\1", unitID),
#           unitID = gsub("(ENSG[0-9]+)\\..+", "\\1", unitID)) %>%
#    dplyr::select(unitID, groupID, Z)

## --------------------------------------------------------------------------
#  phiLogFC %>%
#    separate(groupID, into = c("unitID", "groupID"), sep = "\\.") %>%
#    left_join(exonZ) %>%
#    arrange(desc(abs(B)))

## --------------------------------------------------------------------------
#  unit <- 1
#  myFit <- fitBmeaSingle(csN, bgCelSet, unit, conditions, contrastMatrix,
#                         mcmcParam = mcmcParam)

## --------------------------------------------------------------------------
#  myFit$logFC %>%
#    as.data.frame() %>%
#    rownames_to_column("Comparison") %>%
#    ggplot(aes(`50%`, 1)) +
#    geom_point() +
#    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`)) +
#    geom_vline(xintercept = 0, colour = "blue", linetype = 2) +
#    labs(x = "95% CPI",
#         y = c()) +
#    facet_wrap(~Comparison) +
#    theme_bw() +
#    theme(axis.text.y = element_blank(),
#          axis.title.y = element_blank())

## --------------------------------------------------------------------------
#  myFit$phiLogFC$BVsA %>%
#    as.data.frame() %>%
#    rownames_to_column("groupID") %>%
#    dplyr::select(groupID, `2.5%`, median = `50%`, `97.5%`) %>%
#    left_join(exonZ) %>%
#    mutate(groupID = str_extract(groupID, "_[0-9]+"),
#           candidate = Z > 50) %>%
#    ggplot(aes(median, groupID, colour = candidate))+
#    geom_point() +
#    geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`)) +
#    geom_vline(xintercept = 0, colour = "red", linetype = 2) +
#    scale_colour_manual(values = c("grey50", "green")) +
#    labs(x = "95% CPI",
#         colour = "Z > 50") +
#    theme_bw()

## ----sessionInfo, echo=FALSE, eval=TRUE------------------------------------
pander(sessionInfo())

