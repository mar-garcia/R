# dt_preparation

dt_preparation <- function(polarity = c("POS", "NEG"),
                           class, mznoise){
  require(xcms)
  
  load(paste0("data/RData/IS_data_XCMS_", polarity, ".RData"))
  
  data <- featureValues(xdata, method = "sum", value = "into")
  data[is.na(data)] <- 0
  
  features <- data.frame(featureDefinitions(xdata))
  features$polarity <- polarity
  
  class.levels <- levels(factor(pData(xdata)[, class]))
  for(i in seq(length(class.levels))){
    tmp1 <- rowMeans(data[, pData(xdata)[, class] == class.levels[i]])
    features <- cbind(features, tmp1)
    colnames(features)[ncol(features)] <- paste0("mean_", class.levels[i])
    
    tmp2 <- rowMeans(data[, pData(xdata)[, class] != class.levels[i]])
    features <- cbind(features, tmp2)
    colnames(features)[ncol(features)] <- paste0("mean_", class.levels[i], "_NO")
    
    tmp3 <- tmp1 / tmp2
    tmp3[tmp3 == "Inf"] <- tmp1[tmp3 == "Inf"]
    features <- cbind(features, tmp3)
    colnames(features)[ncol(features)] <- paste0("prop_", class.levels[i])
  }
  
  features$mean_max <- apply(
    features[,paste0("mean_", class.levels)], 1, max)
  
  # Leave noisy peaks for last
  features$mean_max2 <- features$mean_max
  tmp <- min(features$mean_max)
  for(i in seq(length(mznoise))){
    idx <- which((features$mzmed > mznoise[i] - 0.01) & 
                   (features$mzmed < mznoise[i] + 0.01))
    if(length(idx) > 0){
      features$mean_max2[idx] <- features$mean_max2[idx] / tmp
    }
  }
  
  features <- features[order(features$mean_max2, decreasing = T), ]
  rownames(features) <- gsub("FT", substring(polarity, 1, 1), rownames(features))
  rownames(data) <- gsub("FT", substring(polarity, 1, 1), rownames(data))
  colnames(data) <- gsub(paste0("_", polarity, "_FS.mzData"), "", colnames(data))
  
  datax <- list(
    xdata = xdata,
    data = data,
    features = features,
    polarity = polarity
  )
  
  return(datax)
}


# ft_grouping

ft_grouping <- function(datax){
  
  require(CluMSID)
  require(CompoundDb)
  require(Rdisop)
  
  neutral_losses <- data.frame(
    formula = c("H2O", "C5H12NO5P", "C8H15NO4", "C9H17O10P", "C16H32O2", 
                "C19H38O4"),
    massdiff = NA
  )
  for(i in seq(nrow(neutral_losses))){
    neutral_losses$massdiff[i] <- getMolecule(neutral_losses$formula[i])$exactmass
  }
  
  xdata <- datax[["xdata"]]
  data <- datax[["data"]]
  features <- datax[["features"]]
  
  load(paste0("../tissues/data/RData/MS2_library_", datax[["polarity"]], ".RData"))
  
  features$FG <- NA
  features$MS2 <- NA
  
  y.n <- 0
  
  for(y in seq(nrow(features))){
    if(is.na(features$FG[y])){ # continue only if feature "y" still has not 
      # assigned to any "FG"
      y.n <- y.n + 1 
      y.features <- features[is.na(features$FG), ] # subset the feature  matrix taking only
      # those still have not assigned to any "FG"
      y.ft <- which(rownames(y.features) == rownames(features)[y])
      
      # RT range: get the co-eluting features
      y.features <- y.features[(y.features$rtmed > (y.features$rtmed[y.ft] - 10)) & 
                                 ( y.features$rtmed < (y.features$rtmed[y.ft] + 10)), ]
      if(nrow(y.features) > 1){
        
        # intensity correlation
        y.data <- data[rownames(data) %in% rownames(y.features), ]
        y.cor <- cor(t(y.data), y.data[rownames(features)[y],])
        y.features <- merge(y.features,y.cor, by = "row.names")
        colnames(y.features)[ncol(y.features)] <- "corr_int"
        rownames(y.features) <- y.features$Row.names
        y.features <- y.features[,-1]
        y.features <- y.features[y.features$corr_int > 0.7, ]
        
        # peak-shape correlation
        y.xdata <- filterFile(xdata, which.max(data[rownames(features)[y],]))
        y.chr <- chromatogram(y.xdata,
                              mz = features$mzmed[y] + 0.01 * c(-1, 1),
                              rt = features$rtmed[y] + 10 * c(-1,1),
                              aggregationFun = "max")
        y.cor <- c()
        for(i in seq(nrow(y.features))){
          y.chr2 <- chromatogram(y.xdata,
                                 mz = y.features$mzmed[i] + 0.01 * c(-1, 1),
                                 rt = y.features$rtmed[i] + 10 * c(-1,1),
                                 aggregationFun = "max")
          y.cor <- c(y.cor, correlate(y.chr[[1]], y.chr2[[1]]))
          names(y.cor)[length(y.cor)] <- rownames(y.features)[i]
        }
        y.cor[is.na(y.cor)] <- 0
        y.features <- merge(y.features,y.cor, by = "row.names")
        colnames(y.features)[ncol(y.features)] <- "corr_ps"
        rownames(y.features) <- y.features$Row.names
        y.features <- y.features[,-1]
        
        # MS2 associations
        y.features <- y.features[y.features$corr_ps > 0.5, ]
        y.features$MS2x <- NA
        for(i in seq(nrow(y.features))){
          if(y.features$polarity[i] == "POS"){
            smbl <- "+"
          } else if(y.features$polarity[i] == "NEG"){
            smbl <- "-"
          }
          
          y.ms2 <- findFragment(ms2list, y.features$mzmed[i])
          y.ms2 <- getSpectrum(y.ms2, "rt", y.features$rtmed[i], rt.tol = 10)
          for(j in seq(nrow(y.features))){
            if(length(y.ms2) == 1){
              tmp <- abs(y.ms2@precursor - y.features$mzmed[j]) < 0.01
            } else if(length(y.ms2) > 1){
              tmp <- length(getSpectrum(y.ms2, "precursor", 
                                        y.features$mzmed[j], mz.tol = 0.01)) > 0
            } else if(length(y.ms2) == 0){
              tmp <- FALSE
            }
            if(is.na(y.features$MS2x[i])){
              y.features$MS2x[i] <- tmp
            } else {
              y.features$MS2x[i] <- paste0(y.features$MS2x[i], "; ", tmp)
            }
            
          }
        }
        y.features$MS2x <- grepl("TRUE", y.features$MS2x)
        
        for(i in which(y.features$MS2x)){
          y.ms2 <- findFragment(ms2list, y.features$mzmed[i])
          y.ms2 <- getSpectrum(y.ms2, "rt", y.features$rtmed[i], rt.tol = 10)
          if(length(y.ms2) == 1){
            y.prec <- y.ms2@precursor
          } else if(length(y.ms2) > 1){
            y.prec <- c()
            for(j in seq(length(y.ms2))){
              y.prec <- c(y.prec, y.ms2[[j]]@precursor)
            }
            y.prec <- matchWithPpm(y.features$mzmed, y.prec, ppm = 10)
            y.prec <- lapply(y.prec, length)
            y.prec <- y.features$mzmed[which.max(unlist(y.prec))] 
          }
          y.nl <- y.prec - y.features$mzmed[i]
          y.nl <- neutral_losses$formula[
            (neutral_losses$massdiff > (y.nl - 0.01)) & 
              (neutral_losses$massdiff < (y.nl + 0.01))]
          if(length(y.nl) > 0){
            if(is.na(y.features$MS2[i])){
              y.features$MS2[i] <- paste0("[", round(y.prec), "-", y.nl, "]", smbl)
            } else {
              y.features$MS2[i] <- paste(
                y.features$MS2[i], paste0("[", round(y.prec), "-", y.nl, "]", 
                                          smbl), collapse = "; ")
            }
          }
        }
        
        y.features <- y.features[
          y.features$corr_ps > 0.7 |
            (y.features$corr_ps > 0.5 & y.features$MS2x), ]  
      }
      features$FG[rownames(features) %in% rownames(y.features)] <- y.n
      features[rownames(y.features), "MS2"] <- y.features[, "MS2"]
    }
  }
  
  features$FG <- paste0(substr(datax[["polarity"]], 1, 1), 
                        formatC(features$FG, width = nchar(nrow(features)), 
                                flag = "0"))
  
  datax <- list(
    xdata = xdata,
    data = data,
    features = features,
    polarity = polarity
  )
  
  return(datax)
}



# isotopologues

isotopologues <- function(features){
  
  features <- features[order(features$mean_max2, decreasing = T), ]
  features$isotopes <- NA
  
  z.polarities <- c("POS", "NEG")
  
  for(z in seq(length(z.polarities))){
    if(z.polarities[z] == "POS"){smbl <- "+"
    } else if(z.polarities[z] == "NEG"){smbl <- "-"}
    z.features <- features[features$polarity == z.polarities[z], ]
    z.FG <- levels(factor(z.features$FG))
    z.FG <- z.FG[!grepl("NA", z.FG)]
    for(y in seq(length(z.FG))){
      y.features <- z.features[z.features$FG == z.FG[y], ]
      for(i in seq(nrow(y.features))){
        if(is.na(y.features$isotopes[i])){
          for(j in seq(4)){
            y.features$isotopes[
              unlist(matchWithPpm((y.features$mzmed[i] + 1.003355*j), 
                                  y.features$mzmed, ppm = 10))
            ] <- paste0("(", j, ")13C[", round(y.features$mzmed[i]), "]", smbl)
          }
          y.features$isotopes[
            unlist(matchWithPpm((y.features$mzmed[i] - 1.007276), 
                                y.features$mzmed, ppm = 10))
          ] <- paste0("M[", round(y.features$mzmed[i]), "]", smbl)
          if(length(grep(round(y.features$mzmed[i]), y.features$isotopes)) > 0){
            y.features$isotopes[i] <- paste0("[", round(y.features$mzmed)[i], "]", smbl)
          } 
        }
      }
      features[rownames(y.features), "isotopes"] <- y.features[, "isotopes"]
    }
  }
  
  return(features)
}



# fg_grouping

fg_grouping <- function(data, features){
  
  require(Rdisop)
  
  massneg <- c(-1.007276, getMolecule("Cl")$exactmass, 
               c(-1.007276 + getMolecule("HCOOH")$exactmass))
  names(massneg) <- c("[M-H]-", "[M+Cl]-", "[M-H+HCOOH]-")
  masspos <- c(1.007276, getMolecule("NH4")$exactmass, getMolecule("Na")$exactmass,
               getMolecule("K")$exactmass, getMolecule("C2H8N")$exactmass)
  names(masspos) <- c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[M+C2H8N]+")
  massdif <- outer(masspos, massneg, "-")
  colnames(massdif) <- names(massneg)
  rownames(massdif) <- names(masspos) 

  
  features$FGx <- NA
  features$annotation <- NA
  y.n <- 0
  for(y in seq(nrow(features))){
    if(is.na(features$FGx[y])){
      y.features <- features[is.na(features$FGx), ]
      y.features <- y.features[(y.features$FG == features$FG[y]) |
                                 y.features$polarity != features$polarity[y],]
      y.features <- y.features[(y.features$rtmed > (features$rtmed[y] - 10)) & 
                                 ( y.features$rtmed < (features$rtmed[y] + 10)), ]
      
      #### FG = 1 ------------------------------------------------------------
      if(length(levels(factor(y.features$FG))) == 1){
        if(nrow(y.features) == 1){
          if(y.features$polarity == "NEG"){
            y.features$annotation <- "[M-H]-"
          } else if(y.features$polarity == "POS"){
            y.features$annotation <- "[M+H]+"
          }
          y.features$FGx <- y.features$FG
          tmp <- y.features
        } else {
          massx <- get(paste0("mass", tolower(y.features$polarity[1])))
          massx <- outer(massx, massx, "-")
          rownames(massx) <- names(get(paste0("mass", tolower(y.features$polarity[1]))))
          colnames(massx) <- names(get(paste0("mass", tolower(y.features$polarity[1]))))
          
          tmp <- y.features[!grepl("13C|M", y.features$isotopes),]
          y.dif <- outer(tmp$mzmed, tmp$mzmed, "-")
          rownames(y.dif) <- round(tmp$mzmed, 4)
          colnames(y.dif) <- round(tmp$mzmed, 4)
          
          y.dif[y.dif <= 0] <- NA
          massx[massx <= 0] <- NA
          
          for(i in seq(nrow(y.dif))){
            for(j in seq(ncol(y.dif))){
              for(a in seq(nrow(massx))){
                for(b in seq(ncol(massx))){
                  if(!is.na(y.dif[i,j]) & !is.na(massx[a,b])){
                    if(dplyr::between(y.dif[i,j], massx[a,b] - 0.01, 
                                      massx[a,b] + 0.01)){
                      tmp$annotation[rownames(y.dif)[i] == 
                                       round(tmp$mzmed, 4)] <- rownames(massx)[a]
                      tmp$annotation[colnames(y.dif)[j] == 
                                       round(tmp$mzmed, 4)] <- colnames(massx)[b]
                    }
                  }
                }
              }
            }
          }
        }
        features[rownames(tmp), "annotation"] <- tmp[, "annotation"]
        features$FGx[
          rownames(features) %in% 
            rownames(features)[features$FG %in% 
                                 unique(y.features$FG)]] <- y.features$FG[1]
        
        #### FG = 2 ----------------------------------------------------------
      } else if(length(levels(factor(y.features$FG))) == 2){
        y.n <- y.n + 1 
        if(cor(data[rownames(y.features)[which(y.features$polarity == "POS")[1]],],
               data[rownames(y.features)[which(y.features$polarity == "NEG")[1]],]
        ) > 0.7){
          tmp <- y.features[!grepl("13C|M", y.features$isotopes),]
          y.dif <- outer(tmp$mzmed[tmp$polarity == "POS"],  
                         tmp$mzmed[tmp$polarity == "NEG"], "-")
          rownames(y.dif) <- round(tmp$mzmed[tmp$polarity == "POS"], 4)
          colnames(y.dif) <- round(tmp$mzmed[tmp$polarity == "NEG"], 4)
          for(i in seq(nrow(y.dif))){
            for(j in seq(ncol(y.dif))){
              for(a in seq(nrow(massdif))){
                for(b in seq(ncol(massdif))){
                  if(dplyr::between(y.dif[i,j], massdif[a,b] - 0.01, 
                                    massdif[a,b] + 0.01)){
                    tmp$annotation[rownames(y.dif)[i] == 
                                     round(tmp$mzmed, 4)] <- rownames(massdif)[a]
                    tmp$annotation[colnames(y.dif)[j] == 
                                     round(tmp$mzmed, 4)] <- colnames(massdif)[b]
                  }
                }
              }
            }
          }
          if(sum(!is.na(tmp$annotation)) > 0){
            features[rownames(tmp), "annotation"] <- tmp[, "annotation"]
            features$FGx[
              rownames(features) %in% 
                rownames(features)[features$FG %in% unique(y.features$FG)]] <- paste0(
                  "X", formatC(y.n, width = nchar(nrow(features)), flag = "0"))
          } else {
            y.features <- y.features[y.features$FG == y.features$FG[1], ]
            massx <- get(paste0("mass", tolower(y.features$polarity[1])))
            massx <- outer(massx, massx, "-")
            rownames(massx) <- names(get(paste0("mass", tolower(y.features$polarity[1]))))
            colnames(massx) <- names(get(paste0("mass", tolower(y.features$polarity[1]))))
            
            tmp <- y.features[!grepl("13C|M", y.features$isotopes),]
            y.dif <- outer(tmp$mzmed, tmp$mzmed, "-")
            rownames(y.dif) <- round(tmp$mzmed, 4)
            colnames(y.dif) <- round(tmp$mzmed, 4)
            
            y.dif[y.dif <= 0] <- NA
            massx[massx <= 0] <- NA
            
            for(i in seq(nrow(y.dif))){
              for(j in seq(ncol(y.dif))){
                for(a in seq(nrow(massx))){
                  for(b in seq(ncol(massx))){
                    if(!is.na(y.dif[i,j]) & !is.na(massx[a,b])){
                      if(dplyr::between(
                        y.dif[i,j], massx[a,b] - 0.01, massx[a,b] + 0.01)){
                        tmp$annotation[rownames(y.dif)[i] == 
                                         round(tmp$mzmed, 4)] <- rownames(massx)[a]
                        tmp$annotation[colnames(y.dif)[j] == 
                                         round(tmp$mzmed, 4)] <- colnames(massx)[b]
                      }
                    }
                  }
                }
              }
            }
            features[rownames(tmp), "annotation"] <- tmp[, "annotation"]
            features$FGx[rownames(features) %in% rownames(y.features)] <- y.features$FG[1]
          }
        } else {
          y.features <- y.features[y.features$FG == y.features$FG[1], ]
          massx <- get(paste0("mass", tolower(y.features$polarity[1])))
          massx <- outer(massx, massx, "-")
          rownames(massx) <- names(get(paste0("mass", tolower(y.features$polarity[1]))))
          colnames(massx) <- names(get(paste0("mass", tolower(y.features$polarity[1]))))
          
          tmp <- y.features[!grepl("13C|M", y.features$isotopes),]
          y.dif <- outer(tmp$mzmed, tmp$mzmed, "-")
          rownames(y.dif) <- round(tmp$mzmed, 4)
          colnames(y.dif) <- round(tmp$mzmed, 4)
          
          y.dif[y.dif <= 0] <- NA
          massx[massx <= 0] <- NA
          
          for(i in seq(nrow(y.dif))){
            for(j in seq(ncol(y.dif))){
              for(a in seq(nrow(massx))){
                for(b in seq(ncol(massx))){
                  if(!is.na(y.dif[i,j]) & !is.na(massx[a,b])){
                    if(dplyr::between(y.dif[i,j], massx[a,b] - 0.01, massx[a,b] + 0.01)){
                      tmp$annotation[rownames(y.dif)[i] == 
                                       round(tmp$mzmed, 4)] <- rownames(massx)[a]
                      tmp$annotation[colnames(y.dif)[j] == 
                                       round(tmp$mzmed, 4)] <- colnames(massx)[b]
                    }
                  }
                }
              }
            }
          }
          features[rownames(tmp), "annotation"] <- tmp[, "annotation"]
          features$FGx[
            rownames(features) %in% 
              rownames(features)[features$FG %in% unique(
                y.features$FG)]] <- y.features$FG[1]
        }
        
        #### FG > 2 ----------------------------------------------------------
      } else if(length(levels(factor(y.features$FG))) > 2){
        y.n <- y.n + 1 
        FT <- c()
        for(i in seq(length(levels(factor(y.features$FG))))){
          FT <- c(FT, rownames(y.features)[y.features$FG == levels(factor(
            y.features$FG))[i]][which.max(y.features$mean_max[
              y.features$FG == levels(factor(y.features$FG))[i]])])
        }
        y.features <- y.features[y.features$FG %in% y.features$FG[
          rownames(y.features) %in% names(cor(t(data[FT,]))[
            ,rownames(y.features)[1]] > 0.7)], ]
        tmp <- y.features[!grepl("13C|M", y.features$isotopes),]
        y.dif <- outer(tmp$mzmed[tmp$polarity == "POS"],  
                       tmp$mzmed[tmp$polarity == "NEG"], "-")
        rownames(y.dif) <- round(tmp$mzmed[tmp$polarity == "POS"], 4)
        colnames(y.dif) <- round(tmp$mzmed[tmp$polarity == "NEG"], 4)
        for(i in seq(nrow(y.dif))){
          for(j in seq(ncol(y.dif))){
            for(a in seq(nrow(massdif))){
              for(b in seq(ncol(massdif))){
                if(dplyr::between(
                  y.dif[i,j], massdif[a,b] - 0.01, massdif[a,b] + 0.01)){
                  tmp$annotation[rownames(y.dif)[i] == 
                                   round(tmp$mzmed, 4)] <- rownames(massdif)[a]
                  tmp$annotation[colnames(y.dif)[j] == 
                                   round(tmp$mzmed, 4)] <- colnames(massdif)[b]
                }
              }
            }
          }
        }
        y.features[rownames(tmp), "annotation"] <- tmp[, "annotation"]
        tmp <- tmp[tmp$polarity != y.features$polarity[1],]
        tmp <- tmp[!is.na(tmp$annotation), ]
        if(length(unique(tmp$FG)) > 1){
          tmp2 <- unique(tmp$FG[duplicated(tmp$FG)])
          if(length(tmp2) == 1){
            y.features <- y.features[(y.features$FG == y.features$FG[1]) |
                                       (y.features$FG == tmp2),]
          } else if(length(tmp2) > 1){
            tmp <- tmp[tmp$FG %in% tmp2, ]
            tmp2 <- tmp$FG[unlist(matchWithPpm((features$mzmed[y]-1.007276*2), 
                                               tmp$mzmed, ppm = 10))]
            y.features <- y.features[(y.features$FG == y.features$FG[1]) |
                                       (y.features$FG == tmp2),]
          }
          
        } else if(length(unique(tmp$FG)) == 0){
          y.features <- y.features[(y.features$FG == y.features$FG[1]), ]
          massx <- get(paste0("mass", tolower(y.features$polarity[1])))
          massx <- outer(massx, massx, "-")
          rownames(massx) <- names(
            get(paste0("mass", tolower(y.features$polarity[1]))))
          colnames(massx) <- names(
            get(paste0("mass", tolower(y.features$polarity[1]))))
          
          tmp <- y.features[!grepl("13C|M", y.features$isotopes),]
          y.dif <- outer(tmp$mzmed, tmp$mzmed, "-")
          rownames(y.dif) <- round(tmp$mzmed, 4)
          colnames(y.dif) <- round(tmp$mzmed, 4)
          
          y.dif[y.dif <= 0] <- NA
          massx[massx <= 0] <- NA
          
          for(i in seq(nrow(y.dif))){
            for(j in seq(ncol(y.dif))){
              for(a in seq(nrow(massx))){
                for(b in seq(ncol(massx))){
                  if(!is.na(y.dif[i,j]) & !is.na(massx[a,b])){
                    if(dplyr::between(
                      y.dif[i,j], massx[a,b] - 0.01, massx[a,b] + 0.01)){
                      tmp$annotation[rownames(y.dif)[i] == 
                                       round(tmp$mzmed, 4)] <- rownames(massx)[a]
                      tmp$annotation[colnames(y.dif)[j] == 
                                       round(tmp$mzmed, 4)] <- colnames(massx)[b]
                    }
                  }
                }
              }
            }
          }
          y.features[rownames(tmp), "annotation"] <- tmp[, "annotation"]
        } else {
          y.features <- y.features[
            (y.features$FG == y.features$FG[1]) |
              (y.features$FG == unique(tmp$FG[!is.na(tmp$annotation)])),
          ]
        } # close "if(length(unique(tmp$FG))...)"
        
        features[rownames(y.features), "annotation"] <- y.features[, "annotation"]
        features$FGx[
          rownames(features) %in% 
            rownames(features)[features$FG %in% unique(y.features$FG)]] <- paste0(
              "X", formatC(y.n, width = nchar(nrow(features)), flag = "0"))
        
      } else{print(y)}
    }
  }
  
  return(features)
}


# identification

identification <- function(features, cmps, rt_d = 60, ppm_d = 10){
  
  require(Rdisop)
  require(CluMSID)
  
  neutral_losses <- data.frame(
    formula = c("H2O", "C5H12NO5P", "C8H15NO4", "C9H17O10P", "C16H32O2", 
                "C19H38O4"),
    massdiff = NA
  )
  for(i in seq(nrow(neutral_losses))){
    neutral_losses$massdiff[i] <- getMolecule(neutral_losses$formula[i])$exactmass
  }
  
  load("../tissues/data/RData/MS2_library_POS.RData")
  ms2list_pos <- ms2list
  load("../tissues/data/RData/MS2_library_NEG.RData")
  ms2list_neg <- ms2list
  
  dt_fg <- data.frame(
    FG = unique(features$FGx),
    mass = NA,
    RT = NA,
    POS = NA,
    NEG = NA
  )
  
  features$assignation <- NA
  
  for(y in seq(nrow(dt_fg))){
    y.features <- features[features$FGx == dt_fg$FG[y], ]
    
    dt_fg$RT[y] <- mean(y.features$rtmed)
    
    massneut <- c()
    for(i in seq(length(which(!is.na(y.features$annotation))))){
      massneut <- c(massneut, 
                    unlist(mz2mass(
                      y.features$mzmed[!is.na(y.features$annotation)][i],
                      y.features$annotation[!is.na(y.features$annotation)][i])))
    }
    massneut <- massneut[!is.na(massneut)]
    if(length(massneut) == 0){
      if(y.features$polarity[which.max(y.features$mean_max)] == "POS"){
        massneut <- y.features$mzmed[which.max(y.features$mean_max)] - 1.007276
      } else if(y.features$polarity[which.max(y.features$mean_max)] == "NEG"){
        massneut <- y.features$mzmed[which.max(y.features$mean_max)] + 1.007276
      }
    } else {
      if((max(massneut) - min(massneut)) > 1){
        tmp <- as.numeric(names(which.max(table(round(massneut)))))
        massneut <- massneut[(massneut > (tmp - 1)) & (massneut < (tmp + 1))]
      }
    }
    dt_fg$mass[y] <- mean(massneut, na.rm = TRUE)
    
    idx <- which(duplicated(y.features$annotation[!is.na(y.features$annotation)]))
    if(length(idx) > 0){
      tmp <- unique(y.features$annotation[!is.na(y.features$annotation)][idx])
      tmp2 <- y.features[!is.na(y.features$annotation), ]
      for(j in seq(length(tmp))){
        tmp3 <- tmp2[tmp2$annotation == tmp[j], ]
        if(grepl("13C", tmp[j])){
          for(i in seq(nrow(tmp3))){
            if(length(unlist(matchWithPpm(
              unlist(mass2mz(mean(massneut), 
                             gsub(".*C", "", tmp3$annotation[i]))) +
              as.numeric(substr(tmp3$annotation[i], 2, 2))*1.003355,
              tmp3$mzmed[i], ppm = 10))) == 0){
              y.features$annotation[rownames(y.features) == rownames(tmp3)[i]] <- NA
            }
          }
        } else {
          for(i in seq(nrow(tmp3))){
            if(length(unlist(matchWithPpm(
              unlist(mass2mz(mean(massneut), tmp3$annotation[i])), tmp3$mzmed[i], 
              ppm = 10))) == 0){
              y.features$annotation[rownames(y.features) == rownames(tmp3)[i]] <- NA
            }
          }
        } # close if(!grepl("13C", tmp[j]))
      }
    } # close if(length(idx) > 0)
    
    for(i in seq(nrow(y.features))){
      if(is.na(y.features$annotation)[i]){
        if(y.features$polarity[i] == "NEG"){
          y.pol <- -1
        } else if(y.features$polarity[i] == "POS"){
          y.pol <- 1
        }
        if(y.features$mzmed[i] < (mean(massneut) - 2)){
          if(y.features$polarity[i] == "POS"){smbl <- "+"
          } else if(y.features$polarity[i] == "NEG"){smbl <- "-"}
          ms2list <- get(paste0("ms2list_", tolower(y.features$polarity[i])))
          y.ms2 <- findFragment(ms2list, y.features$mzmed[i])
          y.ms2 <- getSpectrum(y.ms2, "rt", y.features$rtmed[i], rt.tol = 15)
          if(length(y.ms2) > 0){
            if(length(y.ms2) > 1){
              y.prec <- c()
              for(j in seq(length(y.ms2))){
                y.prec <- c(y.prec, y.ms2[[j]]@precursor)
              }
              if((max(y.prec) - min(y.prec)) > 1){
                tmp <- data.frame(table(round(y.prec)))
                tmp <- tmp[order(tmp$Freq, decreasing = T), ]
                tmp$Var1 <- as.numeric(as.character(tmp$Var1))
                tmp$present <- FALSE
                for(j in seq(nrow(tmp))){
                  tmp$present[j] <- length(unlist(matchWithPpm(mean(
                    y.prec[(y.prec > (tmp$Var1[j] - 1)) & 
                             (y.prec < (tmp$Var1[j] + 1))]),
                    y.features$mzmed[y.features$polarity == 
                                       y.features$polarity[i]], ppm = 10))) > 0
                }
                tmp <- tmp[tmp$present, ]
                if(nrow(tmp) > 0){
                  y.prec <- mean(y.prec[(y.prec > (tmp$Var1[1] - 1)) & 
                                          (y.prec < (tmp$Var1[1] + 1))])
                } else {
                  y.prec <- NA
                }
              } else {y.prec <- mean(y.prec)
              } # close "if((max(y.prec) - min(y.prec)) > 1)"
            } else if(length(y.ms2) == 1){
              y.prec <- y.ms2@precursor
            }
            if(!is.na(y.prec)){
              y.nl <- mean(y.prec) - y.features$mzmed[i]
              y.nl <- neutral_losses$formula[
                (neutral_losses$massdiff > (y.nl - 0.01)) & 
                  (neutral_losses$massdiff < (y.nl + 0.01))]
              if(length(y.nl) == 1){
                y.features$MS2[i] <-paste0("[", round(mean(y.prec)), 
                                           "-", y.nl, "]", smbl)
              } else if(length(y.nl) == 0){
                y.features$MS2[i] <- paste0("[", round(mean(y.prec)), 
                                            "-X]", smbl)
              } else{
                print(paste(y, i))
              } # close  if(length(y.nl) == ...)
            } # close "if(!is.na(y.prec))°
          } # close if(length(y.ms2) > 0)
        } else {
          tmp <- names(unlist(matchWithPpm(
            y.features$mzmed[i], unlist(mass2mz(
              mean(massneut, na.rm = TRUE), adducts(polarity = y.pol))), 
            ppm = 10)))
          if(length(tmp) > 0){
            y.features$annotation[i] <- tmp[1]
          }
        }# close if(!(y.features$mzmed[i] < (mean(massneut) - 2)))
        #} # close "if length(tmp) > 0"
      } # close if(is.na(y.features$annotation[i]))
      
      if(!is.na(y.features$annotation[i])){
        if(!is.na(y.features$isotopes[i])){
          idx <- grep(gsub(".*\\[", "", y.features$isotopes[i]),
                      y.features$isotopes)
          idx <- idx[y.features$polarity[idx] == y.features$polarity[i]]
          if(!is.na(y.features$isotopes[idx][1]) & 
             !is.na(y.features$annotation[idx][1])){
            y.features$assignation[idx] <- gsub(
              paste0("\\", y.features$isotopes[idx][1]), y.features$annotation[i],
              y.features$isotopes[idx])
          } else {
            y.features$assignation[i] <- y.features$annotation[i]
          }
        } else {
          y.features$assignation[i] <- y.features$annotation[i]
        }
      } # close if(!is.na(y.features$annotation[i]))
    } # close "for(i in seq(nrow(y.features)))"
    
    y.features$assignation <- gsub("\\++", "+", y.features$assignation)
    
    idx <- is.na(y.features$assignation) & !is.na(y.features$MS2)
    y.features$assignation[idx] <- y.features$MS2[is.na(y.features$assignation) &
                                                    !is.na(y.features$MS2)]
    
    idx <- is.na(y.features$assignation) & !is.na(y.features$isotopes)
    y.features$assignation[idx] <- y.features$isotopes[idx]
    
    y.features$assignation <- gsub("\\(1)13C", "13C", y.features$assignation)
    y.features$assignation <- gsub("M\\[M\\+H\\]\\+", "\\[M\\]\\+",
                                   y.features$assignation)
    y.features$assignation <- gsub("M\\[M\\+NH4\\]\\+", "\\[M\\]\\+ +NH3",
                                   y.features$assignation)
    
    y.features$assignation2 <- y.features$assignation
    y.features$assignation2[!grepl("M", y.features$assignation2)] <- NA
    annx <- c(ann, levels(factor(y.features$assignation2))[!levels(
      factor(y.features$assignation2)) %in% ann])
    y.features$assignation2 <- factor(y.features$assignation2, levels = annx)
    y.features <- y.features[order(y.features$assignation2, y.features$mzmed), ]
    idx <- y.features$polarity == "POS"
    dt_fg$POS[y] <- paste(paste(round(y.features$mzmed[idx], 4),
                                y.features$assignation[idx]), collapse = "; ")
    idx <- y.features$polarity == "NEG"
    dt_fg$NEG[y] <- paste(paste(round(y.features$mzmed[idx], 4),
                                y.features$assignation[idx]), collapse = "; ")
    
    features[rownames(y.features), "assignation"] <- y.features[, "assignation"]
    
  } #close "for(y in seq(nrow(dt_fg)))"
  dt_fg$POS <- gsub(" NA", "", dt_fg$POS)
  dt_fg$NEG <- gsub(" NA", "", dt_fg$NEG)
  
  dt_fg$compound <- NA
  for(i in seq(nrow(dt_fg))){
    i.cmps <- cmps[unlist(matchWithPpm(dt_fg$mass[i], cmps$mass, ppm = ppm_d)), ]
    if(nrow(i.cmps) > 0){
      i.cmps <- i.cmps[(i.cmps$RT > (dt_fg$RT[i] - rt_d)) & 
                         (i.cmps$RT < (dt_fg$RT[i] + rt_d)), ]
      if(nrow(i.cmps) > 0){
        dt_fg$compound[i] <- paste(i.cmps$ID, collapse = "; ")
      } else {
        dt_fg$compound[i] <- paste("Unknown", sprintf("%.4f", dt_fg$mass[i]))
      }
    } else {
      dt_fg$compound[i] <- paste("Unknown", sprintf("%.4f", dt_fg$mass[i]))
    }
    
  }
  
  dt_fg <- dt_fg[order(dt_fg$RT), c("FG", "compound", "RT", "NEG", "POS")]
  dt_fg$RT <- sprintf("%.2f", round(dt_fg$RT/60, 2))
  
  
  dt_fg$unknown <- grepl("Unknown", dt_fg$compound) | (dt_fg$compound == "")
  dt_fg$duplicated <- dt_fg$compound %in% dt_fg$compound[duplicated(dt_fg$compound)]
  dt_fg$background <- dt_fg$FG %in% unique(
    features$FGx[features$mean_max2 < min(apply(
      features[,c("mean_solv", "mean_blank", "mean_STDmix")], 1, max))])
  dt_fg$color <- 0
  dt_fg$color[!dt_fg$unknown] <- 1
  dt_fg$color[dt_fg$duplicated] <- 2
  dt_fg$color[dt_fg$background] <- 3
  dt_fg$color[!dt_fg$unknown & !dt_fg$duplicated] <- 1
  dt_fg$compound[dt_fg$background & dt_fg$unknown] <- gsub(
    "Unknown", "Background", dt_fg$compound[dt_fg$background & dt_fg$unknown])
  
  datax <- list(
    features = features,
    dt_fg = dt_fg
  )
  
  return(datax)
}


