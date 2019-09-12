all <- read.csv("~/Box Sync/R_Scripts&Data/20171218JAP_iceberg/R-output/allMASTER.csv", header = TRUE, stringsAsFactors = FALSE)
pt <- read.csv("~/Box Sync/R_Scripts&Data/MolMassCalc/20190822JAP_periodictable.csv", header = TRUE, stringsAsFactors = FALSE)
pt$monoisotopic.mass <- as.numeric(substring(pt$monoisotopic.mass, 1, 6))
adducts <- read.csv("~/Box Sync/R_Scripts&Data/MolMassCalc/20190822JAP_adducts.csv", header = TRUE, stringsAsFactors = FALSE)
f <- list()

#extracts formula, parses element and number of atoms into lists
for (i in 1:nrow(all)) {
  f.extract <- function(formula) { 
    # pattern to match the initial chemical 
    # assumes chemical starts with an upper case and optional lowercase followed by zero or more digits
    first <- "^([[:upper:]][[:lower:]]?)([0-9]*).*"
    last <- "^[[:upper:]][[:lower:]]?[0-9]*(.*)" # inverse of above to remove the initial chemical 
    result <- list()
    extract <- formula
    # repeat as long as there is data 
    while ((start <- nchar(extract)) > 0) {
      chem <- sub(first, '\\1 \\2', extract)
      extract <- sub(last, '\\1', extract) 
      # if the number of characters is the same, then there was an error 
      if (nchar(extract) == start){
        warning("Invalid formula:", formula)
        return(NULL) 
      }
      # append to the list
      result[[length(result) + 1L]] <- strsplit(chem, ' ')[[1]]
    }
    result
  }
  f[[i]] <- f.extract(all[i,"Neutral.Formula"])
}

# j is formula in f, k is element in formula
for (j in 1:length(f)) {
  for (k in 1:length(f[[j]])) {
    a <- which((f[[j]][[k]][1] == pt$symbol), arr.ind = TRUE)
    if (is.na(f[[j]][[k]][2]) | f[[j]][[k]][2] == '') {
      f[[j]][[k]][2] <- 1
    } else f[[j]][[k]][2] <- f[[j]][[k]][2]
    f[[j]][[k]][3] <- as.numeric(f[[j]][[k]][2])*as.numeric(pt[a, "monoisotopic.mass"])
  }
}

# # j is formula in f, k is element in formula_mass defect
# for (j in 1:length(f)) {
#   for (k in 1:length(f[[j]])) {
#     a <- which((f[[j]][[k]][1] == pt$symbol), arr.ind = TRUE)
#     if (is.na(f[[j]][[k]][2])) {
#       f[[j]][[k]][2] <- 1
#     } else f[[j]][[k]][2] <- f[[j]][[k]][2]
#     f[[j]][[k]][3] <- as.numeric(f[[j]][[k]][2])*pt[a, "accepted.defect"]
#   }
# }

#calculates neutral molecular mass
for (j in 1:length(f)) {
  for (k in 1:length(f[[j]])) {
    all[j, "expected.mass"] <- Reduce(sum, as.numeric(sapply(f[[j]], "[[", 3)))
  }
}

#calculates m/z
for (m in 1:nrow(all)) {
  a <- which(all[m, "Ion.Species"] == adducts$species, arr.ind = TRUE)
  all[m, "expected.mz"] <- as.numeric((all[m, "expected.mass"] * adducts[a, "m"] + adducts[a, "adduct"]))/abs(all[m, "Charge"])
}

# #calculates m/z_mass defect
# for (m in 1:nrow(all)) {
#   a <- which(all[m, "Ion.Species"] == adducts$species, arr.ind = TRUE)
#   all[m, "expected.mz"] <- as.numeric((all[m, "expected.mass"] * adducts[a, "m"] + adducts[a, "adduct.defect"]))/abs(all[m, "Charge"])
# }

#match_0.1% Da window
for (m in 1:nrow(all)) {
  if (data.table::between(all[m, "mz"], as.numeric(all[m, "expected.mz"])-0.001*as.numeric(all[m, "expected.mz"]), as.numeric(all[m, "expected.mz"])+0.001*as.numeric(all[m, "expected.mz"]))) {
    all[m, "mz.check"] <- "pass"
  } else {all[m, "mz.check"] <- "fail"}
}

# #match_half Da window
# for (m in 1:nrow(all)) {
#   if (data.table::between(all[m, "mz"], as.numeric(all[m, "expected.mz"])-0.5, as.numeric(all[m, "expected.mz"])+0.5)) {
#     all[m, "mz.check"] <- "pass"
#   } else {all[m, "mz.check"] <- "fail"}
# }

#fail subset
fail <- subset(all[which(all$mz.check == "fail"),])
fail <- fail[, c("Compound", "Neutral.Formula", "CAS", "InChi", "mz", "Charge", "Ion.Species", "Sources", "expected.mass", "expected.mz", "mz.check")]
fail <- fail[order(fail$Sources, fail$Neutral.Formula),]
write.csv(fail, "/Users/JAPicache/Box Sync/R_Scripts&Data/fail.csv")
