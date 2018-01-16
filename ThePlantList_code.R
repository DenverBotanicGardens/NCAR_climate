function (splist, genus = NULL, species = NULL, infrasp = NULL, 
          infra = TRUE, corr = TRUE, diffchar = 2, max.distance = 1, 
          version = "1.1", encoding = "UTF-8", author = TRUE, drop.lower.level = FALSE, 
          file = "", silent = TRUE, repeats = 6) 
{
  splist2 <- NULL
  try(splist2 <- splist, silent = TRUE)
  if (!is.null(splist2) && (!is.null(genus) || !is.null(species) || 
                            !is.null(infrasp))) {
    stop("Argument 'splist' incompatible with arguments 'genus' and 'species'")
  }
  else if (is.null(splist2) && ((is.null(genus) && !is.null(species)) || 
                                (!is.null(genus) && is.null(species)))) {
    stop("Arguments 'genus' and 'species' must be provided")
  }
  else if (is.null(splist2) && !is.null(genus) && !is.null(species)) {
    if (infra == TRUE && !is.null(infrasp)) {
      splist <- paste(genus, species, infrasp)
    }
    else if (infra == FALSE || is.null(infrasp)) {
      splist <- paste(genus, species)
    }
  }
  TPLck2 <- function(d) {
    a <- NULL
    if (silent == FALSE) {
      print(paste("Checking", as.character(d), "in The Plant List"))
    }
    counter <- 0
    a <- try(TPLck(sp = d, infra = infra, corr = corr, diffchar = diffchar, 
                   max.distance = max.distance, version = version, 
                   encoding = encoding, author = author, drop.lower.level = drop.lower.level), 
             silent = FALSE)
    while (class(a) == "try-error" && counter < repeats) {
      a <- try(TPLck(sp = d, infra = infra, corr = corr, 
                     diffchar = diffchar, max.distance = max.distance, 
                     version = version, encoding = encoding, author = author, 
                     drop.lower.level = drop.lower.level), silent = TRUE)
      counter <- counter + 1
    }
    if (class(a) == "try-error") {
      a <- data.frame(Taxon = d, Genus = NA, Hybrid.marker = NA, 
                      Species = NA, Abbrev = NA, Infraspecific.rank = NA, 
                      Infraspecific = NA, Authority = NA, ID = NA, 
                      Plant.Name.Index = NA, TPL.version = NA, Taxonomic.status = NA, 
                      Family = NA, New.Genus = NA, New.Hybrid.marker = NA, 
                      New.Species = NA, New.Infraspecific.rank = NA, 
                      New.Infraspecific = NA, New.Authority = NA, 
                      New.ID = NA, New.Taxonomic.status = NA, Typo = NA, 
                      WFormat = NA, Higher.level = NA, Date = NA, 
                      stringsAsFactors = FALSE)
    }
    invisible(a)
  }
  if (length(splist) < 5) {
    results <- do.call("rbind", lapply(splist, TPLck2))
  }
  else {
    op <- pbapply::pboptions() #when adding progress bar, supports parallel processing
    pbapply::pboptions(type = "txt")
    results <- do.call("rbind", pbapply::pblapply(splist, 
                                                  TPLck2))
    pbapply::pboptions(op)
  }
  if (infra == FALSE) {
    results <- results[, -c(4, 13)]
  }
  if (nchar(file) > 0) {
    write.csv(results, file = file, row.names = FALSE)
  }
  return(results)
}


###################

function (sp, infra = TRUE, corr = TRUE, diffchar = 2, max.distance = 1, 
          version = "1.1", encoding = "UTF-8", author = TRUE, drop.lower.level = FALSE) 
{
  sp_orig <- sp
  sp <- gsub(paste0("\\s+|", intToUtf8(160)), " ", as.character(sp))
  sp <- gsub("(?! )\\(", " \\(", sp, perl = TRUE)
  hybrid <- ifelse(grepl(paste(paste0(" ", intToUtf8(1093)), 
                               "nothosp\\.", "grex", "[xX]", intToUtf8(215), sep = " | "), 
                         sp), TRUE, FALSE)
  if (hybrid == TRUE) {
    sp <- gsub(paste(paste0(" ", intToUtf8(1093)), "nothosp\\.", 
                     "grex", "[xX]", intToUtf8(215), sep = " | "), " ", 
               sp)
  }
  gen_hybrid <- ifelse(grepl(paste("^[xX]", intToUtf8(1093), 
                                   intToUtf8(215), sep = " |^"), sp), TRUE, FALSE)
  if (gen_hybrid == TRUE) {
    sp <- gsub(paste("^[xX]", intToUtf8(1093), intToUtf8(215), 
                     sep = " |^"), "", sp)
  }
  match.higher.level <- FALSE
  abbr <- NA
  vec <- c(" aff(\\.| )", " aggr?(\\.| |$)", "(^| )cf(\\.| )", 
           " group( |$)", ",? ?nom\\. inval\\.", ",? ?nom\\. cons\\. prop\\.", 
           ",? ?nom\\. conserv\\.", ",? ?nom\\. cons?\\.", " p\\. ?p\\.", 
           " s\\. ?l\\.?$", " s\\. ?l\\. ", " s(ens)?[\\.u] ?lat[\\.u]$", 
           "sensu latu", " s\\. ?s\\.?$", " s\\. ?s\\.? ", " s\\. ?str\\.?$", 
           " s\\. ?str\\.? ", "sensu str\\.", "sensu strictu")
  remov <- NA
  for (i in 1:length(vec)) {
    remov <- if (grepl(vec[i], sp)) {
      na.omit(c(remov, gsub(paste0(".*(", vec[i], ").*"), 
                            "\\1", sp)))
    }
    else {
      remov
    }
    sp <- gsub(vec[i], " ", sp)
  }
  remov <- unique(trimws(remov))
  rm(vec)
  Abbrev <- ifelse(any(is.na(remov)), NA, paste(remov, collapse = ","))
  vec0 <- c(" cultivar(\\.| )", " cv(\\.| )", " nothossp(\\.| )", 
            " nothosubsp(\\.| )", " nothovar(\\.| )", " subfo?(\\.| )", 
            " subvar(\\.| )", " f[ao]?(\\.| ) ?\\b(?!(var\\.|subsp\\.|f\\.|$))", 
            " fma(\\.| )", " forma?(\\.| )", " var(\\.| )", " ssp(\\.| )", 
            " subsp(\\.| )")
  if (drop.lower.level == FALSE) {
    for (j in 1:length(vec0)) {
      abbr <- ifelse(grepl(vec0[j], sp, perl = TRUE), 
                     trimws(gsub(paste0(".*(", vec0[j], ").*"), "\\1", 
                                 sp, perl = TRUE)), abbr)
      sp <- gsub(vec0[j], " ", sp, perl = TRUE)
    }
  }
  else {
    for (j in 1:length(vec0)) {
      abbr <- if (grepl(vec0[j], sp, perl = TRUE)) {
        na.omit(c(abbr, gsub(paste0(".*(", vec0[j], 
                                    ").*"), "\\1", sp, perl = TRUE)))
      }
      else {
        abbr
      }
    }
    abbr <- unique(trimws(abbr))
    if (any(grepl("^subsp\\.?$|^ssp\\.?$", abbr)) & any(grepl("^var\\.?$", 
                                                              abbr))) {
      sp <- gsub("( subsp| ssp| var)\\.(?! )", "\\1\\. ", 
                 sp, perl = TRUE)
      split <- unlist(strsplit(sp, "\\s+"))
      subsp_pos <- grep("^subsp\\.?$|^ssp\\.?$", split)
      var_pos <- grep("^var\\.?$", split)
      if (split[2] == split[subsp_pos + 1]) {
        sp <- paste(split[-(subsp_pos:var_pos)], collapse = " ")
        abbr <- abbr[-grep("^subsp\\.?$|^ssp\\.?$", 
                           abbr)]
      }
      else {
        sp <- unlist(strsplit(sp, "var\\.?"))[1]
        sp <- gsub(" subsp\\.?| ssp\\.?", "", sp)
        abbr <- abbr[-grep("^var\\.?$|^f[ao]?\\.?$|^fma\\.?$|^forma?\\.?$", 
                           abbr)]
      }
      rm(split, subsp_pos, var_pos)
    }
    if (any(grepl("^subsp\\.?$|^ssp\\.?$", abbr)) & any(grepl("^f[ao]?\\.?$|^fma\\.?$|^forma?\\.?$", 
                                                              abbr))) {
      sp <- gsub("( subsp| ssp| f[ao]?| fma| forma?)\\.(?! )", 
                 "\\1\\. ", sp, perl = TRUE)
      split <- unlist(strsplit(sp, "\\s+"))
      subsp_pos <- grep("^subsp\\.?$|^ssp\\.?$", split)
      f_pos <- grep("^f[ao]?\\.?$|^fma\\.?$|^forma?\\.?$", 
                    split)
      if (split[2] == split[subsp_pos + 1]) {
        sp <- paste(split[-(subsp_pos:f_pos)], collapse = " ")
        abbr <- abbr[-grep("^subsp\\.?$|^ssp\\.?$", 
                           abbr)]
      }
      else {
        sp <- unlist(strsplit(sp, "f[ao]?\\.?|fma\\.?|forma?\\.?"))[1]
        sp <- gsub(" subsp\\.?| ssp\\.?", "", sp)
        abbr <- abbr[-grep("^f[ao]?\\.?$|^fma\\.?$|^forma?\\.?$", 
                           abbr)]
      }
      rm(split, subsp_pos, f_pos)
    }
    if (any(grepl("^var\\.?$", abbr)) & any(grepl("^f[ao]?\\.?$|^fma\\.?$|^forma?\\.?$", 
                                                  abbr))) {
      sp <- gsub("( var| f[ao]?| fma| forma?)\\.(?! )", 
                 "\\1\\. ", sp, perl = TRUE)
      split <- unlist(strsplit(sp, "\\s+"))
      var_pos <- grep("^var\\.?$", split)
      f_pos <- grep("^f[ao]?\\.?$|^fma\\.?$|^forma?\\.?$", 
                    split)
      if (split[2] == split[var_pos + 1]) {
        sp <- paste(split[-(var_pos:f_pos)], collapse = " ")
        abbr <- abbr[-grep("^var\\.?$", abbr)]
      }
      else {
        sp <- unlist(strsplit(sp, "f[ao]?\\.?|fma\\.?|forma?\\.?"))[1]
        sp <- gsub(" var\\.", "", sp)
        abbr <- abbr[-grep("^f[ao]?\\.?$|^fma\\.?$|^forma?\\.?$", 
                           abbr)]
      }
      rm(split, var_pos, f_pos)
    }
    sp <- gsub(paste0(vec0, collapse = "|"), " ", sp, perl = TRUE)
    if (length(abbr) > 1) {
      warning(paste0(sp_orig, ": Multiple infraspecific abbreviations occur: '", 
                     paste0(abbr, collapse = "', '"), "'. Only '", 
                     abbr[length(abbr)], "' will be used."))
      abbr <- abbr[length(abbr)]
    }
  }
  sp <- trimws(gsub(" +", " ", sp))
  Rank <- abbr
  if (!is.na(abbr)) {
    abbr1 <- ifelse(grepl("^f[ao]?\\.?$|^fma\\.?$|^forma?\\.?$", 
                          abbr), "f.", abbr)
    abbr1 <- ifelse(grepl("^var\\.?$", abbr), "var.", abbr1)
    abbr1 <- ifelse(grepl("^subsp\\.?$|^ssp\\.?$", abbr), 
                    "subsp.", abbr1)
  }
  genus <- unlist(strsplit(sp, " "))[1]
  species <- unlist(strsplit(sp, " "))[2]
  if (corr == TRUE) {
    species <- tolower(species)
  }
  if (author == FALSE) {
    infrasp <- ifelse(length(unlist(strsplit(sp, " "))) > 
                        2, unlist(strsplit(sp, " "))[3], "")
    auth <- ""
    if (corr == TRUE) {
      infrasp <- tolower(infrasp)
    }
  }
  else {
    spparts <- unlist(strsplit(sp, " "))
    if (length(spparts) > 2) {
      authorstrings <- c("and", "anon\\.", "auctt?(\\.| )", 
                         "div\\.", "emend\\.", "mult\\.", "non", "plur\\.", 
                         "des?", "et", "ex\\.?", "al\\.\\)?", "f\\.\\)?", 
                         "fil\\.\\)?", "hort\\.", "sensu", "van", "von", 
                         "v\\.", "d'")
      infrasp_pos <- grep(paste0("^(?!^", paste0(authorstrings, 
                                                 collapse = "$|^"), ")[a-z]"), spparts[-c(1:2)], 
                          perl = TRUE) + 2
      if (length(infrasp_pos) > 1) {
        if (length(unique(spparts[infrasp_pos])) > 1) {
          warning(paste0(sp_orig, ": The infraspecific epithet was not unambiguously matched"))
        }
        infrasp_pos <- infrasp_pos[1]
      }
      infrasp <- ifelse(length(infrasp_pos) > 0, spparts[infrasp_pos], 
                        "")
      auth <- ifelse(nchar(infrasp) == 0, paste(spparts[3:length(spparts)], 
                                                collapse = " "), ifelse(infrasp_pos == length(spparts) && 
                                                                          length(spparts) > 3, paste(spparts[3:(length(spparts) - 
                                                                                                                  1)], collapse = " "), ifelse(length(spparts) == 
                                                                                                                                                 3, "", ifelse(infrasp_pos < length(spparts), 
                                                                                                                                                               paste(spparts[(infrasp_pos + 1):length(spparts)], 
                                                                                                                                                                     collapse = " "), ""))))
      auth_sp <- ifelse(nchar(infrasp) > 0 && infrasp_pos > 
                          3 && infrasp_pos < length(spparts), paste(spparts[3:(infrasp_pos - 
                                                                                 1)], collapse = " "), "")
      auth_infra <- ifelse(nchar(infrasp) > 0 && infrasp_pos < 
                             length(spparts), TRUE, FALSE)
      if (grepl(paste0("^(?!^", paste0(authorstrings, 
                                       collapse = "|^"), ")[a-z]"), auth, perl = TRUE)) {
        warning(paste0(sp_orig, ": The author was probably incorrectly matched."))
      }
      rm(authorstrings, infrasp_pos)
    }
    else {
      infrasp <- ""
      auth <- ""
    }
  }
  vv <- ifelse(version == "1.0", "", version)
  searchstring <- paste("http://www.theplantlist.org/tpl", 
                        vv, "/search?q=", genus, "+", species, "&csv=true", 
                        sep = "")
  table.sp <- try(read.csv(searchstring, header = TRUE, sep = ",", 
                           fill = TRUE, colClasses = "character", as.is = TRUE, 
                           encoding = encoding), silent = TRUE)
  if (class(table.sp) == "try-error") {
    stop("Cannot read TPL website.")
  }
  Genus <- genus
  Species <- species
  Infraspecific <- infrasp
  New.Hybrid.marker <- ""
  marker <- FALSE
  marker.infra <- FALSE
  Plant.Name.Index.final <- NA
  if (is.null(table.sp)) {
    Family <- ""
    Taxonomic.status <- ""
    Plant.Name.Index <- FALSE
    New.Genus <- Genus
    New.Species <- Species
    New.Infraspecific <- Infraspecific
    New.Infraspecific.rank <- Rank
    New.Authority <- ""
    New.Taxonomic.status <- ""
    Typo <- FALSE
    WFormat <- TRUE
    ID <- ""
    New.ID <- ""
  }
  else if (!is.null(table.sp)) {
    k <- dim(table.sp)[2]
    z <- dim(table.sp)[1]
    if (k == 1) {
      Family <- ""
      Plant.Name.Index <- FALSE
      Taxonomic.status <- ""
      New.Genus <- Genus
      New.Species <- Species
      New.Infraspecific <- Infraspecific
      New.Infraspecific.rank <- Rank
      New.Authority <- ""
      New.Taxonomic.status <- ""
      Typo <- FALSE
      WFormat <- FALSE
      ID <- ""
      New.ID <- ""
    }
    else if (k > 1) {
      if (z == 0) {
        Family <- ""
        Plant.Name.Index <- FALSE
        Taxonomic.status <- ""
        New.Genus <- Genus
        New.Species <- Species
        New.Infraspecific <- Infraspecific
        New.Infraspecific.rank <- Rank
        New.Authority <- ""
        New.Taxonomic.status <- ""
        Typo <- FALSE
        WFormat <- FALSE
        ID <- ""
        New.ID <- ""
      }
      else if (z > 1) {
        if (length(unique(paste(table.sp$Genus, table.sp$Species))) > 
            1 && corr == TRUE) {
          spx <- length(agrep(species, "sp", max.distance = 0)) + 
            length(agrep(species, "sp.", max.distance = 0))
          mf <- c(as.character(1:1000))
          is.mf <- agrep(species, mf, max.distance = list(deletions = 1), 
                         value = TRUE)
          cck <- agrep(species, table.sp$Species, value = TRUE, 
                       max.distance = max.distance)
          ddf <- abs(nchar(cck) - nchar(species))
          if (length(cck) > 0) {
            cck <- cck[ddf == min(ddf)]
            ddf <- abs(nchar(cck) - nchar(species))
          }
          levs <- length(unique(cck))
          if (levs > 1) {
            warning(paste("The specific epithet of", 
                          sp_orig, "could not be matched, and multiple corrections are possible."))
          }
          spstring <- ifelse(length(cck) > 0, cck[1], 
                             "")
          if (length(cck) == 0) {
            spec1 <- character()
            spec2 <- character()
            if (grepl("is$", species)) {
              spec1 <- gsub("is$", "e", species)
            }
            if (grepl("e$", species)) {
              spec1 <- gsub("e$", "is", species)
            }
            if (grepl("ii$", species)) {
              spec1 <- gsub("ii$", "i", species)
            }
            if (grepl("(?<!i)i$", species, perl = TRUE)) {
              spec1 <- gsub("i$", "ii", species)
            }
            if (grepl("us$", species)) {
              spec1 <- gsub("us$", "a", species)
              spec2 <- gsub("us$", "um", species)
            }
            if (grepl("um$", species)) {
              spec1 <- gsub("um$", "a", species)
              spec2 <- gsub("um$", "us", species)
            }
            if (grepl("a$", species)) {
              spec1 <- gsub("a$", "um", species)
              spec2 <- gsub("a$", "us", species)
            }
            spec1_match <- if (length(spec1) > 0) 
              unique(try(grep(paste0("^", spec1, "$"), 
                              table.sp$Species, value = TRUE), silent = TRUE))
            else character()
            spec2_match <- if (length(spec2) > 0) 
              unique(try(grep(paste0("^", spec2, "$"), 
                              table.sp$Species, value = TRUE), silent = TRUE))
            else character()
            if (length(spec1_match) + length(spec2_match) > 
                0) {
              if (length(spec2_match) == 0) {
                spstring <- spec1_match
              }
              else if (length(spec1_match) == 0 && length(spec2_match == 
                                                          1)) {
                spstring <- spec2_match
              }
              else {
                spstring <- spec1_match
                warning(paste("The specific epithet of", 
                              sp_orig, "could not be matched, and multiple corrections are possible."))
              }
            }
          }
          if ((length(is.mf) == 0 && length(cck) > 0 && 
               ddf <= diffchar && levs == 1 && spx == 0) | 
              (exists("spec1_match") && nchar(spstring) > 
               0)) {
            searchstring <- paste("http://www.theplantlist.org/tpl", 
                                  vv, "/search?q=", genus, "+", spstring, 
                                  "&csv=true", sep = "")
            table.sp <- try(read.csv(searchstring, header = TRUE, 
                                     sep = ",", fill = TRUE, colClasses = "character", 
                                     as.is = TRUE, encoding = encoding), silent = TRUE)
            if (class(table.sp) == "try-error") {
              stop("Cannot read TPL website.")
            }
            k <- dim(table.sp)[2]
            z <- dim(table.sp)[1]
            marker <- TRUE
          }
        }
        if (length(unique(paste(table.sp$Genus, table.sp$Species))) > 
            1) {
          Family <- ""
          Plant.Name.Index <- FALSE
          Taxonomic.status <- ""
          New.Genus <- Genus
          New.Species <- Species
          New.Infraspecific <- Infraspecific
          New.Infraspecific.rank <- Rank
          New.Authority <- ""
          New.Taxonomic.status <- ""
          Typo <- FALSE
          WFormat <- FALSE
          ID <- ""
          New.ID <- ""
        }
        if (length(unique(paste(table.sp$Genus, table.sp$Species))) == 
            1) {
          grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, 
                        value = TRUE, fixed = TRUE)
          ngrep <- nchar(grep1)
          if ((length(ngrep) == 0 || abs(ngrep - nchar(infrasp)) > 
               0) && corr == TRUE && !is.na(infrasp) && 
              nchar(infrasp) > 0) {
            cck.infra <- agrep(infrasp, table.sp$Infraspecific.epithet, 
                               value = TRUE, max.distance = max.distance)
            ddf.infra <- abs(nchar(cck.infra) - nchar(infrasp))
            if (length(cck.infra) > 0) {
              cck.infra <- cck.infra[ddf.infra == min(ddf.infra)]
              ddf.infra <- abs(nchar(cck.infra) - nchar(infrasp))
            }
            levs <- length(unique(cck.infra))
            if (length(cck.infra) > 0 && levs == 1) {
              infrasp <- unique(cck.infra)
              grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, 
                            value = TRUE, fixed = TRUE)
              ngrep <- nchar(grep1)
              marker.infra <- TRUE
            }
          }
          nominal.infra <- ifelse(species == Infraspecific, 
                                  TRUE, FALSE)
          Plant.Name.Index.final <- NA
          if (infra == TRUE && ((length(grep(infrasp, 
                                             table.sp$Infraspecific.epithet, fixed = TRUE)) > 
                                 0 && any(abs(ngrep - (nchar(infrasp))) == 
                                          0)) || (length(grep(paste0("^", Infraspecific, 
                                                                     "$"), table.sp$Infraspecific.epithet, fixed = FALSE)) == 
                                                  0 && nominal.infra == TRUE))) {
            if (length(grep(paste0("^", Infraspecific, 
                                   "$"), table.sp$Infraspecific.epithet, 
                            fixed = FALSE)) == 0 && nominal.infra == 
                TRUE) {
              table.sp <- table.sp[table.sp$Infraspecific.epithet == 
                                     "", ]
              infrasp <- Infraspecific
              match.higher.level <- ifelse(nrow(table.sp) > 
                                             0, TRUE, FALSE)
            }
            else {
              table.sp <- table.sp[table.sp$Infraspecific.epithet == 
                                     infrasp, ]
            }
          }
          else if (infra == FALSE || length(grep(infrasp, 
                                                 table.sp$Infraspecific.epithet, fixed = TRUE)) == 
                   0 || abs(ngrep - (nchar(infrasp))) > 0) {
            if (nchar(infrasp) == 0 && !any(table.sp$Infraspecific.epithet == 
                                            "") && length(grep(Species, table.sp$Infraspecific.epithet)) == 
                0) {
              table.sp <- table.sp[table.sp$Infraspecific.epithet == 
                                     "", ]
              Plant.Name.Index.final <- FALSE
            }
            else if (infra == FALSE) {
              table.sp <- table.sp[table.sp$Infraspecific.epithet == 
                                     "", ]
              table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
              match.higher.level <- TRUE
            }
            else if ((table.sp$Infraspecific.epithet == 
                      "" || is.na(table.sp$Infraspecific.epithet)) == 
                     FALSE && length(grep(Species, table.sp$Infraspecific.epithet)) > 
                     0) {
              table.sp <- table.sp[table.sp$Infraspecific.epithet == 
                                     Species, ]
              table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
            }
            else if (length(grep(infrasp, table.sp$Infraspecific.epithet, 
                                 fixed = TRUE)) == 0 && nominal.infra == 
                     FALSE) {
              table.sp <- table.sp[table.sp$Infraspecific.epithet == 
                                     "", ]
              table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
              Plant.Name.Index.final <- FALSE
              match.higher.level <- ifelse(nrow(table.sp) > 
                                             0, TRUE, FALSE)
            }
          }
          k <- dim(table.sp)[2]
          z <- dim(table.sp)[1]
          if (nchar(auth) > 0) {
            if (match.higher.level == TRUE && auth_sp != 
                "") {
              auth_match <- sapply(table.sp$Authorship, 
                                   adist, auth_sp)
              auth_current <- auth_sp
              auth_current <- ifelse(any(grepl("^ex\\.?$", 
                                               spparts)), gsub(".+ ex (.*)", "\\1", 
                                                               auth_sp), auth_current)
              auth_current <- ifelse(grepl("\\(.+\\) (.*)", 
                                           auth), gsub("\\(.+\\) (.*)", "\\1", 
                                                       auth_sp), auth_current)
            }
            else {
              auth_match <- sapply(table.sp$Authorship, 
                                   adist, auth)
              auth_current <- auth
              auth_current <- ifelse(any(grepl("^ex\\.?$", 
                                               spparts)), gsub(".+ ex (.*)", "\\1", 
                                                               auth), auth_current)
              auth_current <- ifelse(grepl("\\(.+\\) (.*)", 
                                           auth), gsub("\\(.+\\) (.*)", "\\1", 
                                                       auth), auth_current)
            }
            auth_current_match <- sapply(table.sp$Authorship, 
                                         adist, auth_current)
            if (!any(auth_current_match == 0) && any(grepl("\\(.+\\)", 
                                                           table.sp$Authorship))) {
              auth_current_match <- sapply(gsub("\\(.+\\) (.*)", 
                                                "\\1", table.sp$Authorship), adist, 
                                           auth_current)
            }
            if (!any(auth_current_match == 0) && any(grepl(" ex ", 
                                                           table.sp$Authorship))) {
              auth_current_match <- sapply(gsub(".* ex (.*)", 
                                                "\\1", table.sp$Authorship), adist, 
                                           auth_current)
            }
            if (length(auth_match) > 0) {
              bestmatch <- min(as.numeric(auth_match), 
                               na.rm = TRUE)
              bestmatch <- which(auth_match == bestmatch)
              bestmatch_auth_current <- min(as.numeric(auth_current_match), 
                                            na.rm = TRUE)
              bestmatch_auth_current <- which(auth_current_match == 
                                                bestmatch_auth_current)
              if (!any(auth_match < 2) && any(auth_current_match == 
                                              0)) {
                bestmatch <- bestmatch_auth_current
              }
              if (length(unique(names(bestmatch))) > 
                  1) {
                if (!(match.higher.level == TRUE && 
                      (!is.na(Plant.Name.Index.final) && 
                       Plant.Name.Index.final == FALSE))) {
                  warning(paste(sp_orig, ": has multiple synonyms with equally matching author names; the first author will be selected."))
                }
                bestmatch <- bestmatch[1]
              }
            }
          }
          if (z == 0) {
            Family <- ""
            Plant.Name.Index <- FALSE
            Taxonomic.status <- ""
            New.Genus <- Genus
            New.Species <- Species
            New.Infraspecific <- Infraspecific
            New.Infraspecific.rank <- Rank
            New.Authority <- ""
            New.Taxonomic.status <- ""
            Typo <- ifelse(marker == TRUE || marker.infra == 
                             TRUE, TRUE, FALSE)
            WFormat <- FALSE
            ID <- ""
            New.ID <- ""
          }
          else {
            if (nchar(auth) > 0) {
              if (gsub("\\s+", "", chartr(".", " ", 
                                          names(bestmatch))) != gsub("\\s+", "", 
                                                                     chartr(".", " ", gsub(" et ", " &", 
                                                                                           gsub(" and ", " &", auth)))) && !((nchar(infrasp) > 
                                                                                                                              0 && nchar(auth_sp) == 0 && auth_infra == 
                                                                                                                              FALSE) || (nchar(infrasp) > 0 && nchar(auth_sp) > 
                                                                                                                                         0 && infra == FALSE)) && match.higher.level == 
                  FALSE) {
                warning(paste0(sp_orig, ": The input author, '", 
                               auth, "', does not exactly match the author of the selected taxon, '", 
                               table.sp$Authorship[bestmatch], "'."))
              }
              if ((gsub("\\s+", "", chartr(".", " ", 
                                           names(bestmatch))) != gsub("\\s+", "", 
                                                                      chartr(".", " ", gsub(" et ", " &", 
                                                                                            gsub(" and ", " &", auth)))) && auth_infra == 
                   FALSE) && match.higher.level == TRUE) {
                warning(paste0(sp_orig, ": The input author, '", 
                               auth, "', does not exactly match the author of the selected taxon, '", 
                               table.sp$Authorship[bestmatch], "'."))
              }
              if ((gsub("\\s+", "", chartr(".", " ", 
                                           names(bestmatch))) != gsub("\\s+", "", 
                                                                      chartr(".", " ", gsub(" et ", " &", 
                                                                                            gsub(" and ", " &", auth_sp)))) && 
                   auth_infra == TRUE && nchar(auth_sp) > 
                   0) && match.higher.level == TRUE) {
                warning(paste0(sp_orig, ": The input author, '", 
                               auth_sp, "', does not exactly match the author of the selected taxon, '", 
                               table.sp$Authorship[bestmatch], "'."))
              }
              if ((nchar(infrasp) > 0 && auth_infra == 
                   TRUE) || (nchar(infrasp) == 0 && auth_infra == 
                             FALSE)) {
                table.sp <- table.sp[bestmatch, ]
              }
            }
            if (z > 1) {
              if (!is.na(abbr) && any(table.sp$Infraspecific.rank == 
                                      abbr1)) {
                table.sp <- table.sp[table.sp$Infraspecific.rank == 
                                       abbr1, ]
              }
              if (hybrid == TRUE && any(table.sp$Species.hybrid.marker != 
                                        "")) {
                table.sp <- table.sp[table.sp$Species.hybrid.marker != 
                                       "", ]
              }
              if (hybrid == FALSE && any(table.sp$Species.hybrid.marker == 
                                         "")) {
                table.sp <- table.sp[table.sp$Species.hybrid.marker == 
                                       "", ]
              }
            }
            if (any(table.sp$Taxonomic.status.in.TPL == 
                    "Accepted")) {
              table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL == 
                                     "Accepted", ]
              Plant.Name.Index <- TRUE
              Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
              Family <- table.sp$Family[1]
              New.Genus <- table.sp$Genus[1]
              New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
              New.Species <- table.sp$Species[1]
              if (infra == T && length(grep(infrasp, 
                                            table.sp$Infraspecific.epithet, fixed = TRUE)) > 
                  0) {
                New.Infraspecific <- table.sp$Infraspecific.epithet[1]
                New.Infraspecific.rank <- table.sp$Infraspecific.rank[1]
              }
              else if (infra == F || length(grep(infrasp, 
                                                 table.sp$Infraspecific.epithet, fixed = TRUE)) == 
                       0) {
                New.Infraspecific <- ""
                New.Infraspecific.rank <- ""
              }
              New.Authority <- table.sp$Authorship[1]
              New.Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
              Typo <- ifelse(marker == TRUE || marker.infra == 
                               TRUE, TRUE, FALSE)
              WFormat <- FALSE
              ID <- table.sp[1, 1]
              New.ID <- ID
            }
            else if (any(table.sp$Taxonomic.status.in.TPL %in% 
                         c("Synonym", "Misapplied"))) {
              if (sum(table.sp$Taxonomic.status.in.TPL == 
                      "Synonym") > 0) {
                table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL == 
                                       "Synonym", ]
              }
              else if (sum(table.sp$Taxonomic.status.in.TPL == 
                           "Misapplied") > 0) {
                table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL == 
                                       "Misapplied", ]
              }
              if (sum(table.sp$Confidence.level == "H") > 
                  0) {
                table.sp <- table.sp[table.sp$Confidence.level == 
                                       "H", ]
              }
              else if (sum(table.sp$Confidence.level == 
                           "M") > 0) {
                table.sp <- table.sp[table.sp$Confidence.level == 
                                       "M", ]
              }
              if (any(!table.sp$Nomenclatural.status.from.original.data.source %in% 
                      c("Illegitimate", "Invalid")) && any(table.sp$Nomenclatural.status.from.original.data.source != 
                                                           "")) {
                table.sp <- table.sp[!table.sp$Nomenclatural.status.from.original.data.source %in% 
                                       c("Illegitimate", "Invalid"), ]
                warning(paste(sp_orig, "has more than one valid synonym; illegitimate/invalid names were avoided."))
              }
              if (nrow(table.sp) > 1) {
                warning(paste(sp_orig, "has more than one valid synonym; the first entry was selected."))
              }
              table.sp.id <- table.sp[1, 1]
              ID <- table.sp.id
              at <- try(readLines(paste("http://www.theplantlist.org/tpl", 
                                        vv, "/record/", table.sp.id, sep = ""), 
                                  encoding = encoding))
              if (class(at) == "try-error") {
                stop("Cannot read TPL website.")
              }
              if (sum(table.sp$Taxonomic.status.in.TPL == 
                      "Synonym") > 0) {
                if (version == "1.1") {
                  az <- "<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of"
                }
                else if (version == "1.0") {
                  az <- "<p>This name is a <a href=\"/about/#synonym\">synonym</a> of"
                }
              }
              else if (sum(table.sp$Taxonomic.status.in.TPL == 
                           "Misapplied") > 0) {
                if (version == "1.1") {
                  az <- "<p>In the past this name has been <a href=\"/1.1/about/#misapplied\">erroneously used</a> to refer to"
                }
                else if (version == "1.0") {
                  az <- "<p>In the past this name has been <a href=\"/about/#misapplied\">erroneously used</a> to refer to"
                }
              }
              n <- pmatch(az, at)
              nsen <- at[n]
              nsen <- unlist(strsplit(unlist(strsplit(nsen, 
                                                      split = ">")), "<"))
              gen_index <- grep("class=\"genus\"", nsen) + 
                1
              sp_index <- grep("class=\"species\"", 
                               nsen) + 1
              if (grepl("synonym", at[n])) {
                if (version == "1.1") {
                  tpl_id <- gsub("<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of <a href=\"(.+)\"><span class=\"name\">.+", 
                                 "\\1", at[n])
                }
                else if (version == "1.0") {
                  tpl_id <- gsub("<p>This name is a <a href=\"/about/#synonym\">synonym</a> of <a href=\"(.+)\"><i class=\"genus\">.+", 
                                 "\\1", at[n])
                }
              }
              if (grepl("erroneously used", at[n])) {
                if (version == "1.1") {
                  tpl_id <- gsub("<p>In the past this name has been <a href=\"/1.1/about/#misapplied\">erroneously used</a> to refer to <a href=\"(.+)\"><span class=\"name\">.+", 
                                 "\\1", at[n])
                }
                else if (version == "1.0") {
                  tpl_id <- gsub("<p>In the past this name has been <a href=\"/about/#misapplied\">erroneously used</a> to refer to <a href=\"(.+)\"><i class=\"genus\">.+", 
                                 "\\1", at[n])
                }
              }
              if (version == "1.1") {
                searchstring <- paste("http://www.theplantlist.org/tpl", 
                                      vv, "/search?q=", tpl_id, "&csv=true", 
                                      sep = "")
              }
              else if (version == "1.0") {
                searchstring <- paste("http://www.theplantlist.org/tpl", 
                                      vv, "/search?q=", nsen[gen_index], 
                                      "+", nsen[sp_index], "&csv=true", 
                                      sep = "")
              }
              kup <- length(grep("^var\\.$", nsen)) + 
                length(grep("^subsp\\.$", nsen)) + length(grep("^f\\.$", 
                                                               nsen))
              if (infra == T && kup > 0) {
                infrasp <- nsen[grep("class=\\\"infraspe\\\"", 
                                     nsen) + 1]
                abbr_index <- grep("class=\\\"infraspr\\\"", 
                                   nsen) + 1
              }
              else if (kup == 0) {
                infrasp <- ""
              }
              if (sum(table.sp$Taxonomic.status.in.TPL == 
                      "Synonym") > 0) {
                Taxonomic.status <- "Synonym"
              }
              else if (sum(table.sp$Taxonomic.status.in.TPL == 
                           "Misapplied") > 0) {
                Taxonomic.status <- "Misapplied"
              }
              table.sp <- try(read.csv(searchstring, 
                                       header = TRUE, sep = ",", fill = TRUE, 
                                       colClasses = "character", as.is = TRUE, 
                                       encoding = encoding), silent = TRUE)
              if (class(table.sp) == "try-error") {
                stop("Cannot read TPL website.")
              }
              colnames(table.sp) <- gsub("^X.U.FEFF.ID$", 
                                         "ID", colnames(table.sp))
              if (nrow(table.sp) > 1) {
                table.sp <- table.sp[table.sp$ID == 
                                       tpl_id, ]
              }
              Plant.Name.Index <- TRUE
              Family <- table.sp$Family[1]
              New.Genus <- table.sp$Genus[1]
              New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
              New.Species <- table.sp$Species[1]
              New.Infraspecific <- table.sp$Infraspecific.epithet[1]
              New.Infraspecific.rank <- table.sp$Infraspecific.rank[1]
              New.Authority <- table.sp$Authorship[1]
              New.Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
              Typo <- ifelse(marker == TRUE || marker.infra == 
                               TRUE, TRUE, FALSE)
              WFormat <- FALSE
              New.ID <- table.sp[1, 1]
            }
            else if (any(table.sp$Taxonomic.status.in.TPL == 
                         "Unresolved")) {
              Plant.Name.Index <- TRUE
              Taxonomic.status <- "Unresolved"
              Family <- table.sp$Family[1]
              New.Genus <- table.sp$Genus[1]
              New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
              New.Species <- table.sp$Species[1]
              New.Infraspecific <- table.sp$Infraspecific.epithet[1]
              New.Infraspecific.rank <- table.sp$Infraspecific.rank[1]
              New.Authority <- table.sp$Authorship[1]
              New.Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
              Typo <- ifelse(marker == TRUE || marker.infra == 
                               TRUE, TRUE, FALSE)
              WFormat <- FALSE
              ID <- table.sp[1, 1]
              New.ID <- ID
            }
          }
        }
      }
      else if (z == 1) {
        if (nchar(infrasp) > 0 && table.sp$Infraspecific.epithet == 
            "") {
          nominal.infra <- ifelse(species == Infraspecific, 
                                  TRUE, FALSE)
          match.higher.level <- TRUE
          Plant.Name.Index.final <- ifelse(nominal.infra == 
                                             TRUE, TRUE, FALSE)
          if (nchar(auth) > 0) {
            if (gsub("\\s+", "", chartr(".", " ", table.sp$Authorship)) != 
                gsub("\\s+", "", chartr(".", " ", gsub(" et ", 
                                                       " &", gsub(" and ", " &", auth_sp))))) {
              warning(paste0(sp_orig, ": The input author, '", 
                             auth_sp, "', does not exactly match the author of the selected taxon, '", 
                             table.sp$Authorship, "'."))
            }
          }
        }
        if (species != table.sp$Species && corr == TRUE) {
          cck <- agrep(species, table.sp$Species, value = TRUE, 
                       max.distance = max.distance)
          ddf <- abs(nchar(cck) - nchar(species))
          if (length(cck) > 0 && ddf <= diffchar) {
            marker <- TRUE
          }
        }
        if (nchar(auth) > 0 && !(nchar(infrasp) > 0 && 
                                 table.sp$Infraspecific.epithet == "") && ((species != 
                                                                            table.sp$Species && corr == TRUE && length(cck) > 
                                                                            0) || (species == table.sp$Species))) {
          if (gsub("\\s+", "", chartr(".", " ", table.sp$Authorship)) != 
              gsub("\\s+", "", chartr(".", " ", gsub(" et ", 
                                                     " &", gsub(" and ", " &", auth))))) {
            warning(paste0(sp_orig, ": The input author, '", 
                           auth, "', does not exactly match the author of the selected taxon, '", 
                           table.sp$Authorship, "'."))
          }
        }
        if (is.na(table.sp$Taxonomic.status.in.TPL) || 
            (species != table.sp$Species && corr == FALSE) || 
            (species != table.sp$Species && corr == TRUE && 
             (length(cck) == 0 || (length(cck) > 0 && 
                                   ddf > diffchar)))) {
          Taxonomic.status <- ""
          Plant.Name.Index <- FALSE
          Family <- ""
          New.Genus <- Genus
          New.Species <- Species
          New.Infraspecific <- Infraspecific
          New.Infraspecific.rank <- Rank
          New.Authority <- ""
          New.Taxonomic.status <- ""
          Typo <- ifelse(marker == TRUE || marker.infra == 
                           TRUE, TRUE, FALSE)
          WFormat <- FALSE
          ID <- ""
          New.ID <- ""
        }
        else if (table.sp$Taxonomic.status.in.TPL == 
                 "Synonym" || table.sp$Taxonomic.status.in.TPL == 
                 "Misapplied") {
          table.sp.id <- table.sp[1, 1]
          ID <- table.sp.id
          at <- try(readLines(paste("http://www.theplantlist.org/tpl", 
                                    vv, "/record/", table.sp.id, sep = ""), 
                              encoding = encoding))
          if (class(at) == "try-error") {
            stop("Cannot read TPL website.")
          }
          if (sum(table.sp$Taxonomic.status.in.TPL == 
                  "Synonym") > 0) {
            if (version == "1.1") {
              az <- "<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of"
            }
            else if (version == "1.0") {
              az <- "<p>This name is a <a href=\"/about/#synonym\">synonym</a> of"
            }
          }
          else if (sum(table.sp$Taxonomic.status.in.TPL == 
                       "Misapplied") > 0) {
            if (version == "1.1") {
              az <- "<p>In the past this name has been <a href=\"/1.1/about/#misapplied\">erroneously used</a> to refer to"
            }
            else if (version == "1.0") {
              az <- "<p>In the past this name has been <a href=\"/about/#misapplied\">erroneously used</a> to refer to"
            }
          }
          n <- pmatch(az, at)
          nsen <- at[n]
          nsen <- unlist(strsplit(unlist(strsplit(nsen, 
                                                  split = ">")), "<"))
          gen_index <- grep("class=\"genus\"", nsen) + 
            1
          sp_index <- grep("class=\"species\"", nsen) + 
            1
          if (grepl("synonym", at[n])) {
            if (version == "1.1") {
              tpl_id <- gsub("<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of <a href=\"(.+)\"><span class=\"name\">.+", 
                             "\\1", at[n])
            }
            else if (version == "1.0") {
              tpl_id <- gsub("<p>This name is a <a href=\"/about/#synonym\">synonym</a> of <a href=\"(.+)\"><i class=\"genus\">.+", 
                             "\\1", at[n])
            }
          }
          if (grepl("erroneously used", at[n])) {
            if (version == "1.1") {
              tpl_id <- gsub("<p>In the past this name has been <a href=\"/1.1/about/#misapplied\">erroneously used</a> to refer to <a href=\"(.+)\"><span class=\"name\">.+", 
                             "\\1", at[n])
            }
            else if (version == "1.0") {
              tpl_id <- gsub("<p>In the past this name has been <a href=\"/about/#misapplied\">erroneously used</a> to refer to <a href=\"(.+)\"><i class=\"genus\">.+", 
                             "\\1", at[n])
            }
          }
          if (version == "1.1") {
            searchstring <- paste("http://www.theplantlist.org/tpl", 
                                  vv, "/search?q=", tpl_id, "&csv=true", 
                                  sep = "")
          }
          else if (version == "1.0") {
            searchstring <- paste("http://www.theplantlist.org/tpl", 
                                  vv, "/search?q=", nsen[gen_index], "+", 
                                  nsen[sp_index], "&csv=true", sep = "")
          }
          kup <- length(grep("^var\\.$", nsen)) + length(grep("^subsp\\.$", 
                                                              nsen)) + length(grep("^f\\.$", nsen))
          if (infra == T && kup > 0) {
            infrasp <- nsen[grep("class=\\\"infraspe\\\"", 
                                 nsen) + 1]
            abbr_index <- grep("class=\\\"infraspr\\\"", 
                               nsen) + 1
          }
          else if (kup == 0) {
            infrasp <- ""
          }
          if (sum(table.sp$Taxonomic.status.in.TPL == 
                  "Synonym") > 0) {
            Taxonomic.status <- "Synonym"
          }
          else if (sum(table.sp$Taxonomic.status.in.TPL == 
                       "Misapplied") > 0) {
            Taxonomic.status <- "Misapplied"
          }
          table.sp <- try(read.csv(searchstring, header = TRUE, 
                                   sep = ",", fill = TRUE, colClasses = "character", 
                                   as.is = TRUE, encoding = encoding), silent = TRUE)
          if (class(table.sp) == "try-error") {
            stop("Cannot read TPL website.")
          }
          colnames(table.sp) <- gsub("^X.U.FEFF.ID$", 
                                     "ID", colnames(table.sp))
          if (nrow(table.sp) > 1) {
            table.sp <- table.sp[table.sp$ID == tpl_id, 
                                 ]
          }
          Plant.Name.Index <- TRUE
          Family <- table.sp$Family[1]
          New.Genus <- table.sp$Genus[1]
          New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
          New.Species <- table.sp$Species[1]
          New.Infraspecific <- table.sp$Infraspecific.epithet[1]
          New.Infraspecific.rank <- table.sp$Infraspecific.rank[1]
          New.Authority <- table.sp$Authorship[1]
          New.Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
          Typo <- ifelse(marker == TRUE || marker.infra == 
                           TRUE, TRUE, FALSE)
          WFormat <- FALSE
          New.ID <- table.sp[1, 1]
        }
        else if (table.sp$Taxonomic.status.in.TPL == 
                 "Accepted" || table.sp$Taxonomic.status.in.TPL == 
                 "Unresolved") {
          Plant.Name.Index <- TRUE
          Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
          Family <- table.sp$Family[1]
          New.Genus <- table.sp$Genus[1]
          New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
          New.Species <- table.sp$Species[1]
          New.Infraspecific <- table.sp$Infraspecific.epithet[1]
          New.Infraspecific.rank <- table.sp$Infraspecific.rank[1]
          New.Authority <- table.sp$Authorship[1]
          New.Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
          Typo <- ifelse(marker == TRUE || marker.infra == 
                           TRUE, TRUE, FALSE)
          WFormat <- FALSE
          ID <- table.sp[1, 1]
          New.ID <- ID
        }
      }
    }
  }
  Hybrid.marker <- ifelse(hybrid == TRUE, intToUtf8(215), 
                          "")
  Plant.Name.Index <- ifelse(!is.na(Plant.Name.Index.final), 
                             Plant.Name.Index.final, Plant.Name.Index)
  results <- data.frame(Taxon = sp_orig, Genus, Hybrid.marker, 
                        Species, Abbrev = as.character(Abbrev), Infraspecific.rank = as.character(Rank), 
                        Infraspecific, Authority = auth, ID, Plant.Name.Index, 
                        TPL.version = version, Taxonomic.status, Family, New.Genus, 
                        New.Hybrid.marker, New.Species, New.Infraspecific.rank, 
                        New.Infraspecific, New.Authority, New.ID, New.Taxonomic.status, 
                        Typo, WFormat, Higher.level = match.higher.level, Date = Sys.Date(), 
                        stringsAsFactors = FALSE)
}
