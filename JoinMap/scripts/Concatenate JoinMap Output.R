# IWG_Nested_Association_Mapping_Introduction
# Compile Joinmap Output
# Kayla Altendorf 
# 6/16/2020

# load packages
library("dplyr")
library("stringr")
library("ggpubr")
library("ggplot2")
library("tidyr")

#### Step 1: Concatenate the Loc Files ####
# these are the same files that were used to create the maps on 
# an indiviual chromosome basis in JoinMap
# now, since we're going to mapQTL on a per family basis, it makes sense to contactenate them
# per family. we'll import them into mapQTL and markers in common with the map will be used.

# create vectors for iteration
chroms_one_digit <- 1:21
chroms_two_digit <- str_pad(chroms_one_digit, width = 2, side = c("left"), pad = "0")
famID <- c("WGN07", "WGN15", "WGN26", "WGN36", "WGN38", "WGN39", "WGN45", "WGN46", "WGN55", "WGN63")
type <- c("u", "d") # two types - no distorted (u) or distorted markers (d) allowed
round <- 1:3 # rounds of mapping

# create directories for output
dir <- "/users/kraltendorf/desktop/JoinMap Files/concat_files"
dir.create(dir)
dir.create(paste(dir, "/u", sep = "")) # for undistorted
dir.create(paste(dir, "/d", sep = "")) # for distorted

# create lists for output
loc_file <- list()
loc_file_chr <- list()

# read in files, concatenate them, and write output on a per family basis
for (k in 1:length(type)) {
  for (j in 1:length(famID)) {
    for (i in 1:length(chroms_one_digit)) {
      loc <- read.table(paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", type[k], famID[j], "_CON", chroms_two_digit[i], ".loc", sep = ""), skip = 7, fill = TRUE)
      line <- which(grepl("individual", loc$V1)) # find where the individual names start
      loc1 <- slice(loc, 1:line-1) # remove them
      names <- slice(loc, (line+1):nrow(loc)) %>% dplyr::select(V1, V2) # remove comments
      loc_file[[i]] <- loc1
  }
    loc_file_chr <- do.call("rbind", loc_file)
    file <- paste(dir, "/", type[k], "/", famID[j], ".loc", sep = "")
    write.table(paste("name = ", famID[j], sep = "" ), file = file, row.names = F, col.names = F, quote = F, sep = "\t") # add the necessary file components again
    write.table(paste("popt = CP"), file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    write.table(paste("nloc = ", length(which(grepl("Chr", loc_file_chr$V1))), sep = ""), file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    write.table(paste("nind = ", length(names$V1), sep = ""), file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    cat("\n", file = file, append=TRUE)
    write.table(loc_file_chr, file = file, quote = F, row.names = F, col.names = F, append = T)
    cat("\n", file = file, append=TRUE)
    write.table("individual names:", file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    write.table(names, file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  }
}
    

#### Step 2: Concatenate the Map Files ####
# this is the output after creating the maps in JoinMap
# and exporting them individually by chromosome
# again, since we'll be mapping all chromosomes at once,
# it makes sense to contactenate them into one map
# there were two approaches -- one allowing distorted markers ("d" for distorted), 
# and the other not ("u" for undistorted)

# create directories for output
dir <- "/users/kraltendorf/desktop/JoinMap Files/concat_files"

# create list for output
map_file <- list()

# read in map files, concatenate
# if R3 does not exist, default to R2 
# in some cases we weren't able to create a R3 map due to computational time
# if marker count was high

for (k in 1:length(round)) {
  for (i in 1:length(chroms_one_digit)) {
    file <- paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", "CON", chroms_two_digit[i], "_R", round[k], ".map", sep = "")
    if (file.exists(file)) {
      map <- read.table(file, skip = 5, fill = TRUE)
      map_file[[i]] <- map
    }
    else { # if that round does not exist, load in the previous one
      file <- paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", "CON", chroms_two_digit[i], "_R", round[k-1], ".map", sep = "")
      map <- read.table(file, skip = 5, fill = TRUE)
      map_file[[i]] <- map
    }
    file_out <- paste(dir, "/", "R", round[k], ".map", sep = "")
    write.table(paste("group ", chroms_one_digit[i], sep = ""), file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    cat("\n", file = file_out, append = T)
    write.table(map_file[[i]], file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    cat("\n", file = file_out, append = T)
  }
}


# write this out again, but do not add "group" lines. 
# this will make it easier to use in R

for (k in 1:length(round)) {
  for (i in 1:length(chroms_one_digit)) {
    file <- paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", "CON", chroms_two_digit[i], "_R", round[k], ".map", sep = "")
    if (file.exists(file)) {
      map <- read.table(file, skip = 5, fill = TRUE)
      map_file[[i]] <- map
    }
    else { # if that round does not exist, load in the previous one
      file <- paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", "CON", chroms_two_digit[i], "_R", round[k-1], ".map", sep = "")
      map <- read.table(file, skip = 5, fill = TRUE)
    write.table(map_file[[i]], file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  }
  }
}


#### Step 3: Create Inverted Map Files ####
# some maps were inverted during the creation
# too keep them consistent with their physical positions and
# the previous IWG maps, we'll invert them here

# with map designator for use in mapqtl - these are the linkage groups that require inversion
invert <- c(01, 03, 05, 07, 08, 09, 10, 13, 14, 16, 17, 20)

for (k in 1:length(round)) {
  for (i in 1:length(chroms_one_digit)) {
    file <- paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", "CON", chroms_two_digit[i], "_R", round[k], ".map", sep = "")
    if (file.exists(file)) {
      map <- read.table(file, skip = 5, fill = TRUE)
      if (i %in% invert) { # if map number is in invert list, invert it here and rename "map"
       
        # here we take the maximum position on the map, and take the absolute value of every postion minus the max to reverse it
        # then arrange so it's in order
        map <- map %>% mutate(V2 = round(abs(V2 - max(map$V2)), 3)) %>% arrange(V2)
      }
      map_file[[i]] <- map
    }
    else { # if that round does not exist, load in the previous one
      file <- paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", "CON", chroms_two_digit[i], "_R", round[k-1], ".map", sep = "")
      map <- read.table(file, skip = 5, fill = TRUE)
      if (i %in% invert) { # if map number is in invert list, invert it here and rename "map"
        map <- map %>% mutate(V2 = round(abs(V2 - max(map$V2)), 3)) %>% arrange(V2)
      }
      map_file[[i]] <- map
    }
    file_out <- paste(dir, "/", "R", round[k], "_inverted.map", sep = "") # change file name to inverted
    write.table(paste("group ", chroms_one_digit[i], sep = ""), file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    cat("\n", file = file_out, append = T)
    write.table(map_file[[i]], file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    cat("\n", file = file_out, append = T)
  }
}

# write this out again, but do not add "group" lines. 
# this will make it easier to use in R

for (k in 1:length(round)) {
  for (i in 1:length(chroms_one_digit)) {
    file <- paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", "CON", chroms_two_digit[i], "_R", round[k], ".map", sep = "")
    if (file.exists(file)) {
      map <- read.table(file, skip = 5, fill = TRUE)
      if (i %in% invert) { # if map number is in invert list, invert it here and rename "map"
        map <- map %>% mutate(V2 = round(abs(V2 - max(map$V2)), 3)) %>% arrange(V2)
        }
      map_file[[i]] <- map
    }
    else { # if that round does not exist, load in the previous one
      file <- paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", "CON", chroms_two_digit[i], "_R", round[k-1], ".map", sep = "")
      map <- read.table(file, skip = 5, fill = TRUE)
      map_file[[i]] <- map
      if (i %in% invert) { # if map number is in invert list, invert it here and rename "map"
        map <- map %>% mutate(V2 = round(abs(V2 - max(map$V2)), 3)) %>% arrange(V2)
        }
      map_file[[i]] <- map
    }
  }
  file_out <- paste(dir, "/", "R", round[k], "_no_group_inverted.map", sep = "") # change file name to inverted
  write.table(do.call("rbind", map_file), file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = F)
}




#### Step 4: Compare with the Existing Consensus Map ####
# using the previously created IWG consensus map:
# Citation: Kantarski, T., S. Larson, X. Zhang, L. DeHaan, J. Borevitz, J. Anderson, and J. Poland. 2016. 
# Development of the first consensus genetic map of intermediate wheatgrass (Thinopyrum intermedium) 
# using genotyping-by-sequencing. Theor. Appl. Genet.: 1–14.

# read in consensus map
con <- read.csv("/Users/Kraltendorf/Desktop/Flowering Time/Scripts for Github/ConsensusMap.csv", header = T)
con <- con %>% dplyr::select(GBSV2RS, cM) %>% dplyr::rename(CON_cM = cM) # select necessary columns

# read in my map 
my_map <- read.table("/users/kraltendorf/desktop/JoinMap Files/concat_files/R2_no_group_inverted.map", header = F) # round 2

# format a bit
my_map1 <- my_map %>% mutate(GBSV2_Chrom = substring(V1, 4, 5)) %>%
  filter(V2 >= 0.00) %>% # remove negative positions
  mutate(GBSV2RS = paste("S", GBSV2_Chrom, "_", substr(V1, 7, nchar(as.character(V1))), sep = "")) %>%
  mutate(linkage_group = as.factor(as.integer(GBSV2_Chrom))) %>%
  dplyr::select(-GBSV2_Chrom, -V1) %>%
  dplyr::rename(NAM_cM = V2)

# find markers in common
common <- inner_join(my_map1, con, by = "GBSV2RS")

# plot correlations of both maps
# set the panel order
panel.order<- c("1","4","7","10","13","16","19","2","5","8","11",'14','17',"20","3","6","9","12","15","18","21")

ggplot(common, aes(CON_cM, NAM_cM)) + 
  geom_smooth(method = "lm", color = "#0072B2", size = 1) + 
  geom_point(size = 1.5, color = "#616365") + 
  stat_cor(method = "pearson", label.x.npc = c("left"), label.y.npc = c("top"), size = 4) +
  facet_wrap(~factor(linkage_group, levels = panel.order), nrow = 3, ncol = 7) +
  theme_bw() +
  xlab("Consensus Map Position (cM)") +
  ylab("NAM Map Position (cM)")+ 
  theme( 
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 22),
    axis.title.x = element_text(size = 22, color = "black"),
    strip.text.x = element_text(size = 20),
    strip.background.x = element_rect(fill = "white", colour = NA), 
    strip.text.y = element_text(size = 25),
    strip.background.y = element_rect(fill = "white", colour = NA), 
  )

# how many markers were in common?
length(common$NAM_cM) # 1159 markers
1159/21 # would have been on average 55 markers per chrom if I used existing map 
# but that doesn't mean all of them were "informative" in every population. 


#### Step 5: Preparing Files for the DH Model a.k.a. the Two-Way Psuedo Testcross Model ####

# make new directory for these files
dir <- c("/Users/Kraltendorf/Desktop/JoinMap Files/dh_files/")
dir.create(dir)
setwd(dir)

# create lists for output
nnnp_markers <- list()
lmll_markers <- list()
hkhk_markers <- list()

nnnp_markers_all <- list()
lmll_markers_all <- list()
hkhk_markers_all <- list()

k =1 
j = 1
# read in loc files, get a list of markers to keep, and write output on a per family basis
for (k in 1:length(type)) {
  for (j in 1:length(famID)) {
    loc <- read.table(paste("/Users/Kraltendorf/Desktop/JoinMap Files/concat_files/", type[k], "/", famID[j], ".loc", sep = ""), skip = 5, fill = TRUE)
    for (q in 1:nrow(loc)) {
      if (loc$V2[q] == "<lmxll>") { # if marker type is lmxll
        lmll <- slice(loc, q) %>% dplyr::select(V1) # slice out that row, select marker name, or V1
        lmll_markers[[q]] <- lmll # add to the list of lmll_markers
      }
      else if (loc$V2[q] == "<nnxnp>") { # now the same for the following marker types
        nnnp <- slice(loc, q) %>% dplyr::select(V1)
        nnnp_markers[[q]] <- nnnp
      }
      else if (loc$V2[q] == "<hkxhk>") {
        hkhk <- slice(loc, q) %>% dplyr::select(V1)
        hkhk_markers[[q]] <- hkhk
      }
    }
      
  # here is the list of loci names to keep in a list by chromosome
    lmll_markers_all[[j]] <- do.call("rbind", lmll_markers) # rbind all the marker names that are that type
    nnnp_markers_all[[j]] <- do.call("rbind", nnnp_markers)
    hkhk_markers_all[[j]] <- do.call("rbind", hkhk_markers)
    
  }
  lmll_all_families <- do.call("rbind", lmll_markers_all) %>% distinct() # combine marker types across families, and keep distinct ones
  nnnp_all_families <- do.call("rbind", nnnp_markers_all) %>% distinct()
  hkhk_all_families <- do.call("rbind", hkhk_markers_all) %>% distinct()
}

# now that we know which ones to filter/keep, proceed
# edit map files only to include lm or np markers

# create list for output
exclude_lm_by_chrom <- list()
exclude_np_by_chrom <- list()

# load in map files by round, if round 3 is not there, default to 2
for (i in 1:length(round)) {
  file <- paste("/Users/Kraltendorf/Desktop/JoinMap Files/concat_files/R", round[i], "_no_group.map", sep = "")
  if (file.exists(file)) {
      map <- read.table(file)
      }
  else { # if that round does not exist, load in the previous one
    file <- paste("/Users/Kraltendorf/Desktop/JoinMap Files/concat_files/R", round[i-1], "_no_group.map", sep = "")
    map <- read.table(file)
  }
  
  # filter 
  DP <- map %>% filter(! V1 %in% lmll_all_families$V1) # exclude all lms, leaving only markers that are informative in the donor parent
  CP <- map %>% filter(! V1 %in% nnnp_all_families$V1) # exclude all nps, leaving only markers that are informative in the common parent
  
  # create a chromosome column to filter on
  DP_chrom <- DP %>% mutate(chrom = substr(V1, 4,5)) %>% mutate(V1 = paste("DP_", V1, sep = ""))
  CP_chrom <- CP %>% mutate(chrom = substr(V1, 4,5)) %>% mutate(V1 = paste("CP_", V1, sep = ""))
  
  # begin writing map file
  file_out <- paste(dir, "DH_R" , round[i], ".map", sep = "")
  write.table(paste("; DH round ", round[i], sep = ""), file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = F)
  
  # for creating the DONOR PARENT MAP
  # DP iterate through chromosomes, and write out files while removing added chrom col
  for (j in 1:length(chroms_two_digit)) {
    cat("\n", file = file_out, append = T)
    write.table(paste("group ", chroms_one_digit[j], "_DP",  sep = ""), file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    cat("\n", file = file_out, append = T)
    group <- DP_chrom %>% filter(chrom == chroms_two_digit[j]) %>% dplyr::select(-chrom)
    write.table(group, file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    }
  
  # for creating the COMMON PARENT MAP
  # CP iterate through chromosomes, and write out files while removing added chrom col
  for (j in 1:length(chroms_two_digit)) {
    cat("\n", file = file_out, append = T)
    write.table(paste("group ", chroms_one_digit[j], "_CP", sep = ""), file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    cat("\n", file = file_out, append = T)
    group <- CP_chrom %>% filter(chrom == chroms_two_digit[j])  %>% dplyr::select(-chrom)
    write.table(group, file = file_out, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  }
}


#### Step 6: Edit .loc Files for DH Approach ####
# NOTE: family 38 has too few progeny and has to be prepared separately because rownumbers become off
# change line 2:   for (j in 5:5) #length(famID)) {
# change this from 6 to 5: for (q in 2:5) {
# change this from 5 to 4: locus <- slice(loc2, p:(p+4)) 

# create lists for output
loc_file <- list()

# read in files, concatenate them, and write output on a per family basis
output <- list()

# make output directories for u and d (undistorted and distorted)
dir.create(paste(dir, "u", sep = ""))
dir.create(paste(dir, "d", sep = ""))

# convert genotypes to the DH model
for (k in 1:length(type)) {
  for (j in 1:length(famID)) {
    for (i in 1:length(chroms_one_digit)) {
      
      # import and process the loc file by removing header, individuals, and loci 
      loc <- read.table(paste("/users/kraltendorf/desktop/JoinMap Files/chr", chroms_one_digit[i], "/", type[k], famID[j], "_CON", chroms_two_digit[i], ".loc", sep = ""), skip = 7, fill = TRUE)
      line <- which(grepl("individual", loc$V1)) # find where the individual names start
      loc1 <- slice(loc, 1:line-1) # remove them
      names <- slice(loc, (line+1):nrow(loc)) %>% dplyr::select(V1, V2) # remove comments
    
      # iterate through the loc file and edit the genotypes
      loc2 <- as.data.frame(loc1, stringsAsFactors = FALSE)

      for (p in 1:nrow(loc2)) {
        if (loc2$V2[p] %in% c("<lmxll>", "<nnxnp>")) { # if the marker is either lmxll or nnxnp
          locus <- slice(loc2, p:(p+5)) # slice that locus out, plus the following rows of progeny genotypes
          locus_ed <- data.frame(lapply(locus, as.character), stringsAsFactors=FALSE) # make it a dataframe
            for (q in 2:6) { # for each row 
            call <-  unlist(locus_ed[q,]) # unlist the genotypes
            call1 <- call
            call1[call == "ll"] <- "A" # if ll, change to A
            call1[call == "lm"] <- "B" # if lm, change to B
            call1[call == "nn"] <- "B" # etc.
            call1[call == "np"] <- "A"
            call1[call == "hk"] <- "-" # change hkhk markers to missing, missing now is one dash instead of two
            call1[call == "kk"] <- "-"
            call1[call == "hh"] <- "-"
            call1[call == "--"] <- "-"
            locus_ed[q,] <- call1 # add the call back to the locus dataframe 
          }
        
          # extract header to determine segregation type
          first_row <- locus_ed[1,]
          
          # edit the first row so it reflects this marker type, CP or DP
          if (first_row$V2 == "<lmxll>") {
            first_row_ed <- first_row %>% mutate(V1 = paste("CP_", V1, sep = ""))  %>% # if it's a lm marker, edit marker name to CP
              mutate(V2 = V3) %>%
              mutate(V3 = "") %>% # make these columns empty
              mutate(V4 = "") %>%
              mutate(V5 = "") %>%
              mutate(V6 = "") %>%
              mutate(V2 = paste(substr(V2, 1, 2), "}", sep = "")) } 
          if (first_row$V2 == "<nnxnp>") {
            first_row_ed <- first_row %>% mutate(V1 = paste("DP_", V1, sep = "")) %>% # if its a np marker, edit marer name to DP
              mutate(V2 = V3) %>%
              mutate(V3 = "") %>%
              mutate(V4 = "") %>%
              mutate(V5 = "") %>%
              mutate(V6 = "") %>%
              mutate(V2 = paste("{", substr(V2, 3, 4), sep = ""))
          }
          locus_ed[1,] <- first_row_ed
          output[[p]] <- locus_ed # send locus to output
        }
      }
    loc2 <- do.call("rbind", output)
    
    output <- list() # empty this list after it's complete to avoid problems on next iteration
    loc_file[[i]] <- loc2 # per chromosome
    
    # provide an update to the console for monitoring progress
    print(paste(Sys.time(), "type", type[k], famID[j], "complete for chromosome ", chroms_one_digit[i]), sep = "")
    }
    
  # rbind all chromosomes together
  loc_file_chr <- do.call("rbind", loc_file)
  loc_file <- list() # empty this list after its complete
  file <- paste(dir, type[k], "/", famID[j], "_dh.loc", sep = "")
  write.table(paste("name = ", famID[j], sep = "" ), file = file, row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(paste("popt = DH"), file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  write.table(paste("nloc = ", length(which(grepl("Chr", loc_file_chr$V1))), sep = ""), file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  write.table(paste("nind = ", length(names$V1), sep = ""), file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  cat("\n", file = file, append=TRUE)
  write.table(loc_file_chr, file = file, quote = F, row.names = F, col.names = F, append = T)
  cat("\n", file = file, append=TRUE)
  write.table("individual names:", file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  write.table(names, file = file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  print(paste(Sys.time(), "type", type[k], famID[j], "all chromosomes complete" ), sep = "")
  }
}
