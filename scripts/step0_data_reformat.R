#!/opt/homebrew/bin/Rscipt

library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
lst_file <- args[1]
std_file <- args[2]
dat_file <- args[3]

note <- read.csv("notes.csv", header = T)
skip_lst <- note[note$Annotation == "Exclude", ]$Compound
internal_lst <- note[note$InternalStandard == "x", ]$Compound
###############################################################################

lst <- read.table(lst_file, header = F, sep = "\t")
std <- read.table(std_file, header = T, sep = ",", check.names = F)

## add the mapping information for the colums in dat
lst$Order <- as.numeric(rownames(lst))
lst$Order <- as.character(lst$Order - 1)

colnames(lst)[1] <- "Compound"

lst$RT <- paste0("RT.", lst$Order)
lst$Area <- paste0("Area.", lst$Order)
lst$CalCon <- paste0("Calc..Conc..", lst$Order)
lst$FWHM <- paste0("FWHM.", lst$Order)
lst$SN <- paste0("S.N.", lst$Order)

lst$RT[1] <- "RT"
lst$Area[1] <- "Area"
lst$CalCon[1] <- "Calc..Conc."
lst$FWHM[1] <- "FWHM"
lst$SN[1] <- "S.N" 

write.table(lst, "s0_reformat_compound.lst", quote = F, row.names = F, sep = "\t")

## pivot the std table
colnames(std)[1:2] <- c("Compound", "Concentration_in_stock_uM")

std <- std %>% pivot_longer(!c("Compound", "Concentration_in_stock_uM"),
                            values_to = "Concentration", names_to = "Level")

write.table(std, "s0_reformat_std_concentration.txt", quote = F, row.names = F, sep = "\t")
std <- std[, c("Compound", "Level", "Concentration")]

## read in data file and reformat
## Also detect LOD for standard samples

dat <- read.table(dat_file, header = T, sep = ",")

## number of column before "Measurement" and number of measurements per compound

### !!! In this specific data file, a column named "Data File" is added before "Name",  therefore, n_meta is 7 instead of 6. !!!###
n_meta <- 7
n_step <- 4

i <- n_meta + 1

new_dat <- c()
new_dat_mat <- dat[, c("Name", "Type")]
new_internal <- c()
lod <- c()

get_rsq <- function(x, y) {
    ## x as actual [concentration], y as preds [area]
    tss <- sum((x - mean(x)) ^ 2) ## total sum of squares
    regss <- sum((y - mean(y)) ^ 2) ## regression sum of squares
    regss / tss
}

while (i < ncol(dat)) {
    this_dat <- dat[, c(1:n_meta, i:(i + n_step))]

    # get the compound name of this compound
    ## CAREFUL: Here the order of result of each compound matters!
    ## For now, the order is: FWHM, Area, RT, S/N, Calc.Conc 
    ## Update on 2023.07.11
    ## Due to the changes in order, here cmp_id will be n_meta + 2
    cmp_id <- colnames(this_dat)[n_meta + 2]
    cmp <- lst[lst$Area == cmp_id, ]$Compound

    if (cmp %in% skip_lst){
        i <- i + n_step + 1
        next
    }

    colnames(this_dat)[(n_meta+1) : ncol(this_dat)] <- c("FWHM", "Area", "RT", "SN_ratio", "Calculated_Concentration")
    
    this_dat$Compound <- cmp
    this_dat$Time <- as.numeric(rownames(this_dat))

    this_dat$Level[is.na(this_dat$Level)] <- 0
    this_dat$Level <- as.character(this_dat$Level)
    this_dat[is.na(this_dat)] <- 0
    
    print(i)
    print(cmp)    

    new_dat <- rbind(new_dat, this_dat)
   
    if (cmp %in% internal_lst) {
        i <- i + n_step + 1
        next
    }

    ## take the Calculated_Concentration column and form a new matrix for pqn normalization.
    new_dat_mat <- cbind(new_dat_mat, this_dat$Area)
    colnames(new_dat_mat)[ncol(new_dat_mat)] <- cmp

    this_dat1 <- this_dat[this_dat$Type == "Cal", ]
    this_dat_std <- unique(merge(this_dat1, std, by = c("Compound", "Level")))
   
    preds <- unlist(this_dat_std$Area)
    
    #print(cmp)
    #print(preds)

    actual <- unlist(this_dat_std$Concentration)

    preds[is.na(preds)] <- 0
    preds <- log(preds + 1)

    initial_rsq <- get_rsq(actual, preds)
    len <- length(preds)
    # using the last 20% data point for end rsq calculation
    len2 <- round(len/5)
    end_rsq <- get_rsq(actual[(len - len2) : len], preds[(len - len2) : len])

	## get the data point where the rsq changed more than 50%
	for (j in 1:(len - 3)) {
	    if (preds[j] == 0 ) next
	
	    this_p <- preds[j : length(preds)]
	    this_a <- actual[j : length(actual)]
	
	    this_rsq <- get_rsq(this_a, this_p)
	
#	    print(paste(this_rsq, initial_rsq,end_rsq))
	
	    if (abs(this_rsq - initial_rsq) / abs(initial_rsq - end_rsq) > 0.5) {
	        lod <- rbind(lod, c(cmp, actual[j]))
	        break
	    }
	}

    i <- i + n_step + 1
}

new_dat_mat <- new_dat_mat[new_dat_mat$Type == "Sample", ]
rownames(new_dat_mat) <- new_dat_mat$Name
new_dat_mat <- new_dat_mat[, -c(1:2)]

## keep only "samples" in the matrix
#sample_lst <- dat[dat$Type == "Sample"]$Name
#new_dat_mat <- new_dat_mat[rownames(new_dat_mat) %in% sample_lst, ]
new_dat_mat <- as.matrix(new_dat_mat)

write.table(new_dat_mat, "s0_reformat_data_matrix.txt", quote = F, sep = "\t")

write.table(new_dat, "s0_reformat_data.txt", quote = F, sep = "\t", row.names = F)

colnames(lod) <-c("Compound", "LOD")
write.table(lod, "s0_std_LOD.txt", quote = F, sep = "\t", row.names = F)
