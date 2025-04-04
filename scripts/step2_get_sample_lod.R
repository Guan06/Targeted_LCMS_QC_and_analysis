
library(ggplot2)
library(dplyr)
 
dat <- read.table("./s0_reformat_data.txt", header = T, sep = "\t")
lod <- read.table("s0_std_LOD.txt", header = T, sep = "\t")

## Detect if there is different volumn for standards amd samples
## If yes, add column Concentration_volumn_adjust

nv <- length(unique(dat$Vol.))
if (nv > 1) {
    dat$Concentration_volumn_adjust <- dat$Calculated_Concentration/dat$Vol.
}

## read in standard file and get the average area of lower LOD
std <- read.table("s0_reformat_std_concentration.txt", header = T, sep = "\t")
std <- std[, c("Compound", "Level", "Concentration")]

lod_std <- merge(lod, std)
lod_std <- lod_std[lod_std$LOD == lod_std$Concentration, ]

## Get the area of standards
dat_std <- dat[dat$Type == "Cal", ]
dat_std <- dat_std[, c("Compound", "Area", "Level")]

dat_std <- dat_std %>% group_by(Compound, Level) %>%
    summarise(Area_mean = mean(Area))

lod_std_dat <- merge(dat_std, lod_std, by = c("Compound", "Level"))
lod <- lod_std_dat[, c("Compound", "LOD", "Area_mean")]

## give 1e-5 to Blanck samples and 1e+03 to samples
dat1 <- dat[dat$Type == "Cal", ]

## get the Area of the standard with highest concentration
lod2 <- dat1[dat1$Level == 8, ]
lod2 <- as.data.frame(lod2 %>% group_by(Compound) %>%
                      summarise(Max_Area = mean(Area)))

dat3 <- dat[dat$Type == "Sample", ]

dat3 <- merge(dat3, lod)
dat3 <- merge(dat3, lod2)

dat3$Within_detection <- ifelse(dat3$Area < dat3$Area_mean, "Lower",
                ifelse(dat3$Area > dat3$Max_Area, "Upper", "Yes"))

### normalize by internal standards
# Normalize the concentration by internal standard
note <- read.csv("notes.csv", header = T)
skip_lst <- note[note$Annotation == "Exclude", ]$Compound
internal_lst <- note[note$InternalStandard == "x", ]$Compound

## Get the mean Area for the internal standards -> IN SAMPLES [2023.11.17]
## Get the CV value (sd / mean) and chose the smallest one

## CAREFUL: take into account of different volumn of samples too!
#### 2023.11.17 NOTES: in QC samples, no internal standards were added!

inter_df <- c()
n <- length(internal_lst)

if (n > 1) {
    for (i in 1:n){
        this_internal <- internal_lst[i]
        this_internal_dat <- dat[dat$Compound == this_internal, ]
        
        ## only keep Samples
        this_internal_dat <- this_internal_dat[this_internal_dat$Type == "Sample", ]
        ## filter out QC samples (which do NOT contains internal standards)
        no_qc <- !sapply(this_internal_dat$Name, function(x){ grepl("QC", x)})
        this_internal_dat <- this_internal_dat[no_qc, ]

        this_df <- this_internal_dat$Area
        names(this_df) <- this_internal_dat$Name

        inter_df <- rbind(inter_df, this_df)
        rownames(inter_df)[i] <- this_internal
    }

	library(Hmisc)
	inter_cor <- rcorr(t(inter_df))
	cor_mat <- inter_cor$r
	p_mat <- inter_cor$P
	write.table(cor_mat, "s2_internal_cor.txt", quote = F, sep = "\t")
	write.table(p_mat, "s2_internal_cor_P.txt", quote = F, sep = "\t")

	cv <- apply(inter_df, 1, function(x) sd(x)/mean(x))
	print(cv)

    internal_compound <- names(which.min(cv))
} else {
    internal_compound <- internal_lst
}

## According to suggestion from Stephan, we will use Amoxicilin for normalization
internal_compound <- "Amoxicilin"

print(paste("The internal compound used for normalization is", internal_compound))

internal <- dat[(dat$Compound == internal_compound) & (dat$Type == "Sample"), ]
#internal$Factor <- internal$Area / mean(internal$Area)

internal <- as.data.frame(internal %>% group_by(Vol.) %>%
                            mutate(Factor = Area / mean(Area)))

re_normalize <- internal[, c("Name", "Factor")]

dat3 <- merge(dat3, re_normalize, all = T)

## CAREFUL: take into account of different volumn of samples too!
if (nv > 1) {
    dat3$Calculated_concentration_reNorm <- dat3$Concentration_volumn_adjust / dat3$Factor
} else {
    dat3$Calculated_concentration_reNorm <- dat3$Calculated_Concentration / dat3$Factor
}

### normalize samples by pqn, use the first sample as the "reference" sample
mat <- read.table("s0_reformat_data_matrix.txt", header = T, sep = "\t")
#
### change the choise of reference sample in the future
ref <- mat[3, ]
print("Reference sample for pqn normalization: ")
print(rownames(mat)[3])
coe <- apply(mat, 1, function(x) median(as.numeric(x/ref)))
#
reNorm_df <- data.frame(Name = names(coe), Factor_pqn = coe)
#

## CAREFUL: take into account of different volumn of samples too!
dat4 <- merge(dat3, reNorm_df)
if (nv > 1) {
    dat4$Calculated_concentration_reNorm_pqn <- dat4$Concentration_volumn_adjust / dat4$Factor_pqn
} else {
    dat4$Calculated_concentration_reNorm_pqn <- dat4$Calculated_Concentration / dat4$Factor_pqn
}

write.table(dat4, "s2_reformat_data_with_LOD_reNorm.txt",
            quote = F, sep = "\t", row.names = F)

