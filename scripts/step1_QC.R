
source("~/Desktop/Cam_2022/20220301_ATP_ADP_data_processing/plot_settings.R")
library(ggplot2)
library(dplyr)

std <- read.table("s0_reformat_std_concentration.txt", header = T, sep = "\t")
dat <- read.table("s0_reformat_data.txt", header = T, sep = "\t")
lst <- read.table("s0_reformat_compound.lst", header = T, sep = "\t")
lod <- read.table("s0_std_LOD.txt", header = T, sep = "\t")

## give 1e-5 to Blanck samples and 1e+03 to samples
dat1 <- dat[dat$Type == "Cal", ]
std <- std[, c("Compound", "Level", "Concentration")]

dat_std <-  merge(dat1, std, by = c("Compound", "Level"))
#dat_std$SN_bi <- ifelse(dat_std$SN_ratio > 10, "Yes", "No")
#dat_std <- merge(dat_std, lod)

dat2 <- dat[dat$Type == "Blank", ]
dat2$Concentration <- 0.0001

dat3 <- dat[dat$Type == "Sample", ]
dat3$Concentration <- 1000

dat_std <- rbind(dat_std, dat2)
dat_std <- rbind(dat_std, dat3)

dat_std <- merge(dat_std, lod)
dat_std$SN_bi <- ifelse(dat_std$SN_ratio > 10, "Yes", "No")

p <-  ggplot(dat_std, aes(x = Concentration, y = Area)) +
        facet_wrap(~ Compound, scales = "free") +
        geom_vline(aes(xintercept = LOD), linetype = "dotted",
                   color = "salmon", linewidth = 0.6) +
        geom_point(aes(shape = SN_bi, color = Type)) +
        scale_shape_manual(values = c(1, 2)) +
        scale_x_continuous(trans = 'log10') +
        scale_y_continuous(trans = 'log10') +
        main_theme +
        theme_minimal() +
        labs(x = "Concentration", y = "Area") +
        theme(legend.position = "top") +
        scale_color_manual(values = c("plum", "#F0CD42", "#56B4E9"))

ggsave("s1_std_curve.pdf", width = 16, height = 16)

dat_std <- dat_std %>% group_by(Compound) %>% 
    mutate(Mean = mean(RT), Upper = mean(RT) * 1.05, Lower = mean(RT) * 0.95)

## plot the retention time 
p <- ggplot(dat_std, aes(x = Time, y = RT)) +
    facet_wrap(~ Compound, scales = "free") +
    geom_point(aes(color = Type), alpha = 0.9) +
    geom_hline(aes(yintercept = Mean), linetype = "solid",
               color = "gray", linewidth = 0.8) +
    geom_hline(aes(yintercept = Upper), linetype = "dashed",
                color = "gray", linewidth = 0.8) +
    geom_hline(aes(yintercept = Lower), linetype = "dashed",
                color = "gray", linewidth = 0.8) +
    main_theme +
    labs(x = "Measure time", y = "Retention time") +
    theme(legend.position = "top") +
    scale_color_manual(values = c("plum", "#F0CD42", "#56B4E9"))

ggsave("s1_rt.pdf", width = 20,  height = 16)

write.table(dat_std, "s1_data_QC.txt", quote = F, sep = "\t", row.names = F)

