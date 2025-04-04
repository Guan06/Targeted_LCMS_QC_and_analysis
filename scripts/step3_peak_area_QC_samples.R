source("~/Desktop/Cam_2022/20220301_ATP_ADP_data_processing/plot_settings.R")

library(ggplot2)
library(dplyr)

dat <- read.table("./s0_reformat_data.txt", header = T, sep = "\t")
qc <- dat[sapply(dat$Name, function(x){ grepl("QC", x)}), ]
qc <- qc %>% group_by(Compound) %>% 
    mutate(SD = sd(Area), Mean = mean(Area),
            Upper = mean(Area) * 1.25, Lower = mean(Area) * 0.75)

qc$Warning <- ifelse(abs(qc$Area - qc$Mean) > qc$Mean * 0.25, "Yes", "No")
qc$Deviation <- abs(qc$Area - qc$Mean) / qc$Mean

p <- ggplot(qc, aes(x = Time, y = Area)) +
    facet_wrap(~ Compound, scales = "free") +
    geom_point(aes(color = Warning)) + geom_line(color = "gray") + 
    main_theme +
    geom_hline(aes(yintercept = Mean), linetype = "solid",
               color = "gold", linewidth = 0.8) +
    geom_hline(aes(yintercept = Upper), linetype = "dashed",
                color = "gold", linewidth = 0.8) +
    geom_hline(aes(yintercept = Lower), linetype = "dashed",
                color = "gold", linewidth = 0.8) +
    scale_color_manual(values = c("gray", "salmon")) +
    labs(x = "Measure time", y = "Area") +
    theme(legend.position = "top") 

ggsave("s3_area_over_time_facet.pdf", p, height = 10, width = 18)

p2 <- ggplot(qc, aes(x = Time, y = log(Area))) +
  geom_point(aes(color = Compound, shape = Warning), size = 1.8) +
  geom_line(aes(group = Compound, color = Compound)) +
  scale_color_manual(values = getOI(length(unique(dat$Compound)))) +
  scale_shape_manual(values = c(1, 16)) +
  main_theme +
  labs(x = "Measure time", y = "Area (log transformed)") +
  theme(legend.position = "right") 

ggsave("s3_area_over_time_all.pdf", p2, height = 6, width = 10)

write.table(qc, "s3_QC_samples_area_over_time_SD_Mean.txt", quote = F,
            row.names = F, sep = "\t")
