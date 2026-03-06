setwd("/Users/qcluke/Desktop/Chapter 1/Analysis/Passive Sampling")

### mesocosm eDNA conc. ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## data
library(tidyverse)

# total eDNA conc.
mesocosm_qubit = read_csv("Qubit data/Chapter1_Mesocosm_Qubit.csv")
str(mesocosm_qubit)

mesocosm_qubit = as.data.frame(mesocosm_qubit) 
mesocosm_qubit$Time = factor (mesocosm_qubit$Time, levels = c("10m","1h","6h","12h","24h"))
mesocosm_qubit$Extraction = factor (mesocosm_qubit$Extraction, levels = c("PW","PS","BT"))
str(mesocosm_qubit)

mesocosm_qubit_statistics = mesocosm_qubit %>%
  group_by(Time, Extraction) %>%
  summarise(
    Conc_Mean = mean(Qubit_concentration),
    Conc_SD = sd(Qubit_concentration),
    Conc_SEM = sd(Qubit_concentration) / sqrt(6),
    .groups = "drop"
  )

mesocosm_qubit_statistics_time = mesocosm_qubit %>%
  group_by(Time) %>%
  summarise(
    Mean = mean(Qubit_concentration),
    SD = sd(Qubit_concentration),
    SEM = sd(Qubit_concentration) / sqrt(18),
    .groups = "drop"
  )

mesocosm_qubit_statistics_extraction = mesocosm_qubit %>%
  group_by(Extraction) %>%
  summarise(
    Mean = mean(Qubit_concentration),
    SD = sd(Qubit_concentration),
    SEM = sd(Qubit_concentration) / sqrt(30),
    .groups = "drop"
  )

mesocosm_qubit_gam_table = mesocosm_qubit %>%
  mutate(Time = case_match(as.character(Time), "10m" ~ 1/6, "1h" ~ 1, "6h" ~ 6,"12h" ~ 12, "24h" ~ 24))
str(mesocosm_qubit_gam_table)

# prokaryotic eDNA conc.
mesocosm_qpcr_16s = read_csv("qPCR data/Chapter1_Mesocosm_16S.csv")
str(mesocosm_qpcr_16s)

mesocosm_qpcr_16s = as.data.frame(mesocosm_qpcr_16s) 
mesocosm_qpcr_16s$Time = factor (mesocosm_qpcr_16s$Time, levels = c("10m","1h","6h","12h","24h"))
mesocosm_qpcr_16s$Extraction = factor (mesocosm_qpcr_16s$Extraction, levels = c("PW","PS","BT"))
str(mesocosm_qpcr_16s)

mesocosm_qpcr_16s_statistics = mesocosm_qpcr_16s %>% 
  group_by(Time, Extraction) %>%
  summarise(
    Mean = mean(Pro_Concentration),
    SD = sd(Pro_Concentration),
    SEM = sd(Pro_Concentration) / sqrt(6),
    .groups = "drop"
  )

mesocosm_qpcr_16s_statistics_time = mesocosm_qpcr_16s %>%
  group_by(Time) %>%
  summarise(
    Mean = mean(Pro_Concentration),
    SD = sd(Pro_Concentration),
    SEM = sd(Pro_Concentration) / sqrt(18),
    .groups = "drop"
  )

mesocosm_qpcr_16s_statistics_extraction = mesocosm_qpcr_16s %>%
  group_by(Extraction) %>%
  summarise(
    Mean = mean(Pro_Concentration),
    SD = sd(Pro_Concentration),
    SEM = sd(Pro_Concentration) / sqrt(30),
    .groups = "drop"
  )

mesocosm_qpcr_16s_gam_table = mesocosm_qpcr_16s %>%
  mutate(Time = case_match(as.character(Time), "10m" ~ 1/6, "1h" ~ 1, "6h" ~ 6,"12h" ~ 12, "24h" ~ 24))
str(mesocosm_qpcr_16s_gam_table)

# microeukaryotic eDNA conc.
mesocosm_qpcr_18s = read_csv("qPCR data/Chapter1_Mesocosm_18S.csv")

mesocosm_qpcr_18s = as.data.frame(mesocosm_qpcr_18s) 
mesocosm_qpcr_18s$Time = factor (mesocosm_qpcr_18s$Time, levels = c("10m","1h","6h","12h","24h"))
mesocosm_qpcr_18s$Extraction = factor (mesocosm_qpcr_18s$Extraction, levels = c("PW","PS","BT"))
str(mesocosm_qpcr_18s)

mesocosm_qpcr_18s_statistics = mesocosm_qpcr_18s %>% 
  group_by(Time, Extraction) %>%
  summarise(
    Mean = mean(Eu_Concentration),
    SD = sd(Eu_Concentration),
    SEM = sd(Eu_Concentration) / sqrt(6),
    .groups = "drop"
  ) 

mesocosm_qpcr_18s_statistics_time = mesocosm_qpcr_18s %>%
  group_by(Time) %>%
  summarise(
    Mean = mean(Eu_Concentration),
    SD = sd(Eu_Concentration),
    SEM = sd(Eu_Concentration) / sqrt(18),
    .groups = "drop"
  )

mesocosm_qpcr_18s_statistics_extraction = mesocosm_qpcr_18s %>%
  group_by(Extraction) %>%
  summarise(
    Mean = mean(Eu_Concentration),
    SD = sd(Eu_Concentration),
    SEM = sd(Eu_Concentration) / sqrt(30),
    .groups = "drop"
  )

mesocosm_qpcr_18s_gam_table = mesocosm_qpcr_18s %>%
  mutate(Time = case_match(as.character(Time), "10m" ~ 1/6, "1h" ~ 1, "6h" ~ 6,"12h" ~ 12, "24h" ~ 24))
str(mesocosm_qpcr_18s_gam_table)





## Stats
library(car)
library(multcompView)
library(mgcv)

# total eDNA conc.
shapiro.test(mesocosm_qubit$Qubit_concentration)   

leveneTest(mesocosm_qubit$Qubit_concentration ~ mesocosm_qubit$Extraction, data = mesocosm_qubit)
leveneTest(mesocosm_qubit$Qubit_concentration ~ mesocosm_qubit$Time, data = mesocosm_qubit) 

summary(aov(Qubit_concentration~Time * Extraction, data=mesocosm_qubit)) 
TukeyHSD(aov(Qubit_concentration~Time * Extraction, data=mesocosm_qubit))
multcompLetters4(aov(Qubit_concentration~Time * Extraction, data=mesocosm_qubit), TukeyHSD(aov(Qubit_concentration~Time * Extraction, data=mesocosm_qubit)))

summary(gam(Qubit_concentration ~ s(Time,bs="cs", k =3), data=mesocosm_qubit_gam_table)) 

# prokaryotic eDNA conc.
shapiro.test(mesocosm_qpcr_16s$Pro_Concentration) 

leveneTest(mesocosm_qpcr_16s$Pro_Concentration ~ mesocosm_qpcr_16s$Extraction, data = mesocosm_qpcr_16s)
leveneTest(mesocosm_qpcr_16s$Pro_Concentration ~ mesocosm_qpcr_16s$Time, data = mesocosm_qpcr_16s) 

summary(aov(Pro_Concentration~Time * Extraction, data=mesocosm_qpcr_16s)) 
TukeyHSD(aov(Pro_Concentration~Time * Extraction, data=mesocosm_qpcr_16s))
multcompLetters4(aov(Pro_Concentration~Time * Extraction, data=mesocosm_qpcr_16s), TukeyHSD(aov(Pro_Concentration~Time * Extraction, data=mesocosm_qpcr_16s)))

summary(gam(Pro_Concentration ~ s(Time,bs="cs", k =3), data=mesocosm_qpcr_16s_gam_table))

# microeukaryotic eDNA conc.
shapiro.test(mesocosm_qpcr_18s$Eu_Concentration) 

leveneTest(mesocosm_qpcr_18s$Eu_Concentration ~ mesocosm_qpcr_18s$Extraction, data = mesocosm_qpcr_18s)
leveneTest(mesocosm_qpcr_18s$Eu_Concentration ~ mesocosm_qpcr_18s$Time, data = mesocosm_qpcr_18s)

summary(aov(Eu_Concentration~Time * Extraction, data=mesocosm_qpcr_18s))  
TukeyHSD(aov(Eu_Concentration~Time * Extraction, data=mesocosm_qpcr_18s))
multcompLetters4(aov(Eu_Concentration~Time * Extraction, data=mesocosm_qpcr_18s), TukeyHSD(aov(Eu_Concentration~Time * Extraction, data=mesocosm_qpcr_18s)))

summary(gam(Eu_Concentration ~ s(Time,bs="cs", k =3), data=mesocosm_qpcr_18s_gam_table))





## plot
library(tidyverse)

# total eDNA conc.
mesocosm_qubit_figure = ggplot(mesocosm_qubit_statistics, aes(x = Time, y = Conc_Mean, fill = Extraction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_errorbar(aes(ymin = Conc_Mean - Conc_SEM, ymax = Conc_Mean + Conc_SEM, color = factor(Extraction)), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_color_manual(values = c("BT" = "#F27970","PS" = "#BB9727","PW"= "#05B9E2")) +
  scale_fill_manual(values = c("BT" = "#F27970", "PS" = "#BB9727", "PW" = "#05B9E2")) + 
  scale_y_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", " 24 hours")) +
  #ggtitle("(A)") + 
  #xlab("Submersion time") + 
  ylab("Total eDNA conc.\n(ng/μL)") +
  geom_text(x="10m", y=30, label="h h gh", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=30, label="h fgh fgh", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=30, label="fhg efg ef",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=30, label="ef de b",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=30, label="cd bc a",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x=0.5, y=20, label = "(A)", size =3.5, color = "black",family = "Helvetica", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qubit_figure.pdf", plot = mesocosm_qubit_figure, width = 8, height = 5, units = "cm", device = "pdf")

mesocosm_qubit_figure_time = ggplot(mesocosm_qubit_statistics_time, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", " 24 hours")) +
  #ggtitle("(A)") + 
  #xlab("Submersion time") + 
  ylab("Total eDNA conc.\n(ng/μL)") +
  geom_text(x="10m", y=20, label="d", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=20, label="d", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=20, label="c",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=20, label="b",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=20, label="a",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x=0.5, y=20, label = "(A)", size =3.5, color = "black",family = "Helvetica", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qubit_figure_time.png", plot = mesocosm_qubit_figure_time, width = 8, height = 5, units = "cm", device = "pdf")

mesocosm_qubit_figure_extraction = ggplot(mesocosm_qubit_statistics_extraction, aes(x = Extraction, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_discrete(label = c("PowerWater", "PowerSoil", "Blood&Tissue")) +
  #ggtitle("(A)") + 
  #xlab("Submersion time") + 
  ylab("Total eDNA conc.\n(ng/μL)") +
  geom_text(x="PW", y=20, label="c", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="PS", y=20, label="b", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="BT", y=20, label="a",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x=0.5, y=20, label = "(D)", size =3.5, color = "black",family = "Helvetica", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 2.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) 

ggsave("mesocosm_qubit_figure_extraction.png", plot = mesocosm_qubit_figure_extraction, width = 8, height = 5, units = "cm", device = "pdf") 

# prokaryotic eDNA conc.
mesocosm_qpcr_16s_figure = ggplot(mesocosm_qpcr_16s_statistics, aes(x = Time, y = Mean, fill = Extraction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM, color = factor(Extraction)), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_color_manual(values = c("BT" = "#F27970","PS" = "#BB9727","PW"= "#05B9E2")) +
  scale_fill_manual(values = c("BT" = "#F27970", "PS" = "#BB9727", "PW" = "#05B9E2")) + 
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3,4)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", "24 hours")) +
  #ggtitle("(B)") + 
  #xlab("Submersion time") + 
  ylab("Prokaryotic eDNA conc.\n(million copies/μL)") +
  geom_text(x="10m", y=4, label="h h h", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=4, label="gh gh g", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=4, label="fg g ef ",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=4, label="de de b",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=4, label="cd bc a",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x=0.5, y=3, label = "(B)", size =3.5, color = "black",family = "Helvetica", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qpcr_16s_figure.pdf", plot = mesocosm_qpcr_figure, width = 8, height = 5, units = "cm", device = "pdf")

mesocosm_qpcr_16s_figure_time = ggplot(mesocosm_qpcr_16s_statistics_time, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", " 24 hours")) +
  #ggtitle("(A)") + 
  #xlab("Submersion time") + 
  ylab("Prokaryotic eDNA yield\n(million copies/μL)") +
  geom_text(x="10m", y=4, label="e", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=4, label="d", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=4, label="c",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=4, label="b",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=4, label="a",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x=0.5, y=4, label = "(B)", size =3.5, color = "black",family = "Helvetica", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qpcr_16s_figure_time.pdf", plot = mesocosm_qpcr_figure_time, width = 8, height = 5, units = "cm", device = "pdf")

mesocosm_qpcr_16s_figure_extraction = ggplot(mesocosm_qpcr_16s_statistics_extraction, aes(x = Extraction, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) +
  scale_x_discrete(label = c("PowerWater", "PowerSoil", "Blood&Tissue")) +
  #ggtitle("(A)") + 
  #xlab("Submersion time") + 
  ylab("Prokaryotic eDNA yield\n(million copies/μL)") +
  geom_text(x="PW", y=4, label="b", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="PS", y=4, label="b", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="BT", y=4, label="a",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x=0.5, y=4, label = "(E)", size =3.5, color = "black",family = "Helvetica", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 2.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) 

ggsave("mesocosm_qpcr_16s_figure_time_extraction.pdf", plot = mesocosm_qpcr_figure_time_extraction, width = 8, height = 5, units = "cm", device = "pdf") 


# microeukaryotic eDNA conc.
mesocosm_qpcr_18s_figure = ggplot(mesocosm_qpcr_18s_statistics, aes(x = Time, y = Mean, fill = Extraction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM, color = factor(Extraction)), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_color_manual(values = c("BT" = "#F27970","PS" = "#BB9727","PW"= "#05B9E2")) +
  scale_fill_manual(values = c("BT" = "#F27970", "PS" = "#BB9727", "PW" = "#05B9E2")) + 
  scale_y_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", "24 hours")) +
  #ggtitle("(C)") + 
  xlab("Submersion time") + 
  ylab("Microeukaryotic eDNA conc.\n(million copies/μL)") +
  geom_text(x="10m", y=3, label="h h h", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=3, label="g g ef", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=3, label="f ef c",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=3, label="ef de b",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=3, label="cd c a",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x=0.5, y=2, label = "(C)", size =3.5, color = "black",family = "Helvetica", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qpcr_18s_figure.pdf", plot = mesocosm_qpcr_18s_figure, width = 8, height = 5, units = "cm", device = "pdf")

mesocosm_qpcr_18s_figure_time = ggplot(mesocosm_qpcr_18s_statistics_time, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", " 24 hours")) +
  #ggtitle("(A)") + 
  xlab("Submersion time") + 
  ylab("Microeukaryotic eDNA yield\n(million copies/μL)") +
  geom_text(x="10m", y=3, label="e", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=3, label="d", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=3, label="c",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=3, label="b",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=3, label="a",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x=0.5, y=3, label = "(C)", size =3.5, color = "black",family = "Helvetica", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qpcr_18s_figure_time.pdf", plot = mesocosm_qpcr_18s_figure_time, width = 8, height = 5, units = "cm", device = "pdf")

mesocosm_qpcr_18s_figure_extraction =  ggplot(mesocosm_qpcr_18s_statistics_extraction, aes(x = Extraction, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 3), breaks = c(0, 1, 2,3)) +
  scale_x_discrete(label = c("PowerWater", "PowerSoil", "Blood&Tissue")) +
  #ggtitle("(A)") + 
  xlab("Extraction method") + 
  ylab("Microeukaryotic eDNA yield\n(million copies/μL)") +
  geom_text(x="PW", y=3, label="c", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="PS", y=3, label="b", size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="BT", y=3, label="a",size =2.5,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x=0.5, y=3, label = "(F)", size =3.5, color = "black",family = "Helvetica", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 2.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) 

ggsave("mesocosm_qpcr_18s_figure_extraction.pdf", plot = mesocosm_qpcr_18s_figure_extraction, width = 8, height = 5, units = "cm", device = "pdf") 

# three eDNA conc.
figure_mesocosm_conc_gam = ggplot(mesocosm_conc_gam_table, aes(x = Time)) +
  geom_point(aes(y = Qubit_concentration), color = "#999999", size = 2, show.legend = FALSE, alpha = 1) +
  geom_point(aes(y = Pro_Concentration_scaled), color = "#BEB8DC", size = 2, show.legend = FALSE, alpha = 1) +
  geom_point(aes(y = Eu_Concentration_scaled), color = "#8ECFC9", size = 2, show.legend = FALSE, alpha = 1) +
  stat_smooth(aes(y = Qubit_concentration), method = "gam", formula = y ~ s(x, bs = "cs", k = 3), color = "#999999", size = 1.5, se = TRUE, level = 0.95) +
  stat_smooth(aes(y = Pro_Concentration_scaled), method = "gam", formula = y ~ s(x, bs = "cs", k = 3), color = "#BEB8DC", size = 1.5, se = TRUE, level = 0.95) +
  stat_smooth(aes(y = Eu_Concentration_scaled), method = "gam", formula = y ~ s(x, bs = "cs", k = 5), color = "#8ECFC9", size = 1.5, se = TRUE, level = 0.95) + 
  scale_y_continuous(name = "Total eDNA conc.\n(ng/µL)", limits = c(0, 40), breaks = c(0, 10, 20, 30, 40), sec.axis = sec_axis(~ . /10,  breaks = c(0, 1, 2,3,4), name = "Prokaryotic/Micoeukaryotic\neDNA conc. (milliion copies/µL)")) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 6, 12, 18, 24)) +
  ggtitle("") + xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=10,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("figure_mesocosm_conc_gam.png", plot = figure_mesocosm_conc_gam, width = 8, height = 8, units = "cm", dpi = 600)



## figure-mesocosm conc.------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)

figure_mesocosm_conc = mesocosm_qubit_figure + mesocosm_qpcr_16s_figure + mesocosm_qpcr_18s_figure + plot_layout(ncol = 1, nrow = 3, byrow = T)
ggsave("figure_mesocosm_conc.pdf", plot = figure_mesocosm_conc, width = 8, height = 15, units = "cm", device = "pdf")

figure_mesocosm_conc_time_extraction = mesocosm_qubit_figure_time + mesocosm_qubit_figure_extraction + mesocosm_qpcr_16s_figure_time + mesocosm_qpcr_16s_figure_extraction + mesocosm_qpcr_18s_figure_time + mesocosm_qpcr_18s_figure_extraction + plot_layout(ncol = 2, nrow = 3, byrow = T)
ggsave("figure_mesocosm_conc_time_extraction.pdf", plot = figure_mesocosm_conc_time_extraction, width = 16, height = 15, units = "cm", device = "pdf")




















### asv and tax tables tidy up 16s ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(vegan)

## import tables
asv_16s = read_tsv("Sequence data/asv_16s.tsv")
dim(asv_16s)

tax_16s = read_tsv("Sequence data/tax_16s.tsv")
dim(tax_16s)

## merge tables
asv_tax_16s = inner_join(tax_16s,asv_16s, by = "Feature ID") 
dim(asv_tax_16s) 

## remove Chloroplast, Mitochondria and unassigned ASVs
asv_tax_16s_filtered = asv_tax_16s_filtered[!str_detect(asv_tax_16s_filtered$Taxon, "Chloroplast"),] 
dim(asv_tax_16s_filtered) 
asv_tax_16s_filtered = asv_tax_16s_filtered[!str_detect(asv_tax_16s_filtered$Taxon, "Mitochondria"),] 
dim(asv_tax_16s_filtered) 
asv_tax_16s_filtered = asv_tax_16s_filtered[!str_detect(asv_tax_16s_filtered$Taxon, "Unassigned"),] 
dim(asv_tax_16s_filtered) 
dim(asv_tax_16s_filtered) 

## decontamination 
asv_tax_16s_PCR_NC = asv_tax_16s_filtered %>% 
  select (`PCR-NEG-1-1`,`PCR-NEG-1-2`,`PCR-NEG-1-3`) %>%
  rowwise() %>%
  mutate(PCR_NC_Max = max(`PCR-NEG-1-1`, `PCR-NEG-1-2`, `PCR-NEG-1-3`)) 

asv_tax_16s_decontamed = asv_tax_16s_filtered %>% 
  select(-"PCR-NEG-1-1", -"PCR-NEG-1-2", -"PCR-NEG-1-3") %>%
  mutate(across(matches("^F-|^S-|^C1|^C2|^C3|^T1|^T2|^T3"), ~ pmax(.-asv_tax_16s_PCR_NC$PCR_NC_Max, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("F-C1"), ~ pmax(.-asv_tax_16s_decontamed$`F-C1-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("F-C2"), ~ pmax(.-asv_tax_16s_decontamed$`F-C2-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("F-C3"), ~ pmax(.-asv_tax_16s_decontamed$`F-C3-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("F-T1"), ~ pmax(.-asv_tax_16s_decontamed$`F-T1-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("F-T2"), ~ pmax(.-asv_tax_16s_decontamed$`F-T2-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("F-T3"), ~ pmax(.-asv_tax_16s_decontamed$`F-T3-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("S-C1"), ~ pmax(.-asv_tax_16s_decontamed$`S-C1-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("S-C2"), ~ pmax(.-asv_tax_16s_decontamed$`S-C2-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("S-C3"), ~ pmax(.-asv_tax_16s_decontamed$`S-C3-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("S-T1"), ~ pmax(.-asv_tax_16s_decontamed$`S-T1-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("S-T2"), ~ pmax(.-asv_tax_16s_decontamed$`S-T2-NC`, 0)))

asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% 
  mutate(across(starts_with("S-T3"), ~ pmax(.-asv_tax_16s_decontamed$`S-T3-NC`, 0)))

## remove NC, ID, confidence columns
asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% select(-contains("-NC"))
dim(asv_tax_16s_decontamed)

asv_tax_16s_decontamed = as.data.frame(asv_tax_16s_decontamed) 
rownames(asv_tax_16s_decontamed) = asv_tax_16s_decontamed$`Feature ID`
asv_tax_16s_decontamed = asv_tax_16s_decontamed[,-1]
asv_tax_16s_decontamed = asv_tax_16s_decontamed[,-2]
dim(asv_tax_16s_decontamed)

## produce rarefaction figure
library(vegan)
rarecurve(t(asv_tax_16s_decontamed[,-1]), step=100, lwd=2, ylab="No. of prokaryptic ASVs", xlab="No. of sequences",label =F)

## Separate to mesocosm and field table, remove no longer existed ASVs
asv_16s_clean_mesocosm = asv_tax_16s_decontamed %>%
  select(starts_with(c("C","T")),-Taxon)
dim(asv_16s_clean_mesocosm) 
asv_16s_clean_mesocosm = asv_16s_clean_mesocosm[which(rowSums(asv_16s_clean_mesocosm)>0),] 
dim(asv_16s_clean_mesocosm)
asv_16s_clean_mesocosm = t(asv_16s_clean_mesocosm)
dim(asv_16s_clean_mesocosm)

asv_16s_clean_field = asv_tax_16s_decontamed %>%
  select(starts_with(c("F","S")))
dim(asv_16s_clean_field)
asv_16s_clean_field = asv_16s_clean_field[which(rowSums(asv_16s_clean_field)>0),] 
dim(asv_16s_clean_field)
asv_16s_clean_field = t(asv_16s_clean_field)
dim(asv_16s_clean_field)

## rarefy ASV tables, remove no longer existed ASVs
asv_16s_clean_mesocosm_rarefy = rrarefy(asv_16s_clean_mesocosm, min(rowSums(asv_16s_clean_mesocosm)))
dim(asv_16s_clean_mesocosm_rarefy) 
asv_16s_clean_mesocosm_rarefy = asv_16s_clean_mesocosm_rarefy[,which(colSums(asv_16s_clean_mesocosm_rarefy)>0)] 
dim(asv_16s_clean_mesocosm_rarefy) 

asv_16s_clean_field_rarefy = rrarefy(asv_16s_clean_field, min(rowSums(asv_16s_clean_field)))
dim(asv_16s_clean_field_rarefy) 
asv_16s_clean_field_rarefy = asv_16s_clean_field_rarefy[,which(colSums(asv_16s_clean_field_rarefy)>0)] 
dim(asv_16s_clean_field_rarefy) 






























### asv and tax tables tidy up 18s ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(vegan)

## import tables
asv_18s = read_tsv("Sequence data/asv_18s.tsv")
dim(asv_18s)

tax_18s = read_tsv("Sequence data/tax_18s.tsv")
dim(tax_18s)

## merge tables
asv_tax_18s = inner_join(tax_18s,asv_18s, by = "Feature ID") 
dim(asv_tax_18s) 

# remove metazoans,land plants,and unassigned ASVs
asv_tax_18s_filtered = asv_tax_18s[!str_detect(asv_tax_18s$Taxon, "Metazoa"),] 
dim(asv_tax_18s_filtered) 
asv_tax_18s_filtered = asv_tax_18s_filtered[!str_detect(asv_tax_18s_filtered$Taxon, "Streptophyta"),] 
dim(asv_tax_18s_filtered) 
asv_tax_18s_filtered = asv_tax_18s_filtered[!str_detect(asv_tax_18s_filtered$Taxon, "Unassigned"),] 
dim(asv_tax_18s_filtered) 
dim(asv_tax_18s_filtered)

## decontamination 
# no ASV existed in 3 PCR negative controls
asv_tax_18s_decontamed = asv_tax_18s_filtered %>% 
  mutate(across(starts_with("F-C1"), ~ pmax(.-asv_tax_18s_filtered$`F-C1-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("F-C2"), ~ pmax(.-asv_tax_18s_decontamed$`F-C2-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("F-C3"), ~ pmax(.-asv_tax_18s_decontamed$`F-C3-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("F-T1"), ~ pmax(.-asv_tax_18s_decontamed$`F-T1-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("F-T2"), ~ pmax(.-asv_tax_18s_decontamed$`F-T2-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("F-T3"), ~ pmax(.-asv_tax_18s_decontamed$`F-T3-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("S-C1"), ~ pmax(.-asv_tax_18s_decontamed$`S-C1-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("S-C2"), ~ pmax(.-asv_tax_18s_decontamed$`S-C2-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("S-C3"), ~ pmax(.-asv_tax_18s_decontamed$`S-C3-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("S-T1"), ~ pmax(.-asv_tax_18s_decontamed$`S-T1-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("S-T2"), ~ pmax(.-asv_tax_18s_decontamed$`S-T2-NC`, 0)))

asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% 
  mutate(across(starts_with("S-T3"), ~ pmax(.-asv_tax_18s_decontamed$`S-T3-NC`, 0)))

## remove NC, ID, confidence columns
asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% select(-contains("-NC"))
dim(asv_tax_18s_decontamed)

asv_tax_18s_decontamed = as.data.frame(asv_tax_18s_decontamed) 
rownames(asv_tax_18s_decontamed) = asv_tax_18s_decontamed$`Feature ID`
asv_tax_18s_decontamed = asv_tax_18s_decontamed[,-1]
asv_tax_18s_decontamed = asv_tax_18s_decontamed[,-2]
dim(asv_tax_18s_decontamed)

## produce rarefaction figure
library(vegan)
rarecurve(t(asv_tax_18s_decontamed[,-1]), step=100, lwd=2, ylab="No. of microeukaryptic ASVs", xlab="No. of sequences",label =F)

## Separate to mesocosm and field table, remove no longer existed ASVs
asv_18s_clean_mesocosm = asv_tax_18s_decontamed %>%
  select(starts_with(c("C","T")),-Taxon)
dim(asv_18s_clean_mesocosm) 
asv_18s_clean_mesocosm = asv_18s_clean_mesocosm[which(rowSums(asv_18s_clean_mesocosm)>0),] 
dim(asv_18s_clean_mesocosm)
asv_18s_clean_mesocosm = t(asv_18s_clean_mesocosm)
dim(asv_18s_clean_mesocosm)

asv_18s_clean_field = asv_tax_18s_decontamed %>%
  select(starts_with(c("F","S")))
dim(asv_18s_clean_field)
asv_18s_clean_field = asv_18s_clean_field[which(rowSums(asv_18s_clean_field)>0),] 
dim(asv_18s_clean_field)
asv_18s_clean_field = t(asv_18s_clean_field)
dim(asv_18s_clean_field)

## rarefy ASV tables, remove no longer existed ASVs
asv_18s_clean_mesocosm_rarefy = rrarefy(asv_18s_clean_mesocosm, min(rowSums(asv_18s_clean_mesocosm)))
dim(asv_18s_clean_mesocosm_rarefy) 
asv_18s_clean_mesocosm_rarefy = asv_18s_clean_mesocosm_rarefy[,which(colSums(asv_18s_clean_mesocosm_rarefy)>0)] 
dim(asv_18s_clean_mesocosm_rarefy) 

asv_18s_clean_field_rarefy = rrarefy(asv_18s_clean_field, min(rowSums(asv_18s_clean_field)))
dim(asv_18s_clean_field_rarefy) 
asv_18s_clean_field_rarefy = asv_18s_clean_field_rarefy[,which(colSums(asv_18s_clean_field_rarefy)>0)] 
dim(asv_18s_clean_field_rarefy) 
























### metadata------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mesocosm_metadata = asv_16s_clean_mesocosm %>% rownames() %>% as_tibble() %>% 
  mutate (Sample = value) %>%
  separate(col = value, into = c("Water", "Time", "ExtractionMethod"), sep = "-", extra = "merge") %>%
  separate(col = Water, into = c("Water", "Tank"), sep = "(?<=\\D)(?=\\d)", extra = "merge") %>%
  select(Sample, Water, Tank, Time, ExtractionMethod)
mesocosm_metadata$Water = factor(mesocosm_metadata$Water)
mesocosm_metadata$Time = factor(mesocosm_metadata$Time, levels = c("10m","1h","6h","12h","24h"))
mesocosm_metadata$ExtractionMethod = factor(mesocosm_metadata$ExtractionMethod)
str(mesocosm_metadata)

field_metadata = asv_16s_clean_field %>% rownames() %>% as_tibble() %>% 
  mutate (Sample = value) %>%
  separate(col = value, into = c("Type", "Water", "Repetition"), sep = "-", extra = "merge") %>%
  separate(col = Water, into = c("Water", "Site"), sep = "(?<=\\D)(?=\\d)", extra = "merge") %>% 
  mutate(Sampling_Water = paste(Type, Water, sep = "_")) %>% 
  select(Sample, Type, Water, Site, Repetition, Sampling_Water) 
field_metadata$Type = factor(field_metadata$Type)
field_metadata$Water = factor(field_metadata$Water)
field_metadata$Site = factor(field_metadata$Site)
field_metadata$Sampling_Water = factor(field_metadata$Sampling_Water)
str(field_metadata)




















### mesocosm diversity ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## data 
library(tidyverse)
library(vegan)
library(picante)

# 16s
mesocosm_16s_richness = estimateR(asv_16s_clean_mesocosm_rarefy)[1,]
mesocosm_16s_richness = data.frame(Sample = names(mesocosm_16s_richness), Richness = mesocosm_16s_richness)
mesocosm_16s_richness = left_join(mesocosm_metadata, mesocosm_16s_richness, by = "Sample")
str(mesocosm_16s_richness)
mesocosm_16s_richness = as.data.frame(mesocosm_16s_richness) 
str(mesocosm_16s_richness)

tree_rooted_16sv4 = read.tree(file = "Sequence data/16sv4-rooted-tree.nwk")
mesocosm_16s_pd = pd(asv_16s_clean_mesocosm_rarefy, tree_rooted_16sv4, include.root = TRUE)
mesocosm_16s_pd = data.frame(Sample = rownames(mesocosm_16s_pd), PD = mesocosm_16s_pd$PD)
mesocosm_16s_pd = left_join(mesocosm_metadata, mesocosm_16s_pd, by = "Sample")
str(mesocosm_16s_pd)
mesocosm_16s_pd = as.data.frame(mesocosm_16s_pd) 
str(mesocosm_16s_pd) 

mesocosm_16s_richness_statistics = mesocosm_16s_richness %>% group_by(Time) %>%
  summarise(
    Mean = mean(Richness),
    SD = sd(Richness),
    SEM = sd(Richness) / sqrt(6),
    .groups = "drop"
  ) 

library(tidyverse)
mesocosm_16s_pd_statistics = mesocosm_16s_pd %>% group_by(Time) %>%
  summarise(
    Mean = mean(PD),
    SD = sd(PD),
    SEM = sd(PD) / sqrt(6),
    .groups = "drop"
  ) 

mesocosm_16s_richness_gam_table = mesocosm_16s_richness %>%
  mutate(Time = case_match(as.character(Time), "10m" ~ 1/6, "1h" ~ 1, "6h" ~ 6,"12h" ~ 12, "24h" ~ 24))
str(mesocosm_16s_richness_gam_table)

mesocosm_16s_pd_gam_table = mesocosm_16s_pd %>%
  mutate(Time = case_match(as.character(Time), "10m" ~ 1/6, "1h" ~ 1, "6h" ~ 6,"12h" ~ 12, "24h" ~ 24))
str(mesocosm_16s_pd_gam_table)


# 18s
mesocosm_18s_richness = estimateR(asv_18s_clean_mesocosm_rarefy)[1,]
mesocosm_18s_richness = data.frame(Sample = names(mesocosm_18s_richness), Richness = mesocosm_18s_richness)
mesocosm_18s_richness = left_join(mesocosm_metadata, mesocosm_18s_richness, by = "Sample")
str(mesocosm_18s_richness)
mesocosm_18s_richness = as.data.frame(mesocosm_18s_richness) 
str(mesocosm_18s_richness)

tree_rooted_18sv9 = read.tree(file = "Sequence data/18sv9-rooted-tree.nwk")
mesocosm_18s_pd = pd(asv_18s_clean_mesocosm_rarefy, tree_rooted_18sv9, include.root = TRUE)
mesocosm_18s_pd = data.frame(Sample = rownames(mesocosm_18s_pd), PD = mesocosm_18s_pd$PD)
mesocosm_18s_pd = left_join(mesocosm_metadata, mesocosm_18s_pd, by = "Sample")
str(mesocosm_18s_pd)
mesocosm_18s_pd = as.data.frame(mesocosm_18s_pd) 
str(mesocosm_18s_pd)

mesocosm_18s_richness_gam_table = mesocosm_18s_richness %>%
  mutate(Time = case_match(as.character(Time), "10m" ~ 1/6, "1h" ~ 1, "6h" ~ 6,"12h" ~ 12, "24h" ~ 24))
str(mesocosm_18s_richness_gam_table)

mesocosm_18s_pd_gam_table = mesocosm_18s_pd %>%
  mutate(Time = case_match(as.character(Time), "10m" ~ 1/6, "1h" ~ 1, "6h" ~ 6,"12h" ~ 12, "24h" ~ 24))
str(mesocosm_18s_pd_gam_table)

mesocosm_18s_richness_statistics = mesocosm_18s_richness %>% group_by(Time) %>%
  summarise(
    Mean = mean(Richness),
    SD = sd(Richness),
    SEM = sd(Richness) / sqrt(6),
    .groups = "drop"
  ) 

library(tidyverse)
mesocosm_18s_pd_statistics = mesocosm_18s_pd %>% group_by(Time) %>%
  summarise(
    Mean = mean(PD),
    SD = sd(PD),
    SEM = sd(PD) / sqrt(6),
    .groups = "drop"
  ) 



## stats
library(car)
library(multcompView)
library(mgcv)

# 16s
shapiro.test(mesocosm_16s_richness$Richness) 
shapiro.test(mesocosm_16s_pd$PD) 

leveneTest(mesocosm_16s_richness$Richness ~ mesocosm_16s_richness$Time, data = mesocosm_16s_richness) 
leveneTest(mesocosm_16s_pd$PD ~ mesocosm_16s_pd$Time, data = mesocosm_16s_pd) 

summary(aov(Richness ~ Time, data = mesocosm_16s_richness)) 
TukeyHSD(aov(Richness ~ Time, data = mesocosm_16s_richness))
summary(aov(PD ~ Time, data = mesocosm_16s_pd)) 
TukeyHSD(aov(PD ~ Time, data = mesocosm_16s_pd))

multcompLetters4(aov(Richness~Time, data=mesocosm_16s_richness), TukeyHSD(aov(Richness~Time, data=mesocosm_16s_richness)))
multcompLetters4(aov(PD~Time, data=mesocosm_16s_pd), TukeyHSD(aov(PD~Time, data=mesocosm_16s_pd)))

summary(gam(Richness ~ s(Time,bs="cs", k =3), data=mesocosm_16s_richness_gam_table)) 
summary(gam(PD ~ s(Time,bs="cs", k =3), data=mesocosm_16s_pd_gam_table)) 

# 18s
shapiro.test(mesocosm_18s_richness$Richness) 
shapiro.test(mesocosm_18s_pd$PD) 

leveneTest(mesocosm_18s_richness$Richness ~ mesocosm_18s_richness$Time, data = mesocosm_18s_richness) 
leveneTest(mesocosm_18s_pd$PD ~ mesocosm_18s_pd$Time, data = mesocosm_18s_pd) 

summary(aov(Richness~Time, data=mesocosm_18s_richness))
TukeyHSD(aov(Richness~Time, data=mesocosm_18s_richness))
summary(aov(PD~Time, data=mesocosm_18s_pd))
TukeyHSD(aov(PD~Time, data=mesocosm_18s_pd))

multcompLetters4(aov(Richness~Time, data=mesocosm_18s_richness), TukeyHSD(aov(Richness~Time, data=mesocosm_18s_richness)))
multcompLetters4(aov(PD~Time, data=mesocosm_18s_pd), TukeyHSD(aov(PD~Time, data=mesocosm_18s_pd)))

summary(gam(Richness ~ s(Time,bs="cs", k =3), data=mesocosm_18s_richness_gam_table)) 
summary(gam(PD ~ s(Time,bs="cs", k =3), data=mesocosm_18s_pd_gam_table))



## plot
library(tidyverse)
# 16s
mesocosm_16s_richness_figure_histogram = ggplot(mesocosm_16s_richness_statistics, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6,fill = "black",color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),color = "black", position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 2000), breaks = c(0, 500, 1000, 1500, 2000)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", "24 hours")) +
  ggtitle("") + 
  xlab("Submersion time") + 
  ylab("Number of prokaryotic ASVs") +
  geom_text(x="10m", y=2000, label="c",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=2000, label="b",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=2000, label="ab",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=2000, label="a",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=2000, label="a",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_16s_richness_figure_histogram.png", plot = mesocosm_16s_richness_figure_histogram, width = 8, height = 5, units = "cm", device = "pdf")

mesocosm_16s_pd_figure_histogram = ggplot(mesocosm_16s_pd_statistics, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6,fill = "black",color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),color = "black", position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 140), breaks = c(0, 20, 40, 60, 80, 100, 100, 120, 140)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", "24 hours")) +
  ggtitle("") + 
  xlab("Submersion time") + 
  ylab("Prokaryotic Faith's PD") +
  geom_text(x="10m", y=130, label="c",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=130, label="b",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=130, label="b",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=130, label="ab",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=130, label="a",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_16s_pd_figure_histogram.pdf", plot = mesocosm_16s_pd_figure_histogram, width = 8, height = 5, units = "cm", device = "pdf")

# 18s
mesocosm_18s_richness_figure_histogram = ggplot(mesocosm_18s_richness_statistics, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6,fill = "black",color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),color = "black", position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 400), breaks = c(0, 100, 200, 300, 400)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", "24 hours")) +
  ggtitle("") + 
  xlab("Submersion time") + 
  ylab("Number of microeukaryotic ASVs") +
  geom_text(x="10m", y=400, label="c",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=400, label="b",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=400, label="ab",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=400, label="ab",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=400, label="a",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_18s_richness_figure_histogram.pdf", plot = mesocosm_18s_richness_figure_histogram, width = 8, height = 5, units = "cm", device = "pdf")

mesocosm_18s_pd_figure_histogram = ggplot(mesocosm_18s_pd_statistics, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6,fill = "black",color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),color = "black", position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", "24 hours")) +
  ggtitle("") + 
  xlab("Submersion time") + 
  ylab("Microeukaryotic Faith's PD") +
  geom_text(x="10m", y=90, label="d",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="1h", y=90, label="c",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="6h", y=90, label="bc",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="12h", y=90, label="ab",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  geom_text(x="24h", y=90, label="a",size =3,color = "black",family = "Helvetica", fontface = "plain") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_18s_pd_figure_histogram.pdf", plot = mesocosm_18s_pd_figure_histogram, width = 8, height = 5, units = "cm", device = "pdf")

# two diveristy
figure_mesocosm_16s_18s_richness_gam = ggplot(mesocosm_16s_18s_richness_gam_table, aes(x = Time)) +
  geom_point(aes(y = Pro), color = "#BEB8DC", size = 2, show.legend = FALSE, alpha = 1) +
  geom_point(aes(y = Eu_scaled), color = "#8ECFC9", size = 2, show.legend = FALSE, alpha = 1) +
  stat_smooth(aes(y = Pro), method = "gam", formula = y ~ s(x, bs = "cs", k = 3), color = "#BEB8DC", size = 1.5, se = TRUE, level = 0.95) +
  stat_smooth(aes(y = Eu_scaled), method = "gam", formula = y ~ s(x, bs = "cs", k = 3), color = "#8ECFC9", size = 1.5, se = TRUE, level = 0.95) +
  #scale_y_continuous(name = "Observed ", limits = c(0, 2400), breaks = c(0, 600, 1200, 1800, 2400)) +
  scale_y_continuous(name = "Prokaryotic", limits = c(0, 2400), breaks = c(0, 600, 1200, 1800, 2400), sec.axis = sec_axis(~ . /3, breaks = c(0, 200, 400,600, 800), name = "Microeukaryotic")) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 6, 12, 18, 24)) +
  ggtitle("") + xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y =element_text(color="black",size=10,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("figure_mesocosm_16s_18s_richness_gam.pdf", plot = figure_mesocosm_16s_18s_richness_gam, width = 8, height = 8, units = "cm", device = "pdf")

figure_mesocosm_16s_18s_pd_gam = ggplot(mesocosm_16s_18s_pd_gam_table, aes(x = Time)) +
  geom_point(aes(y = Pro), color = "#BEB8DC", size = 2, show.legend = FALSE, alpha = 1) +
  geom_point(aes(y = Eu), color = "#8ECFC9", size = 2, show.legend = FALSE, alpha = 1) +
  stat_smooth(aes(y = Pro), method = "gam", formula = y ~ s(x, bs = "cs", k = 3), color = "#BEB8DC", size = 1.5, se = TRUE, level = 0.95) +
  stat_smooth(aes(y = Eu), method = "gam", formula = y ~ s(x, bs = "cs", k = 3), color = "#8ECFC9", size = 1.5, se = TRUE, level = 0.95) +
  #scale_y_continuous(name = "PD ", limits = c(0, 160), breaks = c(0, 40, 80, 120, 160)) +
  scale_y_continuous(name = "Prokaryotic", limits = c(0, 160), breaks = c(0, 40, 80, 120, 160), sec.axis = sec_axis(~ . /1, breaks = c(0, 40, 80, 120, 160), name = "Microeukaryotic")) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 6, 12, 18, 24)) +
  ggtitle("") + xlab("Submersion time (hour)") + ylab("") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=10,family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=10,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5))  +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("figure_mesocosm_16s_18s_pd_gam.pdf", plot = figure_mesocosm_16s_18s_pd_gam, width = 8, height = 8, units = "cm", device = "pdf")



## figure-mesocosm gam------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## data
library(patchwork)

figure_mesocosm_richness_histogram = mesocosm_16s_richness_figure_histogram + mesocosm_18s_richness_figure_histogram + 
                                     mesocosm_16s_pd_figure_histogram + mesocosm_18s_pd_figure_histogram + plot_layout(ncol = 2, nrow = 2, byrow = T)
ggsave("figure_mesocosm_richness_histogram.pdf", plot = figure_mesocosm_richness_histogram, width = 16, height = 10, units = "cm", device = "pdf")


figure_mesocosm_conc_richness_gam= figure_mesocosm_conc_gam + figure_mesocosm_16s_18s_richness_gam + figure_mesocosm_16s_18s_pd_gam + plot_layout(ncol = 1, nrow = 3, byrow = T)
ggsave("figure_mesocosm_conc_richness_gam.pdf", plot = figure_mesocosm_conc_richness_gam, width = 8, height = 19, units = "cm", device = "pdf")





























### field total eDNA conc.  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## data
library(tidyverse)

field_qubit = read_csv("Qubit data/Chapter1_Field_Qubit.csv")

field_qubit = as.data.frame(field_qubit) 
field_qubit$Water = factor (field_qubit$Water, levels = c("Clear","Turbid"))
field_qubit$Site = factor (field_qubit$Site, levels = c("1","2","3"))
field_qubit$Sampling = factor (field_qubit$Sampling, levels = c("Active","Passive"))
str(field_qubit)

field_qubit_statistics = field_qubit %>% 
  group_by(Water, Site, Sampling) %>%
  summarise(
    Mean = mean(Qubit_concentration),
    SD = sd(Qubit_concentration),
    SEM = sd(Qubit_concentration) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 

field_qubit_statistics$Water_Site = factor(field_qubit_statistics$Water_Site)
str(field_qubit_statistics)

## stats
library(car)

shapiro.test(field_qubit_statistics$Mean) 

leveneTest(field_qubit_statistics$Mean ~ field_qubit_statistics$Sampling, data = field_qubit_statistics) 

t.test(subset(field_qubit_statistics, Sampling == "Active")$Mean, subset(field_qubit_statistics, Sampling == "Passive")$Mean, paired = TRUE) 

wilcox.test(field_qubit_statistics$Mean[field_qubit_statistics$Sampling == "Active"], 
            field_qubit_statistics$Mean[field_qubit_statistics$Sampling == "Passive"], paired = TRUE) 

## plot
library(tidyverse)
library(ggsignif)

field_qubit_figure =ggplot(field_qubit_statistics, aes(x=Sampling, y=Mean)) +
  geom_point(size = 3, aes(color = factor(Sampling), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color =factor(Sampling)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("Active" = "#1F77B4", "Passive" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("Active", "Passive")), annotations="*", y_position = 45,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 50), breaks = c(0, 10, 20, 30, 40, 50)) +
  ggtitle("") + xlab("") + ylab("Total eDNA conc. (ng/μL)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "Helvetica", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_qubit_figure.pdf", plot = field_qubit_figure, width = 4, height = 5, units = "cm", device = "pdf") 




















### field prokaryotic eDNA conc. ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## data
library(tidyverse)

field_qpcr_16s = read_csv("qPCR data/Chapter1_Field_16S.csv")

field_qpcr_16s = as.data.frame(field_qpcr_16s) 
field_qpcr_16s$Water = factor (field_qpcr_16s$Water, levels = c("Clear","Turbid"))
field_qpcr_16s$Site = factor (field_qpcr_16s$Site, levels = c("1","2","3"))
field_qpcr_16s$Sampling = factor (field_qpcr_16s$Sampling, levels = c("Active","Passive"))
str(field_qpcr_16s)

field_qpcr_16s_statistics = field_qpcr_16s %>% 
  group_by(Water, Site, Sampling) %>%
  summarise(
    Mean = mean(Pro_Concentration),
    SD = sd(Pro_Concentration),
    SEM = sd(Pro_Concentration) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 

field_qpcr_16s_statistics$Water_Site = factor(field_qpcr_16s_statistics$Water_Site)
str(field_qpcr_16s_statistics)


## stats
library(car)
shapiro.test(field_qpcr_16s_statistics$Mean) 

leveneTest(field_qpcr_16s_statistics$Mean ~ field_qpcr_16s_statistics$Sampling, data = field_qpcr_16s_statistics) 

t.test(subset(field_qpcr_16s_statistics, Sampling == "Active")$Mean, subset(field_qpcr_16s_statistics, Sampling == "Passive")$Mean, paired = TRUE) 

wilcox.test(field_qpcr_16s_statistics$Mean[field_qpcr_16s_statistics$Sampling == "Active"], 
            field_qpcr_16s_statistics$Mean[field_qpcr_16s_statistics$Sampling == "Passive"], paired = TRUE) 

# plot
library(tidyverse)
library(ggsignif)

field_qpcr_16s_figure = ggplot(field_qpcr_16s_statistics, aes(x=Sampling, y=Mean)) +
  geom_point(size = 3, aes(color = factor(Sampling), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color= factor(Sampling)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("Active" = "#1F77B4", "Passive" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("Active", "Passive")), annotations="*", y_position = 83,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 90), breaks = c(0, 30, 60, 90)) +
  ggtitle("") + xlab("") + ylab("Prokaryotic eDNA conc.\n(million copies/μL)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "Helvetica", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_qpcr_16s_figure.pdf", plot = field_qubit_figure, width = 4, height = 5, units = "cm", device = "pdf") 
















### field microeukaryotic eDNA conc. ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## data
library(tidyverse)

field_qpcr_18s = read_csv("qPCR data/Chapter1_Field_18S.csv")

field_qpcr_18s = as.data.frame(field_qpcr_18s) 
field_qpcr_18s$Water = factor (field_qpcr_18s$Water, levels = c("Clear","Turbid"))
field_qpcr_18s$Site = factor (field_qpcr_18s$Site, levels = c("1","2","3"))
field_qpcr_18s$Sampling = factor (field_qpcr_18s$Sampling, levels = c("Active","Passive"))
str(field_qpcr_18s)

field_qpcr_18s_statistics = field_qpcr_18s %>% 
  group_by(Water, Site, Sampling) %>%
  summarise(
    Mean = mean(Eu_Concentration),
    SD = sd(Eu_Concentration),
    SEM = sd(Eu_Concentration) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 

field_qpcr_18s_statistics$Water_Site = factor(field_qpcr_18s_statistics$Water_Site)
str(field_qpcr_18s_statistics)

## stats
library(car)

shapiro.test(field_qpcr_18s_statistics$Mean) 

t.test(subset(field_qpcr_18s_statistics, Sampling == "Active")$Mean, subset(field_qpcr_18s_statistics, Sampling == "Passive")$Mean, paired = TRUE) 

leveneTest(field_qpcr_18s_statistics$Mean ~ field_qpcr_18s_statistics$Sampling, data = field_qpcr_18s_statistics) 

wilcox.test(field_qpcr_18s_statistics$Mean[field_qpcr_18s_statistics$Sampling == "Active"], 
            field_qpcr_18s_statistics$Mean[field_qpcr_18s_statistics$Sampling == "Passive"], paired = TRUE) 

# plot
library(tidyverse)
library(ggsignif)

field_qpcr_18s_figure = ggplot(field_qpcr_18s_statistics, aes(x=Sampling, y=Mean)) +
  geom_point(size = 3, aes(color = factor(Sampling), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color= factor(Sampling)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("Active" = "#1F77B4", "Passive" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("Active", "Passive")), annotations="*", y_position = 5.5,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6)) +
  ggtitle("") + xlab("") + ylab("Microeukaryotic eDNA conc.\n(million copies/μL)") +
  #geom_text(x=0.5, y=5, label = "(C)", size =3, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "Helvetica", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_qpcr_18s_figure.pdf", plot = field_qubit_figure, width = 4, height = 5, units = "cm", device = "pdf") 

















### field prokaryotic diveristy ------------------------------------------------------------------------------------------------------------------------------------------------------------------
## data
library(vegan)
library(picante)
field_16s_richness = estimateR(asv_16s_clean_field_rarefy)[1,]
field_16s_richness = data.frame(Sample = names(field_16s_richness), Richness = field_16s_richness)
field_16s_richness = left_join(field_metadata, field_16s_richness, by = "Sample")
str(field_16s_richness)
field_16s_richness = as.data.frame(field_16s_richness) 
str(field_16s_richness)


tree_rooted_16sv4 = read.tree(file = "Sequence data/16sv4-rooted-tree.nwk")
field_16s_pd = pd(asv_16s_clean_field_rarefy, tree_rooted_16sv4, include.root = TRUE)
field_16s_pd = data.frame(Sample = rownames(field_16s_pd), PD = field_16s_pd$PD)
field_16s_pd = left_join(field_metadata, field_16s_pd, by = "Sample")
str(field_16s_pd)
field_16s_pd = as.data.frame(field_16s_pd) 
str(field_16s_pd)

library(tidyverse)
field_16s_richness_statistics = field_16s_richness %>% group_by(Water, Site, Type) %>%
  summarise(
    Mean = mean(Richness),
    SD = sd(Richness),
    SEM = sd(Richness) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 
str(field_16s_richness_statistics)

field_16s_richness_statistics$Water = factor(field_16s_richness_statistics$Water)
field_16s_richness_statistics$Site = factor(field_16s_richness_statistics$Site)
field_16s_richness_statistics$Type = factor(field_16s_richness_statistics$Type)
field_16s_richness_statistics$Water_Site = factor(field_16s_richness_statistics$Water_Site)
str(field_16s_richness_statistics)

field_16s_pd_statistics = field_16s_pd %>% group_by(Water, Site, Type) %>%
  summarise(
    Mean = mean(PD),
    SD = sd(PD),
    SEM = sd(PD) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 
str(field_16s_pd_statistics)

field_16s_pd_statistics$Water = factor(field_16s_pd_statistics$Water)
field_16s_pd_statistics$Site = factor(field_16s_pd_statistics$Site)
field_16s_pd_statistics$Type = factor(field_16s_pd_statistics$Type)
field_16s_pd_statistics$Water_Site = factor(field_16s_pd_statistics$Water_Site)
str(field_16s_pd_statistics) 


## stats
library(car)

shapiro.test(field_16s_richness_statistics$Mean) 
shapiro.test(field_16s_pd_statistics$Mean) 

leveneTest(field_16s_richness_statistics$Mean ~ field_16s_richness_statistics$Type, data = field_16s_richness_statistics) 
leveneTest(field_16s_pd_statistics$Mean ~ field_16s_pd_statistics$Type, data = field_16s_pd_statistics) 

t.test(subset(field_16s_richness_statistics, Type == "F")$Mean, subset(field_16s_richness_statistics, Type == "S")$Mean, paired = TRUE) 
t.test(subset(field_16s_pd_statistics, Type == "F")$Mean, subset(field_16s_pd_statistics, Type == "S")$Mean, paired = TRUE) #

wilcox.test(field_16s_richness_statistics$Mean[field_16s_richness_statistics$Type == "F"], 
            field_16s_richness_statistics$Mean[field_16s_richness_statistics$Type == "S"], paired = TRUE) 

wilcox.test(field_16s_pd_statistics$Mean[field_16s_pd_statistics$Type == "F"], 
            field_16s_pd_statistics$Mean[field_16s_pd_statistics$Type == "S"], paired = TRUE) 

## plot
library(ggsignif)

field_16s_richness_figure = ggplot(field_16s_richness_statistics, aes(x=Type, y=Mean)) +
  geom_point(size = 3, aes(color = factor(Type), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color = factor(Type)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("F", "S")), annotations="*", y_position = 2300,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 2500), breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
  ggtitle("") + xlab("") + ylab("Prokaryotic richness") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "Helvetica", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_richness_figure.pdf", plot = field_16s_richness_figure, width = 4, height = 5, units = "cm", device = "pdf") 

field_16s_pd_figure = ggplot(field_16s_pd_statistics, aes(x=Type, y=Mean)) +
  geom_point(size = 3, aes(color = factor(Type), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color = factor(Type)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("F", "S")), annotations="*", y_position = 170,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 180), breaks = c(0, 60, 120, 180)) +
  ggtitle("") + xlab("") + ylab("Prokaryotic\nphylogenetic diversity") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "Helvetica", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_pd_figure.pdf", plot = field_16s_pd_figure, width = 4, height = 5, units = "cm", device = "pdf") 





















### field microeukaryotic diveristy ------------------------------------------------------------------------------------------------------------------------------------------------------------------
## data
library(vegan)
library(picante)

field_18s_richness = estimateR(asv_18s_clean_field)[1,]
field_18s_richness = data.frame(Sample = names(field_18s_richness), Richness = field_18s_richness)
field_18s_richness = left_join(field_metadata, field_18s_richness, by = "Sample")
str(field_18s_richness)
field_18s_richness = as.data.frame(field_18s_richness) 
str(field_18s_richness)

tree_rooted_18sv9 = read.tree(file = "Sequence data/18sv9-rooted-tree.nwk")
field_18s_pd = pd(asv_18s_clean_field_rarefy, tree_rooted_18sv9)
field_18s_pd = data.frame(Sample = rownames(field_18s_pd), PD = field_18s_pd$PD)
field_18s_pd = left_join(field_metadata, field_18s_pd, by = "Sample")
str(field_18s_pd)
field_18s_pd = as.data.frame(field_18s_pd) 
str(field_18s_pd)

field_18s_richness_statistics = field_18s_richness %>% group_by(Water, Site, Type) %>%
  summarise(
    Mean = mean(Richness),
    SD = sd(Richness),
    SEM = sd(Richness) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 
str(field_18s_richness_statistics)

field_18s_richness_statistics$Water = factor(field_18s_richness_statistics$Water)
field_18s_richness_statistics$Site = factor(field_18s_richness_statistics$Site)
field_18s_richness_statistics$Type = factor(field_18s_richness_statistics$Type)
field_18s_richness_statistics$Water_Site = factor(field_18s_richness_statistics$Water_Site)
str(field_18s_richness_statistics)

field_18s_pd_statistics = field_18s_pd %>% group_by(Water, Site, Type) %>%
  summarise(
    Mean = mean(PD),
    SD = sd(PD),
    SEM = sd(PD) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 
str(field_18s_pd_statistics)

field_18s_pd_statistics$Water = factor(field_18s_pd_statistics$Water)
field_18s_pd_statistics$Site = factor(field_18s_pd_statistics$Site)
field_18s_pd_statistics$Type = factor(field_18s_pd_statistics$Type)
field_18s_pd_statistics$Water_Site = factor(field_18s_pd_statistics$Water_Site)
str(field_18s_pd_statistics) 



## stats
library(car)

shapiro.test(field_18s_richness_statistics$Mean) 
shapiro.test(field_18s_pd_statistics$Mean) 

leveneTest(field_18s_richness_statistics$Mean ~ field_18s_richness_statistics$Type, data = field_18s_richness_statistics) 
leveneTest(field_18s_pd_statistics$Mean ~ field_18s_pd_statistics$Type, data = field_18s_pd_statistics)

t.test(subset(field_18s_richness_statistics, Type == "F")$Mean, subset(field_18s_richness_statistics, Type == "S")$Mean, paired = TRUE) 
t.test(subset(field_18s_pd_statistics, Type == "F")$Mean, subset(field_18s_pd_statistics, Type == "S")$Mean, paired = TRUE) 

wilcox.test(field_18s_richness_statistics$Mean[field_18s_richness_statistics$Type == "F"], 
            field_18s_richness_statistics$Mean[field_18s_richness_statistics$Type == "S"], paired = TRUE) 

wilcox.test(field_18s_pd_statistics$Mean[field_18s_pd_statistics$Type == "F"], 
            field_18s_pd_statistics$Mean[field_18s_pd_statistics$Type == "S"], paired = TRUE) 



## plot
library(tidyverse)
library(ggsignif)

field_18s_richness_figure = ggplot(field_18s_richness_statistics, aes(x=Type, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Type), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color = factor(Type)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("F", "S")), annotations="*", y_position = 2300,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 2500), breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
  ggtitle("") + xlab("") + ylab("Microeukaryotic richness") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "Helvetica", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_richness_figure.png", plot = field_18s_richness_figure, width = 5, height = 5, units = "cm", dpi = 600)


field_18s_pd_figure = ggplot(field_18s_pd_statistics, aes(x=Type, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Type), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color = factor(Type)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("F", "S")), annotations="*", y_position = 170,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 180), breaks = c(0, 60, 120, 180)) +
  ggtitle("") + xlab("") + ylab("Microeukaryotic\nphylogenetic diversity") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color = "black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "Helvetica", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_pd_figure.pdf", plot = field_18s_pd_figure, width = 4, height = 5, units = "cm", device = "pdf") 
























### Figure-field conc. and diveristy   ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)

figure_field_conc_alpha = field_16s_richness_figure + field_16s_pd_figure + field_18s_richness_figure + field_18s_pd_figure + field_qubit_figure + field_qpcr_16s_figure + field_qpcr_18s_figure + plot_layout(ncol =4, nrow = 2, byrow = T)
ggsave("figure_field_conc_alpha.pdf", plot = figure_field_conc_alpha, width = 18, height =11, units = "cm", device = "pdf") 
























### field 16s beta ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## data
library(tidyverse)
library(vegan)
library(betapart)

field_16s_jaccard = vegdist(asv_16s_clean_field_rarefy, method = "jaccard")
print(field_16s_jaccard)

field_16s_jaccard_between_af_ps = c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "jaccard"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "jaccard"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "jaccard"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "jaccard"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "jaccard"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "jaccard"))))])
field_16s_jaccard_between_af_ps = data.frame(value = field_16s_jaccard_between_af_ps, group = "pro")
summary(field_16s_jaccard_between_af_ps)

field_16s_turnover_between_af_ps = c(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jtu)))])
field_16s_turnover_between_af_ps = data.frame(value = field_16s_turnover_between_af_ps, group = "pro")
summary(field_16s_turnover_between_af_ps)

field_16s_nestedness_between_af_ps = c(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jne)))])
field_16s_nestedness_between_af_ps = data.frame(value = field_16s_nestedness_between_af_ps, group = "pro")
summary(field_16s_nestedness_between_af_ps)

field_16s_jaccard_nmds = metaMDS(field_16s_jaccard, k = 2)
field_16s_jaccard_nmds_table = data.frame(NMDS1 = field_16s_jaccard_nmds$points[,1], NMDS2 = field_16s_jaccard_nmds$points[,2])
field_16s_jaccard_nmds_table$Group = as.factor(field_metadata$Type)
levels(field_16s_jaccard_nmds_table$Group) <- c("F", "S")
str(field_16s_jaccard_nmds_table)


## stats
library(vegan)
library(betapart)

betadisper(field_16s_jaccard, field_metadata$Type, type = "median")
plot(betadisper(field_16s_jaccard, field_metadata$Type, type = "median"))
anova(betadisper(field_16s_jaccard, field_metadata$Type, type = "median")) 
TukeyHSD(betadisper(field_16s_jaccard, field_metadata$Type, type = "median"))
boxplot(betadisper(field_16s_jaccard, field_metadata$Type, type = "median"))

adonis2(field_16s_jaccard ~ field_metadata$Type, permutations = 999) #
adonis2(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C"),], method = "jaccard") ~ field_metadata$Type[field_metadata$Water == "C"], permutations = 999) 
adonis2(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T"),], method = "jaccard") ~ field_metadata$Type[field_metadata$Water == "T"], permutations = 999) 


## plot 
library(tidyverse)

field_16s_jaccard_nmds_figure = ggplot(field_16s_jaccard_nmds_table, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, aes(color = factor(Group), shape = factor(Group)), fill = NA, alpha = 0.5, show.legend = F) +
  stat_ellipse(aes(color = factor(Group), fill = factor(Group)), geom = "polygon", level = 0.85, alpha = 0.2, show.legend = F) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_fill_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 16)) + 
  scale_x_continuous(limits = c(-8, 8), breaks = c(-8, -4, 0, 4, 8)) + 
  scale_y_continuous(limits = c(-6, 6), breaks = c(-6, -3, 0, 3, 6)) +
  ggtitle("") + xlab("NMDS1") + ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_blank()) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("field_16s_jaccard_nmds_figure.pdf", plot = field_16s_jaccard_nmds_figure, width = 8, height = 6, units = "cm", device = "pdf")

field_16s_turnover_between_af_ps_figure = ggplot(field_16s_turnover_between_af_ps, aes(x = group, y = value)) +
  geom_point(position=position_jitter(width = 0.3), size = 1, color = "#BEB8DC", shape = 16, alpha = 0.8, show.legend = F) +
  geom_boxplot(width=0.5, fill = NA, outlier.shape = NA, size = 0.5, color = "#BEB8DC") +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.5, linetype = "dotted", linewidth = 0.3, color = "#BEB8DC") +
  scale_color_manual(values = "#BEB8DC") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("") + xlab("Active-Passive") + ylab("Turnover") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("field_16s_turnover_between_af_ps_figure.pdf", plot = field_16s_turnover_between_af_ps_figure, width = 3, height = 5, units = "cm", device = "pdf")

field_16s_nestedness_between_af_ps_figure = ggplot(field_16s_nestedness_between_af_ps, aes(x = group, y = value)) +
  geom_point(position=position_jitter(width = 0.3), size = 1, color = "#BEB8DC", shape = 16, alpha = 0.8, show.legend = F) +
  geom_boxplot(width=0.5, fill = NA, outlier.shape = NA, size = 0.5, color = "#BEB8DC") +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.5, linetype = "dotted", linewidth = 0.3, color = "#BEB8DC") +
  scale_color_manual(values = "#BEB8DC") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("") + xlab("Active-Passive") + ylab("Nestedness") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("field_16s_nestedness_between_af_ps_figure.pdf", plot = field_16s_nestedness_between_af_ps_figure, width = 3, height = 5, units = "cm", device = "pdf")

field_16s_jaccard_sum_between_af_ps_figure = ggplot(field_16s_jaccard_sum_between_af_ps, aes(x = group, y = value)) +
  geom_point(position=position_jitter(width = 0.3), size = 1, color = "#BEB8DC", shape = 16, alpha = 0.8, show.legend = F) +
  geom_boxplot(width=0.5, fill = NA, outlier.shape = NA, size = 0.5, color = "#BEB8DC") +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.5, linetype = "dotted", linewidth = 0.3, color = "#BEB8DC") +
  scale_color_manual(values = "#BEB8DC") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("") + xlab("Active-Passive") + ylab("Dissimilarity") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("field_16s_jaccard_sum_between_af_ps_figure.pdf", plot = field_16s_jaccard_sum_between_af_ps_figure, width = 3, height = 5, units = "cm", device = "pdf")
















### field 18s beta ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## data
library(tidyverse)
library(vegan)
library(betapart)

field_18s_jaccard = vegdist(asv_18s_clean_field_rarefy, method = "jaccard")
str(field_18s_jaccard)

field_18s_jaccard_between_af_ps = c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "jaccard"))))],
                                    as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "jaccard"))))],
                                    as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "jaccard"))))],
                                    as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "jaccard"))))],
                                    as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "jaccard"))))],
                                    as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "jaccard"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "jaccard")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "jaccard"))))])
field_18s_jaccard_between_af_ps = data.frame(value = field_18s_jaccard_between_af_ps, group = "pro")
summary(field_18s_jaccard_between_af_ps)

field_18s_turnover_between_af_ps = c(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jtu)))],
                                     as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jtu)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jtu))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jtu)))])
field_18s_turnover_between_af_ps = data.frame(value = field_18s_turnover_between_af_ps, group = "eu")
summary(field_18s_turnover_between_af_ps)

field_18s_nestedness_between_af_ps = c(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "pa"), index.family = "jaccard")$beta.jne)))],
                                       as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jne)[grepl("^F", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jne))), grepl("^S", rownames(as.matrix(beta.pair(decostand(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "pa"), index.family = "jaccard")$beta.jne)))])
field_18s_nestedness_between_af_ps = data.frame(value = field_18s_nestedness_between_af_ps, group = "eu")
summary(field_18s_nestedness_between_af_ps)

field_18s_jaccard_nmds = metaMDS(field_18s_jaccard, k = 2)
field_18s_jaccard_nmds_table = data.frame(NMDS1 = field_18s_jaccard_nmds$points[,1], NMDS2 = field_18s_jaccard_nmds$points[,2])
field_18s_jaccard_nmds_table$Group = as.factor(field_metadata$Type)
levels(field_18s_jaccard_nmds_table$Group) <- c("F", "S")
str(field_18s_jaccard_nmds_table)


## stats
library(vegan)

betadisper(field_18s_jaccard, field_metadata$Type, type = "median")
plot(betadisper(field_18s_jaccard, field_metadata$Type, type = "median"))
anova(betadisper(field_18s_jaccard, field_metadata$Type, type = "median")) 
TukeyHSD(betadisper(field_18s_jaccard, field_metadata$Type, type = "median"))
boxplot(betadisper(field_18s_jaccard, field_metadata$Type, type = "median"))

adonis2(field_18s_jaccard ~ field_metadata$Type, permutations = 999) 
adonis2(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C"),], method = "jaccard") ~ field_metadata$Type[field_metadata$Water == "C"], permutations = 999) 
adonis2(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T"),], method = "jaccard") ~ field_metadata$Type[field_metadata$Water == "T"], permutations = 999) 


## plot 
library(tidyverse)

field_18s_jaccard_nmds_figure = ggplot(field_18s_jaccard_nmds_table, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, aes(color = factor(Group), shape = factor(Group)), fill = NA, alpha = 0.5, show.legend = F) +
  stat_ellipse(aes(color = factor(Group), fill = factor(Group)), geom = "polygon", level = 0.88, alpha = 0.2, show.legend = F) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_fill_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 16)) + 
  scale_x_continuous(limits = c(-6, 6), breaks = c(-6, -3, 0, 3, 6)) +
  scale_y_continuous(limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4)) +
  ggtitle("") + xlab("NMDS1") + ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_blank()) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("field_18s_jaccard_nmds_figure.pdf", plot = field_18s_jaccard_nmds_figure, width = 8, height = 6, units = "cm", device = "pdf")

field_18s_turnover_between_af_ps_figure = ggplot(field_18s_turnover_between_af_ps, aes(x = group, y = value)) +
  geom_point(position=position_jitter(width = 0.3), size = 1, color = "#8ECFC9", shape = 16, alpha = 0.8, show.legend = F) +
  geom_boxplot(width=0.5, fill = NA, outlier.shape = NA, size = 0.5, color = "#8ECFC9") +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.5, linetype = "dotted", linewidth = 0.3, color = "#8ECFC9") +
  scale_color_manual(values = "#8ECFC9") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("") + xlab("Active-Passive") + ylab("Turnover") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("field_18s_turnover_between_af_ps_figure.pdf", plot = field_18s_turnover_between_af_ps_figure, width = 3, height = 5, units = "cm", device = "pdf")


field_18s_nestedness_between_af_ps_figure = ggplot(field_18s_nestedness_between_af_ps, aes(x = group, y = value)) +
  geom_point(position=position_jitter(width = 0.3), size = 1, color = "#8ECFC9", shape = 16, alpha = 0.8, show.legend = F) +
  geom_boxplot(width=0.5, fill = NA, outlier.shape = NA, size = 0.5, color = "#8ECFC9") +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.5, linetype = "dotted", linewidth = 0.3, color = "#8ECFC9") +
  scale_color_manual(values = "#8ECFC9") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("") + xlab("Active-Passive") + ylab("Nestedness") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("field_18s_nestedness_between_af_ps_figure.pdf", plot = field_18s_nestedness_between_af_ps_figure, width = 3, height = 5, units = "cm", device = "pdf")

field_18s_jaccard_sum_between_af_ps_figure = ggplot(field_18s_jaccard_sum_between_af_ps, aes(x = group, y = value)) +
  geom_point(position=position_jitter(width = 0.3), size = 1, color = "#8ECFC9", shape = 16, alpha = 0.8, show.legend = F) +
  geom_boxplot(width=0.5, fill = NA, outlier.shape = NA, size = 0.5, color = "#8ECFC9") +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.5, linetype = "dotted", linewidth = 0.3, color = "#8ECFC9") +
  scale_color_manual(values = "#8ECFC9") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("") + xlab("Active-Passive") + ylab("Dissimilarity") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("field_18s_jaccard_sum_between_af_ps_figure.pdf", plot = field_18s_jaccard_sum_between_af_ps_figure, width = 3, height = 5, units = "cm", device = "pdf")













### Figure-field beta   ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)

figure_field_beta = field_16s_jaccard_nmds_figure + field_16s_jaccard_sum_between_af_ps_figure + field_16s_turnover_between_af_ps_figure + field_16s_nestedness_between_af_ps_figure + 
                    field_18s_jaccard_nmds_figure + field_18s_jaccard_sum_between_af_ps_figure + field_18s_turnover_between_af_ps_figure + field_18s_nestedness_between_af_ps_figure + 
                    plot_layout(ncol =4, nrow = 2, byrow = T, widths = c(8, 3, 3, 3))
ggsave("figure_field_beta.pdf", plot = figure_field_beta, width = 17, height =11, units = "cm", device = "pdf") 






























### field 16s relative abundance  --------------------------------------------------------
library(tidyverse)
# calculate relative proportion out of 100
field_16s_rela_abun = asv_16s_clean_field_rarefy / rowSums(asv_16s_clean_field_rarefy) *100
field_16s_rela_abun  = as.data.frame(field_16s_rela_abun)

rowSums(field_16s_rela_abun) 

# change colnames to taxons
colnames(field_16s_rela_abun) = sapply(colnames(field_16s_rela_abun), function(x) {
  if (x %in% names(setNames(asv_tax_16s$Taxon,asv_tax_16s$`Feature ID`))) {
    return(setNames(asv_tax_16s$Taxon,asv_tax_16s$`Feature ID`)[x])
  } else {
    return(x) 
  }
})

head(colnames(field_16s_rela_abun),100)

# define and apply a function to extract target name
extract_name <- function(x) {
  match <- str_extract(x, "(?<=p__)[^;]+") 
  return(match)
}

colnames(field_16s_rela_abun) = sapply(colnames(field_16s_rela_abun), extract_name)

head(colnames(field_16s_rela_abun),100)

# change NA to Others
colnames(field_16s_rela_abun)[is.na(colnames(field_16s_rela_abun))] = "Others"

head(colnames(field_16s_rela_abun),100)

# make colnames unique
colnames(field_16s_rela_abun) = make.unique(colnames(field_16s_rela_abun))

head(colnames(field_16s_rela_abun),100)

# make table as long form
field_16s_rela_abun_clean <- field_16s_rela_abun %>%
  mutate(Sample = rownames(.)) %>% 
  pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")

# generate unique group column
field_16s_rela_abun_clean <- field_16s_rela_abun_clean %>%
  mutate(Group_unique = gsub("\\..*", "", Group))

# sum columns 
field_16s_rela_abun_clean <- field_16s_rela_abun_clean %>%
  group_by(Sample, Group_unique) %>%
  summarise(Value = sum(Value), .groups = "drop")

# make table as short form
field_16s_rela_abun_clean <- field_16s_rela_abun_clean %>%
  pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

field_16s_rela_abun_clean = as.data.frame(field_16s_rela_abun_clean)
rownames(field_16s_rela_abun_clean) = field_16s_rela_abun_clean$Sample
field_16s_rela_abun_clean = field_16s_rela_abun_clean[, -1] 

# check
rowSums(field_16s_rela_abun_clean) 
colnames(field_16s_rela_abun_clean) 

# make top taxa
field_16s_rela_abun_top = field_16s_rela_abun_clean %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:10] %>% data.frame() 
sum(field_16s_rela_abun_top) # 3478.805 out of 3600
rownames(field_16s_rela_abun_top) 

# make below top taxa as others
field_16s_rela_abun_others = rowSums(field_16s_rela_abun_clean[, !(colnames(field_16s_rela_abun_clean) %in% rownames(field_16s_rela_abun_top))])
sum(field_16s_rela_abun_others)
sum(field_16s_rela_abun_top) + sum(field_16s_rela_abun_others) 

# add existed "Others" in top to Others
#field_16s_rela_abun_others = field_16s_rela_abun_others + field_16s_rela_abun_clean[, "Others"]
#sum(field_16s_rela_abun_others) # 115.2874

# produce top+Others table
#field_16s_rela_abun_table = cbind(field_16s_rela_abun_clean[,rownames(field_16s_rela_abun_top)] %>% select(-"Others"), Others = field_16s_rela_abun_others)
field_16s_rela_abun_table = cbind(field_16s_rela_abun_clean[,rownames(field_16s_rela_abun_top)], Others = field_16s_rela_abun_others)
rowSums(field_16s_rela_abun_table) # 100
sum(field_16s_rela_abun_table) # 3600

# produce longer table-active filtration
field_16s_rela_abun_table_af_longer = field_16s_rela_abun_table %>% 
  filter(str_detect(rownames(.), "F-")) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Group") %>% 
  as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, 
                             "F.T1.1" ~ "ES1-1", "F.T1.2" ~ "ES1-2", "F.T1.3" ~ "ES1-3", 
                             "F.T2.1" ~ "ES2-1", "F.T2.2" ~ "ES2-2", "F.T2.3" ~ "ES2-3",
                             "F.T3.1" ~ "ES3-1", "F.T3.2" ~ "ES3-2", "F.T3.3" ~ "ES3-3",
                             "F.C1.1" ~ "CS1-1", "F.C1.2" ~ "CS1-2", "F.C1.3" ~ "CS1-3",
                             "F.C2.1" ~ "CS2-1", "F.C2.2" ~ "CS2-2", "F.C2.3" ~ "CS2-3",
                             "F.C3.1" ~ "CS3-1", "F.C3.2" ~ "CS3-2", "F.C3.3" ~ "CS3-3")) %>%
  mutate(Method = "AF") %>%
  mutate(Sample = factor(Sample, levels = c("ES1-1", "ES1-2", "ES1-3", 
                                            "ES2-1", "ES2-2", "ES2-3", 
                                            "ES3-1", "ES3-2", "ES3-3", 
                                            "CS1-1", "CS1-2", "CS1-3",
                                            "CS2-1", "CS2-2", "CS2-3",
                                            "CS3-1", "CS3-2", "CS3-3"))) %>%
  arrange(Sample)


str(field_16s_rela_abun_table_af_longer)

field_16s_rela_abun_table_af_longer$Group = factor(field_16s_rela_abun_table_af_longer$Group,levels = colnames(field_16s_rela_abun_table))
str(field_16s_rela_abun_table_af_longer)

sum(field_16s_rela_abun_table_af_longer[,"Proportion"]) # 1800

# produce longer table-passive sampling
field_16s_rela_abun_table_ps_longer = field_16s_rela_abun_table %>% 
  filter(str_detect(rownames(.), "S-")) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Group") %>% 
  as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, 
                             "S.T1.1" ~ "ES1-1", "S.T1.2" ~ "ES1-2", "S.T1.3" ~ "ES1-3", 
                             "S.T2.1" ~ "ES2-1", "S.T2.2" ~ "ES2-2", "S.T2.3" ~ "ES2-3",
                             "S.T3.1" ~ "ES3-1", "S.T3.2" ~ "ES3-2", "S.T3.3" ~ "ES3-3",
                             "S.C1.1" ~ "CS1-1", "S.C1.2" ~ "CS1-2", "S.C1.3" ~ "CS1-3",
                             "S.C2.1" ~ "CS2-1", "S.C2.2" ~ "CS2-2", "S.C2.3" ~ "CS2-3",
                             "S.C3.1" ~ "CS3-1", "S.C3.2" ~ "CS3-2", "S.C3.3" ~ "CS3-3")) %>%
  mutate(Method = "AF") %>%
  mutate(Sample = factor(Sample, levels = c("ES1-1", "ES1-2", "ES1-3", 
                                            "ES2-1", "ES2-2", "ES2-3", 
                                            "ES3-1", "ES3-2", "ES3-3", 
                                            "CS1-1", "CS1-2", "CS1-3",
                                            "CS2-1", "CS2-2", "CS2-3",
                                            "CS3-1", "CS3-2", "CS3-3"))) %>%
  arrange(Sample)

str(field_16s_rela_abun_table_ps_longer)

field_16s_rela_abun_table_ps_longer$Group = factor(field_16s_rela_abun_table_ps_longer$Group,levels = colnames(field_16s_rela_abun_table))
str(field_16s_rela_abun_table_ps_longer)

sum(field_16s_rela_abun_table_ps_longer[,"Proportion"]) # 1800

# plot
field_16s_rela_abun_figure_af = ggplot(field_16s_rela_abun_table_af_longer, aes(x = Sample, y = Proportion, fill = Group)) +
  geom_bar(stat="identity", position = "stack",  width = 0.8, color = "black", linewidth =0.2) +
  scale_fill_manual(values = c("#aec7e8", "#ffbb78", "#C8E6C9", "#ff9896", "#c5b0d5","#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5","#e0e0e080")) +
  ggtitle("")+ xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 6)) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain",angle = 90, vjust = 0.5)) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_rela_abun_figure_af.pdf", plot = field_16s_rela_abun_figure_af, width = 8, height = 6, units = "cm", device = "pdf") 

field_16s_rela_abun_figure_ps = ggplot(field_16s_rela_abun_table_ps_longer, aes(x = Sample, y = Proportion, fill = Group)) +
  geom_bar(stat="identity", position = "stack",  width = 0.8, color = "black", linewidth =0.2) +
  scale_fill_manual(values = c("#aec7e8", "#ffbb78", "#C8E6C9", "#ff9896", "#c5b0d5","#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5","#e0e0e080")) +
  ggtitle("")+ xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 6)) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain",angle = 90, vjust = 0.5)) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_rela_abun_figure_pdf", plot = field_16s_rela_abun_figure_ps, width = 8, height = 6, units = "cm", evice = "pdf") 





















### field 18s relative abundance  --------------------------------------------------------
library(tidyverse)
# calculate relative proportion out of 100
field_18s_rela_abun = asv_18s_clean_field_rarefy / rowSums(asv_18s_clean_field_rarefy) *100
field_18s_rela_abun  = as.data.frame(field_18s_rela_abun)

rowSums(field_18s_rela_abun) 

# change colnames to taxons
colnames(field_18s_rela_abun) = sapply(colnames(field_18s_rela_abun), function(x) {
  if (x %in% names(setNames(asv_tax_18s$Taxon,asv_tax_18s$`Feature ID`))) {
    return(setNames(asv_tax_18s$Taxon,asv_tax_18s$`Feature ID`)[x])
  } else {
    return(x) 
  }
})

head(colnames(field_18s_rela_abun),100)

# define and apply a function to extract target name
extract_name <- function(x) {
  match <- str_extract(x, "(?<=dvs__)[^;]+") 
  return(match)
}

colnames(field_18s_rela_abun) = sapply(colnames(field_18s_rela_abun), extract_name)

head(colnames(field_18s_rela_abun),100)

# change NA to Others
colnames(field_18s_rela_abun)[is.na(colnames(field_18s_rela_abun))] <- "Others"

head(colnames(field_18s_rela_abun),100)

# make colnames unique
colnames(field_18s_rela_abun) = make.unique(colnames(field_18s_rela_abun))

head(colnames(field_18s_rela_abun),100)

# make table as longer form
field_18s_rela_abun_clean = field_18s_rela_abun %>%
  mutate(Sample = rownames(.)) %>% 
  pivot_longer(cols = -Sample, names_to = "Group", values_to = "Value")

# generate unique group column
field_18s_rela_abun_clean = field_18s_rela_abun_clean %>%
  mutate(Group_unique = gsub("\\..*", "", Group))

# sum columns 
field_18s_rela_abun_clean = field_18s_rela_abun_clean %>%
  group_by(Sample, Group_unique) %>%
  summarise(Value = sum(Value), .groups = "drop")

# make table as short form
field_18s_rela_abun_clean = field_18s_rela_abun_clean %>%
  pivot_wider(names_from = Group_unique, values_from = Value, values_fill = 0)

field_18s_rela_abun_clean = as.data.frame(field_18s_rela_abun_clean)
rownames(field_18s_rela_abun_clean) = field_18s_rela_abun_clean$Sample
field_18s_rela_abun_clean = field_18s_rela_abun_clean[, -1] 

# check
rowSums(field_18s_rela_abun_clean)
colnames(field_18s_rela_abun_clean) 

# make top taxa
field_18s_rela_abun_top = field_18s_rela_abun_clean %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:11] %>% data.frame() 
sum(field_18s_rela_abun_top) 
rownames(field_18s_rela_abun_top) 

# make below top taxa as others
field_18s_rela_abun_others = rowSums(field_18s_rela_abun_clean[, !(colnames(field_18s_rela_abun_clean) %in% rownames(field_18s_rela_abun_top))])
sum(field_18s_rela_abun_others) 
sum(field_18s_rela_abun_top) + sum(field_18s_rela_abun_others) 

# add existed "Others" in top to Others
field_18s_rela_abun_others = field_18s_rela_abun_others + field_18s_rela_abun_clean[, "Others"]
sum(field_18s_rela_abun_others) 

# produce top+Others table
field_18s_rela_abun_table = cbind(field_18s_rela_abun_clean[,rownames(field_18s_rela_abun_top)] %>% select(-"Others"), Others = field_18s_rela_abun_others)
rowSums(field_18s_rela_abun_table) 
sum(field_18s_rela_abun_table) 

# produce longer table-active filtration
field_18s_rela_abun_table_af_longer = field_18s_rela_abun_table %>% 
  filter(str_detect(rownames(.), "F-")) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Group") %>% 
  as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, 
                             "F.T1.1" ~ "ES1-1", "F.T1.2" ~ "ES1-2", "F.T1.3" ~ "ES1-3", 
                             "F.T2.1" ~ "ES2-1", "F.T2.2" ~ "ES2-2", "F.T2.3" ~ "ES2-3",
                             "F.T3.1" ~ "ES3-1", "F.T3.2" ~ "ES3-2", "F.T3.3" ~ "ES3-3",
                             "F.C1.1" ~ "CS1-1", "F.C1.2" ~ "CS1-2", "F.C1.3" ~ "CS1-3",
                             "F.C2.1" ~ "CS2-1", "F.C2.2" ~ "CS2-2", "F.C2.3" ~ "CS2-3",
                             "F.C3.1" ~ "CS3-1", "F.C3.2" ~ "CS3-2", "F.C3.3" ~ "CS3-3")) %>%
  mutate(Method = "AF") %>%
  mutate(Sample = factor(Sample, levels = c("ES1-1", "ES1-2", "ES1-3", 
                                            "ES2-1", "ES2-2", "ES2-3", 
                                            "ES3-1", "ES3-2", "ES3-3", 
                                            "CS1-1", "CS1-2", "CS1-3",
                                            "CS2-1", "CS2-2", "CS2-3",
                                            "CS3-1", "CS3-2", "CS3-3"))) %>%
  arrange(Sample)

str(field_18s_rela_abun_table_af_longer)

field_18s_rela_abun_table_af_longer$Group = factor(field_18s_rela_abun_table_af_longer$Group,levels = colnames(field_18s_rela_abun_table))
str(field_18s_rela_abun_table_af_longer)

sum(field_18s_rela_abun_table_af_longer[,"Proportion"]) # 1800

# produce longer table-passive sampling
field_18s_rela_abun_table_ps_longer = field_18s_rela_abun_table %>% 
  filter(str_detect(rownames(.), "S-")) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Group") %>% 
  as_tibble() %>% 
  pivot_longer(!Group, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = case_match(Sample, 
                             "S.T1.1" ~ "ES1-1", "S.T1.2" ~ "ES1-2", "S.T1.3" ~ "ES1-3", 
                             "S.T2.1" ~ "ES2-1", "S.T2.2" ~ "ES2-2", "S.T2.3" ~ "ES2-3",
                             "S.T3.1" ~ "ES3-1", "S.T3.2" ~ "ES3-2", "S.T3.3" ~ "ES3-3",
                             "S.C1.1" ~ "CS1-1", "S.C1.2" ~ "CS1-2", "S.C1.3" ~ "CS1-3",
                             "S.C2.1" ~ "CS2-1", "S.C2.2" ~ "CS2-2", "S.C2.3" ~ "CS2-3",
                             "S.C3.1" ~ "CS3-1", "S.C3.2" ~ "CS3-2", "S.C3.3" ~ "CS3-3")) %>%
  mutate(Method = "AF") %>%
  mutate(Sample = factor(Sample, levels = c("ES1-1", "ES1-2", "ES1-3", 
                                            "ES2-1", "ES2-2", "ES2-3", 
                                            "ES3-1", "ES3-2", "ES3-3", 
                                            "CS1-1", "CS1-2", "CS1-3",
                                            "CS2-1", "CS2-2", "CS2-3",
                                            "CS3-1", "CS3-2", "CS3-3"))) %>%
  arrange(Sample)

str(field_18s_rela_abun_table_ps_longer)

field_18s_rela_abun_table_ps_longer$Group = factor(field_18s_rela_abun_table_ps_longer$Group,levels = colnames(field_18s_rela_abun_table))
str(field_18s_rela_abun_table_ps_longer)

sum(field_18s_rela_abun_table_ps_longer[,"Proportion"]) # 1800

# plot
field_18s_rela_abun_figure_af = ggplot(field_18s_rela_abun_table_af_longer, aes(x = Sample, y = Proportion, fill = Group)) +
  geom_bar(stat="identity", position = "stack",  width = 0.8, color = "black", linewidth =0.2) +
  scale_fill_manual(values = c("#aec7e8", "#ffbb78", "#C8E6C9", "#ff9896", "#c5b0d5","#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5","#e0e0e080")) +
  ggtitle("")+ xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 6)) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain",angle = 90, vjust = 0.5)) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_rela_abun_figure_af.pdf", plot = field_18s_rela_abun_figure_af, width = 8, height = 6, units = "cm", device = "pdf") 

field_18s_rela_abun_figure_ps = ggplot(field_18s_rela_abun_table_ps_longer, aes(x = Sample, y = Proportion, fill = Group)) +
  geom_bar(stat="identity", position = "stack",  width = 0.8, color = "black", linewidth =0.2) +
  scale_fill_manual(values = c("#aec7e8", "#ffbb78", "#C8E6C9", "#ff9896", "#c5b0d5","#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5","#e0e0e080")) +
  ggtitle("")+ xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 6)) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain",angle = 90, vjust = 0.5)) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_rela_abun_figure_ps.pdf", plot = field_18s_rela_abun_figure_ps, width = 8, height = 6, units = "cm", device = "pdf") 


legned_16s = ggplot(field_16s_rela_abun_table_af_longer, aes(x = Sample, y = Proportion, fill = Group)) +
  geom_bar(stat="identity", position = "stack",  width = 0.4, color = "black", linewidth =0.2) +
  scale_fill_manual(values = c("#aec7e8", "#ffbb78", "#C8E6C9", "#ff9896", "#c5b0d5","#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5","#e0e0e080")) 

ggsave("legned_16s.pdf", plot = legned_16s, width = 6, height = 10, units = "cm", device = "pdf") 














## SIMPER ------------------------------------------------------------------------------------------------------------------------------------------------------------
## 16s
library(phyloseq)
library(vegan)
library(tidyverse)
field_metadata = field_metadata[,-1]
phyloseq_16s_field = phyloseq(otu_table(t(asv_16s_clean_field_rarefy), taxa_are_rows = TRUE), tax_table(tax_16s_clean_field_rarefy), sample_data(field_metadata))

# Phylum
phyloseq_16s_field_phylum = tax_glom(phyloseq_16s_field, taxrank = "Phylum") 
taxa_names(phyloseq_16s_field_phylum) = tax_table(phyloseq_16s_field_phylum)[, "Phylum"]
summary(simper(as.data.frame(t(otu_table(phyloseq_16s_field_phylum))), field_metadata$Type, permutations = 999))

simper_16s_phylum = data.frame(taxa = c("Pseudomonadota", "Cyanobacteriota", "Bacteroidota", "Actinomycetota", "Thermoplasmatota"),
                               contribution = c("0.301", "0.214", "0.182", "0.069", "0.069"))
simper_16s_phylum$contribution_percent = as.numeric(simper_16s_phylum$contribution) * 100
simper_16s_phylum$taxa = factor(simper_16s_phylum$taxa, levels = rev(simper_16s_phylum$taxa))

# Class
phyloseq_16s_field_class = tax_glom(phyloseq_16s_field, taxrank = "Class") 
taxa_names(phyloseq_16s_field_class) = tax_table(phyloseq_16s_field_class)[, "Class"]
summary(simper(as.data.frame(t(otu_table(phyloseq_16s_field_class))), field_metadata$Type, permutations = 999))

simper_16s_class = data.frame(taxa = c("Gammaproteobacteria", "Alphaproteobacteria", "Bacteroidia", "Cyanobacteriia", "Poseidoniia"),
                              contribution = c("0.357", "0.142", "0.133", "0.124", "0.051"))
simper_16s_class$contribution_percent = as.numeric(simper_16s_class$contribution) * 100
simper_16s_class$taxa = factor(simper_16s_class$taxa, levels = rev(simper_16s_class$taxa))

# Order
phyloseq_16s_field_order = tax_glom(phyloseq_16s_field, taxrank = "Order") 
taxa_names(phyloseq_16s_field_order) = tax_table(phyloseq_16s_field_order)[, "Order"]
summary(simper(as.data.frame(t(otu_table(phyloseq_16s_field_order))), field_metadata$Type, permutations = 999))

simper_16s_order = data.frame(taxa = c("Pseudomonadales", "Cytophagales", "Cyanobacteriales", "Flavobacteriales", "PCC-6307"),
                              contribution = c("0.293", "0.124", "0.052", "0.051", "0.043"))
simper_16s_order$contribution_percent = as.numeric(simper_16s_order$contribution) * 100
simper_16s_order$taxa = factor(simper_16s_order$taxa, levels = rev(simper_16s_order$taxa))

# Family 
phyloseq_16s_field_family = tax_glom(phyloseq_16s_field, taxrank = "Family") 
taxa_names(phyloseq_16s_field_family) = tax_table(phyloseq_16s_field_family)[, "Family"]
summary(simper(as.data.frame(t(otu_table(phyloseq_16s_field_family))), field_metadata$Type, permutations = 999))

simper_16s_family = data.frame(taxa = c("Cellvibrionaceae", "Cyclobacteriaceae", "Coleofasciculaceae", "Cyanobiaceae", "Poseidoniaceae"),
                               contribution = c("0.287", "0.115", "0.049", "0.04", "0.032"))
simper_16s_family$contribution_percent = as.numeric(simper_16s_family$contribution) * 100
simper_16s_family$taxa = factor(simper_16s_family$taxa, levels = rev(simper_16s_family$taxa))



## 18s
library(phyloseq)
library(vegan)
library(tidyverse)
phyloseq_18s_field = phyloseq(otu_table(t(asv_18s_clean_field_rarefy), taxa_are_rows = TRUE), tax_table(tax_18s_clean_field_rarefy), sample_data(field_metadata))

# Subdivision
phyloseq_18s_field_subdivision = tax_glom(phyloseq_18s_field, taxrank = "Subdivision") 
taxa_names(phyloseq_18s_field_subdivision) = tax_table(phyloseq_18s_field_subdivision)[, "Subdivision"]
summary(simper(as.data.frame(t(otu_table(phyloseq_18s_field_subdivision))), field_metadata$Type, permutations = 999))

simper_18s_subdivision = data.frame(taxa = c("Dinoflagellata", "Gyrista", "Ciliophora", "Bigyra", "Chlorophyta"),
                                    contribution = c("0.242", "0.214", "0.196", "0.108", "0.058"))
simper_18s_subdivision$contribution_percent = as.numeric(simper_18s_subdivision$contribution) * 100
simper_18s_subdivision$taxa = factor(simper_18s_subdivision$taxa, levels = rev(simper_18s_subdivision$taxa))

# Class
phyloseq_18s_field_class = tax_glom(phyloseq_18s_field, taxrank = "Class") 
taxa_names(phyloseq_18s_field_class) = tax_table(phyloseq_18s_field_class)[, "Class"]
summary(simper(as.data.frame(t(otu_table(phyloseq_18s_field_class))), field_metadata$Type, permutations = 999))

simper_18s_class = data.frame(taxa = c("Mediophyceae", "Bacillariophyceae", "Dinophyceae", "Phyllopharyngea", "Spirotrichea"),
                              contribution = c("0.135", "0.132", "0.125", "0.105", "0.091"))
simper_18s_class$contribution_percent = as.numeric(simper_18s_class$contribution) * 100
simper_18s_class$taxa = factor(simper_18s_class$taxa, levels = rev(simper_18s_class$taxa))

# Order
phyloseq_18s_field_order = tax_glom(phyloseq_18s_field, taxrank = "Order") 
taxa_names(phyloseq_18s_field_order) = tax_table(phyloseq_18s_field_order)[, "Order"]
summary(simper(as.data.frame(t(otu_table(phyloseq_18s_field_order))), field_metadata$Type, permutations = 999))

simper_18s_order = data.frame(taxa = c("Thalassiosirales", "Suctoria", "Nanomonadea", "Gymnodiniales", "Oligotrichida"),
                              contribution = c("0.119", "0.056", "0.125", "0.055", "0.041"))
simper_18s_order$contribution_percent = as.numeric(simper_18s_order$contribution) * 100
simper_18s_order$taxa = factor(simper_18s_order$taxa, evels = rev(simper_18s_order$taxa))

# Family 
phyloseq_18s_field_family = tax_glom(phyloseq_18s_field, taxrank = "Family") 
taxa_names(phyloseq_18s_field_family) = tax_table(phyloseq_18s_field_family)[, "Family"]
summary(simper(as.data.frame(t(otu_table(phyloseq_18s_field_family))), field_metadata$Type, permutations = 999))

simper_18s_family = data.frame(taxa = c("Thalassiosiraceae", "Acinetidae", "MAST-3", "Bacillariaceae", "Stephanodiscaceae"),
                               contribution = c("0.076", "0.073", "0.056", "0.038", "0.038"))
simper_18s_family$contribution_percent = as.numeric(simper_18s_family$contribution) * 100
simper_18s_family$taxa = factor(simper_18s_family$taxa, levels = rev(simper_18s_family$taxa))


# Plot
simper_16s_phylum_plot= ggplot(simper_16s_phylum, aes(x = contribution_percent, y = taxa)) +
  geom_bar(stat="identity", width = 0.5, fill = "#BEB8DC", color = "#BEB8DC", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  ggtitle("")+ xlab("Contribution (%)") + ylab("Phylum") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("simper_16s_phylum_plot.pdf", plot = simper_16s_phylum_plot, width = 5, height = 5, units = "cm", device = "pdf")


simper_16s_class_plot= ggplot(simper_16s_class, aes(x = contribution_percent, y = taxa)) +
  geom_bar(stat="identity", width = 0.5, fill = "#BEB8DC", color = "#BEB8DC", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  ggtitle("")+ xlab("Contribution (%)") + ylab("Class") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("simper_16s_class_plot.pdf", plot = simper_16s_class_plot, width = 5, height = 5, units = "cm", device = "pdf")

simper_16s_order_plot= ggplot(simper_16s_order, aes(x = contribution_percent, y = taxa)) +
  geom_bar(stat="identity", width = 0.5, fill = "#BEB8DC", color = "#BEB8DC", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  ggtitle("")+ xlab("Contribution (%)") + ylab("Order") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("simper_16s_order_plot.pdf", plot = simper_16s_order_plot, width = 5, height = 5, units = "cm", device = "pdf")

simper_16s_family_plot= ggplot(simper_16s_family, aes(x = contribution_percent, y = taxa)) +
  geom_bar(stat="identity", width = 0.5, fill = "#BEB8DC", color = "#BEB8DC", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  ggtitle("")+ xlab("Contribution (%)") + ylab("Family") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("simper_16s_family_plot.pdf", plot = simper_16s_family_plot, width = 5, height = 5, units = "cm", device = "pdf")



simper_18s_subdivision_plot= ggplot(simper_18s_subdivision, aes(x = contribution_percent, y = taxa)) +
  geom_bar(stat="identity", width = 0.5, fill = "#8ECFC9", color = "#8ECFC9", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) +
  ggtitle("")+ xlab("Contribution (%)") + ylab("Subdivision") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("simper_18s_subdivision_plot.pdf", plot = simper_18s_subdivision_plot, width = 5, height = 5, units = "cm", device = "pdf")

simper_18s_class_plot= ggplot(simper_18s_class, aes(x = contribution_percent, y = taxa)) +
  geom_bar(stat="identity", width = 0.5, fill = "#8ECFC9", color = "#8ECFC9", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) +
  ggtitle("")+ xlab("Contribution (%)") + ylab("Class") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("simper_18s_class_plot.pdf", plot = simper_18s_class_plot, width = 5, height = 5, units = "cm", device = "pdf")

simper_18s_order_plot= ggplot(simper_18s_order, aes(x = contribution_percent, y = taxa)) +
  geom_bar(stat="identity", width = 0.5, fill = "#8ECFC9", color = "#8ECFC9", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) +
  ggtitle("")+ xlab("Contribution (%)") + ylab("Order") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("simper_18s_order_plot.pdf", plot = simper_18s_order_plot, width = 5, height = 5, units = "cm", device = "pdf")

simper_18s_family_plot= ggplot(simper_18s_family, aes(x = contribution_percent, y = taxa)) +
  geom_bar(stat="identity", width = 0.5, fill = "#8ECFC9", color = "#8ECFC9", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) +
  ggtitle("")+ xlab("Contribution (%)") + ylab("Family") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill = NA,color="black", linewidth = 1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 8, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "Helvetica", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 6, family = "Helvetica", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=6, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("simper_18s_family_plot.pdf", plot = simper_18s_family_plot, width = 5, height = 5, units = "cm", device = "pdf")




## Figure-relative abundance and simper------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)
figure_16s_relative_abundance_simper = field_16s_rela_abun_figure_af + simper_16s_phylum_plot + simper_16s_class_plot + 
                                       field_16s_rela_abun_figure_ps + simper_16s_order_plot + simper_16s_family_plot +
                                       plot_layout(ncol = 3, nrow = 2, byrow = T, widths = c(10, 3, 3))
ggsave("figure_16s_relative_abundance_simper.pdf", plot = figure_16s_relative_abundance_simper, width = 18, height = 12, units = "cm", device = "pdf")

figure_18s_relative_abundance_simper = field_18s_rela_abun_figure_af + simper_18s_subdivision_plot + simper_18s_class_plot + 
                                       field_18s_rela_abun_figure_ps + simper_18s_order_plot + simper_18s_family_plot +
                                       plot_layout(ncol = 3, nrow = 2, byrow = T, widths = c(10, 3, 3))
ggsave("figure_18s_relative_abundance_simper.pdf", plot = figure_18s_relative_abundance_simper, width = 18, height = 12, units = "cm", device = "pdf")































### field 16s 18s  mantel -------------------------------------------------------------------
library(vegan)
library(tidyverse)

# average 3 replicates 
asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy), "F-C1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-C1") %>%
  bind_rows(asv_16s_clean_field_rarefy, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "F-C2")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-C2") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "F-C3")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-C3") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "F-T1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-T1") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "F-T2")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-T2") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "F-T3")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-T3") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)


asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "S-C1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-C1") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "S-C2")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-C2") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "S-C3")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-C3") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "S-T1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-T1") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "S-T2")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-T2") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_16s_clean_field_rarefy_mean), "S-T3")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-T3") %>%
  bind_rows(asv_16s_clean_field_rarefy_mean, .)

asv_16s_clean_field_rarefy_mean = tail(asv_16s_clean_field_rarefy_mean,12)
rownames(asv_16s_clean_field_rarefy_mean) = asv_16s_clean_field_rarefy_mean$Sample
asv_16s_clean_field_rarefy_mean = asv_16s_clean_field_rarefy_mean %>%
  select(-Sample)

# 18s
# average 3 replicates 
asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy), "F-C1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-C1") %>%
  bind_rows(asv_18s_clean_field_rarefy, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "F-C2")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-C2") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "F-C3")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-C3") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "F-T1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-T1") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "F-T2")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-T2") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "F-T3")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-T3") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "S-C1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-C1") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "S-C2")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-C2") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "S-C3")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-C3") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "S-T1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-T1") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "S-T2")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-T2") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  filter(str_detect(rownames(asv_18s_clean_field_rarefy_mean), "S-T3")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "S-T3") %>%
  bind_rows(asv_18s_clean_field_rarefy_mean, .)

asv_18s_clean_field_rarefy_mean = tail(asv_18s_clean_field_rarefy_mean,12)
rownames(asv_18s_clean_field_rarefy_mean) = asv_18s_clean_field_rarefy_mean$Sample
asv_18s_clean_field_rarefy_mean = asv_18s_clean_field_rarefy_mean %>%
  select(-Sample)

# environmental factors
env_factor = data.frame(Site = c("C1","C2","C3","T1", "T2", "T3"),
                        Temperature = c(27.77, 28.48, 27.49, 28.71, 27.99, 28.56),
                        Salinity = c(33.93, 34.18,  35.31, 18.78,18.22,18.80),
                        Turbidity = c(2.53, 0.64, 1.30, 57.72, 44.58, 41.96),
                        pH = c(8.09, 8.04, 8.05, 7.96, 8.03, 7.96),
                        DO = c(6.5, 6.5, 6.4, 6.0, 6.4, 6.0),
                        Chlorophyll = c(1.55, 0.77, 1.09, 5.00, 8.67, 4.75))

rownames(env_factor) = env_factor[,1]
env_factor = env_factor[,-1]

# define a function
batch_mantel <- function(species_dist, env_data, permutations = 999) {
  results <- list()
  
  for (i in 1:ncol(env_data)) {
    env_var <- env_data[, i] 
    env_dist <- dist(env_var, method = "euclidean") 
    
    result <- mantel(species_dist, env_dist, method = "pearson", permutations = permutations)
    
    results[[colnames(env_data)[i]]] <- data.frame(
      Variable = colnames(env_data)[i],
      Mantel_r = result$statistic,
      p_value = result$signif
    )
  }
  
  return(do.call(rbind, results))
}

# bray distance
field_16s_bray_af = vegdist(asv_16s_clean_field_rarefy_mean[1:6,], method = "bray")
field_16s_bray_ps = vegdist(asv_16s_clean_field_rarefy_mean[7:12,], method = "bray")

field_18s_bray_af = vegdist(asv_18s_clean_field_rarefy_mean[1:6,], method = "bray")
field_18s_bray_ps = vegdist(asv_18s_clean_field_rarefy_mean[7:12,], method = "bray")

# mantel test
field_16s_mantel_af <- batch_mantel(field_16s_bray_af, env_factor)
field_16s_mantel_ps <- batch_mantel(field_16s_bray_ps, env_factor)

field_18s_mantel_af <- batch_mantel(field_18s_bray_af, env_factor)
field_18s_mantel_ps <- batch_mantel(field_18s_bray_ps, env_factor)


## plot
library(linkET)

field_mantel_af_figure = qcorrplot(correlate(env_factor, method = "pearson"), type = "lower", diag = TRUE) +
  geom_square() + 
  geom_couple(aes(colour = pd, size = rd), data = mantal_combined_plot[1:12,],
              curvature = nice_curvature()) +
  scale_fill_gradientn(
    colours = c("#05B9E2", "#4FC0EC", "#A2A2A288", "#E6B366", "#BB9727"),
    limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_colour_manual(values = c("#F27970", "#A2A2A288")) +
  scale_size_manual(values = c(0.5, 2)) +
  guides(
    colour = guide_legend(title = "Mantel's p", order = 1, override.aes = list(size = 3)),
    size = guide_legend(title = "Mantel's r", order = 2), 
    fill = guide_colorbar(title = "Pearson's r", order = 3, ticks = FALSE)
  ) +
  labs(title = "") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.text.y = element_blank(),
    axis.text.y.right = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("field_mantel_af_figure.pdf", plot = field_mantel_af_figure, width = 8, height = 10, units = "cm", device = "pdf")

field_mantel_ps_figure = qcorrplot(correlate(env_factor, method = "pearson"), type = "lower", diag = TRUE) +
  geom_square() + 
  geom_couple(aes(colour = pd, size = rd), data = mantal_combined_plot[13:24,], 
              curvature = nice_curvature()) +
  scale_fill_gradientn(
    colours = c("#05B9E2", "#4FC0EC", "#F8F8F0", "#E6B366", "#BB9727"),
    limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_colour_manual(values = c("#F27970", "#A2A2A288")) +
  scale_size_manual(values = c(0.5, 2)) +
  guides(
    colour = guide_legend(title = "Mantel's p", order = 1, override.aes = list(size = 3)),
    size = guide_legend(title = "Mantel's r", order = 2), 
    fill = guide_colorbar(title = "Pearson's r", order = 3, ticks = FALSE)
  ) +
  labs(title = "") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.text.y = element_blank(),
    axis.text.y.right = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("field_mantel_ps_figure.pdf", plot = field_mantel_ps_figure, width = 8, height = 8, units = "cm", device = "pdf")
























### field 16s lefse ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(phyloseq)
library(tidyverse)

class (asv_16s_clean_field_rarefy)
class (asv_16s_clean_field_rarefy)

tax_16s_seperate_taxonomy = tax_16s %>%
  separate_rows(Taxon, sep = "; ") %>%
  separate_wider_delim(Taxon, delim = "__", names = c("Rank", "Name"), too_few = "align_start") %>%
  mutate(Rank = str_to_lower(str_sub(Rank, 1, 1)),
         Name = na_if(trimws(Name), "")) %>%
  filter(Rank %in% c("d", "p", "c", "o", "f", "g", "s")) %>%
  group_by(`Feature ID`) %>%
  pivot_wider(names_from = Rank, values_from = Name) %>%
  ungroup() %>%
  right_join(tibble(`Feature ID` = tax_16s$`Feature ID`), by = "Feature ID") %>%
  rename(Domain = d, Phylum = p, Class = c, Order = o, Family = f, Genus = g, Species = s) # separate taxonomy

nrow (tax_16s_seperate_taxonomy)
nrow (tax_16s)

class(tax_16s_seperate_taxonomy)
tax_16s_seperate_taxonomy = as.data.frame(tax_16s_seperate_taxonomy)
class(tax_16s_seperate_taxonomy)
rownames(tax_16s_seperate_taxonomy) = tax_16s_seperate_taxonomy[,"Feature ID"] 
tax_16s_clean_field_rarefy = tax_16s_seperate_taxonomy %>%
  filter(`Feature ID` %in% rownames(t(asv_16s_clean_field_rarefy))) # remove extra taxa

tax_16s_clean_field_rarefy = tax_16s_clean_field_rarefy[,-1]
tax_16s_clean_field_rarefy = tax_16s_clean_field_rarefy[,-1]

class (tax_16s_clean_field_rarefy)
tax_16s_clean_field_rarefy = as.matrix(tax_16s_clean_field_rarefy)
class (tax_16s_clean_field_rarefy)


# prepare sample table-active filtration
field_metadata_for_phyloseq_af = as.data.frame(field_metadata)
rownames(field_metadata_for_phyloseq_af) = field_metadata_for_phyloseq_af[,"Sample"]
field_metadata_for_phyloseq_af = field_metadata_for_phyloseq_af[1:18,c("Type","Water")]
str (field_metadata_for_phyloseq_af)

# prepare sample table-passive sampling
field_metadata_for_phyloseq_ps = as.data.frame(field_metadata)
rownames(field_metadata_for_phyloseq_ps) = field_metadata_for_phyloseq_ps[,"Sample"]
field_metadata_for_phyloseq_ps = field_metadata_for_phyloseq_ps[19:36,c("Type","Water")]
str (field_metadata_for_phyloseq_ps)

# 3 table to 1 phyloseq object
phyloseq_16s_af = phyloseq(otu_table(t(asv_16s_clean_field_rarefy[1:18,]), taxa_are_rows = TRUE), tax_table(tax_16s_clean_field_rarefy), sample_data(field_metadata_for_phyloseq_af))

phyloseq_16s_ps = phyloseq(otu_table(t(asv_16s_clean_field_rarefy[19:36,]), taxa_are_rows = TRUE), tax_table(tax_16s_clean_field_rarefy), sample_data(field_metadata_for_phyloseq_ps))


# stats
library(lefser)
library(mia)

phyloseq_16s_af_summarized = mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseq_16s_af)
colData(phyloseq_16s_af_summarized)

set.seed(1234)
lefse_16s_af = lefser(phyloseq_16s_af_summarized, 
                   groupCol = "Water", 
                   #blockCol = "Type", 
                   kruskal.threshold = 0.05,
                   wilcox.threshold = 0.05,
                   lda.threshold = 2)

lefse_16s_af_table <- lefse_16s_af %>% 
  dplyr::left_join(., (tax_16s_seperate_taxonomy %>%
                         rownames_to_column(var="Names"))) %>%
  mutate(Group = if_else(scores < 0, "C", "T")) %>%
  filter(!is.na(Species)) %>%  
  distinct(Genus, .keep_all = TRUE) %>% 
  filter(Confidence >=0.99 ) %>% 
  #filter(scores < -2.5 | scores > 2.5) %>% 
  arrange(scores)

lefse_16s_af_table = lefse_16s_af_table[-13,]

# passive sampling
phyloseq_16s_ps_summarized = mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseq_16s_ps)
colData(phyloseq_16s_ps_summarized)

set.seed(1234)
lefse_16s_ps = lefser(phyloseq_16s_ps_summarized, 
                      groupCol = "Water", 
                      #blockCol = "Type", 
                      kruskal.threshold = 0.05,
                      wilcox.threshold = 0.05,
                      lda.threshold = 2)

lefse_16s_ps_table <- lefse_16s_ps %>% 
  dplyr::left_join(., (tax_16s_seperate_taxonomy %>%
                         rownames_to_column(var="Names"))) %>%
  mutate(Group = if_else(scores < 0, "C", "T")) %>%
  filter(!is.na(Species)) %>%  
  distinct(Genus, .keep_all = TRUE) %>% 
  filter(Confidence >=0.99 ) %>% 
  #filter(scores < -2.5 | scores > 2.5) %>% 
  arrange(scores)



## plot
library(tidyverse)
lefse_16s_af_figure = ggplot(lefse_16s_ps_table, aes(x=scores, y=reorder(Species, +scores))) +
  geom_col(aes(fill = reorder(Group, -scores), color = reorder(Group, -scores)), width=0.5) +
  scale_fill_manual(values=c("#BB9727", "#05B9E2")) + 
  scale_color_manual(values=c("#BB9727", "#05B9E2")) + 
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  ggtitle("16S_AF")+ xlab("LDA score (log10)") + ylab("") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.6, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 0, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=0, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.3)) +
  theme(axis.ticks.y = element_blank()) 

ggsave("lefse_16s_af_figure.pdf", plot = lefse_16s_af_figure, width = 5, height = 5, units = "cm", device = "pdf")

lefse_16s_ps_figure = ggplot(lefse_16s_af_table, aes(x=scores, y=reorder(Species, +scores))) +
  geom_col(aes(fill = reorder(Group, -scores), color = reorder(Group, -scores)), width=0.5) +
  scale_fill_manual(values=c("#BB9727", "#05B9E2")) + 
  scale_color_manual(values=c("#BB9727", "#05B9E2")) + 
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  ggtitle("16S_PS")+ xlab("LDA score (log10)") + ylab("") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.6, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 0, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=0, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.3)) +
  theme(axis.ticks.y = element_blank()) 

ggsave("lefse_16s_ps_figure.pdf", plot = lefse_16s_ps_figure, width = 5, height = 5, units = "cm", device = "pdf")






















### field 18s lefse ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(phyloseq)
library(tidyverse)

class (asv_18s_clean_field_rarefy)
class (asv_18s_clean_field_rarefy)


tax_18s_seperate_taxonomy = tax_18s %>%
  separate(Taxon, 
           into = c("Domain", "Supergroup", "Division", "Subdivision", 
                    "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";", 
           remove = FALSE,  
           fill = "right")  

tax_18s_seperate_taxonomy = tax_18s_seperate_taxonomy %>%
  mutate(across(c(Domain, Supergroup, Division, Subdivision, 
                  Class, Order, Family, Genus, Species), 
                ~gsub("^[a-z]+__", "", .)))

class(tax_18s_seperate_taxonomy)
tax_18s_seperate_taxonomy = as.data.frame(tax_18s_seperate_taxonomy)
class(tax_18s_seperate_taxonomy)
rownames(tax_18s_seperate_taxonomy) = tax_18s_seperate_taxonomy[,"Feature ID"] 
tax_18s_clean_field_rarefy = tax_18s_seperate_taxonomy %>%
  filter(`Feature ID` %in% rownames(t(asv_18s_clean_field_rarefy))) 

tax_18s_clean_field_rarefy = tax_18s_clean_field_rarefy[,-1]
tax_18s_clean_field_rarefy = tax_18s_clean_field_rarefy[,-1]
tax_18s_clean_field_rarefy = tax_18s_clean_field_rarefy[, -10]

class (tax_18s_clean_field_rarefy)
tax_18s_clean_field_rarefy = as.matrix(tax_18s_clean_field_rarefy)
class (tax_18s_clean_field_rarefy)


# prepare sample table-active filtration
field_metadata_for_phyloseq_af = as.data.frame(field_metadata)
rownames(field_metadata_for_phyloseq_af) = field_metadata_for_phyloseq_af[,"Sample"]
field_metadata_for_phyloseq_af = field_metadata_for_phyloseq_af[1:18,c("Type","Water")]
str (field_metadata_for_phyloseq_af)

# prepare sample table-passive sampling
field_metadata_for_phyloseq_ps = as.data.frame(field_metadata)
rownames(field_metadata_for_phyloseq_ps) = field_metadata_for_phyloseq_ps[,"Sample"]
field_metadata_for_phyloseq_ps = field_metadata_for_phyloseq_ps[19:36,c("Type","Water")]
str (field_metadata_for_phyloseq_ps)

# 3 table to 1 phyloseq object
phyloseq_18s_af = phyloseq(otu_table(t(asv_18s_clean_field_rarefy[1:18,]), taxa_are_rows = TRUE), tax_table(tax_18s_clean_field_rarefy), sample_data(field_metadata_for_phyloseq_af))

phyloseq_18s_ps = phyloseq(otu_table(t(asv_18s_clean_field_rarefy[19:36,]), taxa_are_rows = TRUE), tax_table(tax_18s_clean_field_rarefy), sample_data(field_metadata_for_phyloseq_ps))


# lefser
library(lefser)
library(mia)

phyloseq_18s_af_summarized = mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseq_18s_af)
colData(phyloseq_18s_af_summarized)

set.seed(1234)
lefse_18s_af = lefser(phyloseq_18s_af_summarized, 
                      groupCol = "Water", 
                      #blockCol = "Type", 
                      kruskal.threshold = 0.05,
                      wilcox.threshold = 0.05,
                      lda.threshold = 2)

lefse_18s_af_table = lefse_18s_af %>% 
  dplyr::left_join(., (tax_18s_seperate_taxonomy %>%
                         rownames_to_column(var="Names"))) %>%
  mutate(Group = if_else(scores < 0, "C", "T")) %>%
  filter(!is.na(Species)) %>%  
  distinct(Genus, .keep_all = TRUE) %>% 
  filter(Confidence >=0.99 ) %>% 
  #filter(scores < -2.5 | scores > 2.5) %>% 
  arrange(scores)


phyloseq_18s_ps_summarized = mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseq_18s_ps)
colData(phyloseq_18s_ps_summarized)

set.seed(1234)
lefse_18s_ps = lefser(phyloseq_18s_ps_summarized, 
                      groupCol = "Water", 
                      #blockCol = "Type", 
                      kruskal.threshold = 0.05,
                      wilcox.threshold = 0.05,
                      lda.threshold = 2)

lefse_18s_ps_table = lefse_18s_ps %>% 
  dplyr::left_join(., (tax_18s_seperate_taxonomy %>%
                         rownames_to_column(var="Names"))) %>%
  mutate(Group = if_else(scores < 0, "C", "T")) %>%
  filter(!is.na(Species)) %>%  
  distinct(Genus, .keep_all = TRUE) %>% 
  filter(Confidence >=0.99 ) %>% 
  #filter(scores < -2.5 | scores > 2.5) %>% 
  arrange(scores)

# plot
lefse_18s_af_figure = ggplot(lefse_18s_ps_table, aes(x=scores, y=reorder(Species, +scores))) +
  geom_col(aes(fill = reorder(Group, -scores), color = reorder(Group, -scores)), width=0.5) +
  scale_fill_manual(values=c("#BB9727", "#05B9E2")) + 
  scale_color_manual(values=c("#BB9727", "#05B9E2")) + 
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  ggtitle("18S_AF")+ xlab("LDA score (log10)") + ylab("") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.6, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 0, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=0, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.3)) +
  theme(axis.ticks.y = element_blank()) 

ggsave("lefse_18s_af_figure.pdf", plot = lefse_18s_af_figure, width = 4, height = 6, units = "cm", device = "pdf")

lefse_18s_ps_figure = ggplot(lefse_18s_af_table, aes(x=scores, y=reorder(Species, +scores))) +
  geom_col(aes(fill = reorder(Group, -scores), color = reorder(Group, -scores)), width=0.5) +
  scale_fill_manual(values=c("#BB9727", "#05B9E2")) + 
  scale_color_manual(values=c("#BB9727", "#05B9E2")) + 
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  ggtitle("18S_PS")+ xlab("LDA score (log10)") + ylab("") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.6, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "Helvetica", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 0, family = "Helvetica", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "Helvetica", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=0, family = "Helvetica", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.3)) +
  theme(axis.ticks.y = element_blank())          

ggsave("lefse_18s_ps_figure.pdf", plot = lefse_18s_ps_figure, width = 4, height = 6, units = "cm", device = "pdf")













### Figure-mantel and lefse --------------------------------------------------------------------------------------------------------------
library(patchwork)

figure_field_mantel = field_mantel_af_figure + field_mantel_ps_figure + plot_layout(ncol = 1, nrow = 2, byrow = T)
figure_field_lefse = lefse_16s_af_figure + lefse_18s_af_figure + lefse_16s_ps_figure + lefse_18s_ps_figure + plot_layout(ncol = 2, nrow = 2, byrow = T)

ggsave("figure_field_mantel.pdf", plot = figure_field_mantel, width = 15, height =15, units = "cm", device = "pdf")
ggsave("figure_field_lefse.pdf", plot = figure_field_lefse, width = 10, height =10, units = "cm", device = "pdf")
