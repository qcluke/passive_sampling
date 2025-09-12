setwd("/Users/macofqc/Desktop/Chapter 1/Analysis/Passive Sampling")

### mesocosm qubit ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# import table
library(tidyverse)
mesocosm_qubit = read_csv("Qubit data/Chapter1_Mesocosm_Qubit.csv")
str(mesocosm_qubit)

# clean table
mesocosm_qubit = as.data.frame(mesocosm_qubit) 
mesocosm_qubit$Time = factor (mesocosm_qubit$Time, levels = c("10m","1h","6h","12h","24h"))
mesocosm_qubit$Extraction = factor (mesocosm_qubit$Extraction, levels = c("PW","PS","BT"))
str(mesocosm_qubit)

# Shapiro-Wilk test (normality) 
shapiro.test(mesocosm_qubit$Qubit_concentration)   

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(mesocosm_qubit$Qubit_concentration ~ mesocosm_qubit$Extraction, data = mesocosm_qubit)
leveneTest(mesocosm_qubit$Qubit_concentration ~ mesocosm_qubit$Time, data = mesocosm_qubit) 

# ANOVA + Tukey HSD
library(multcompView)
summary(aov(Qubit_concentration~Time, data=mesocosm_qubit)) # p < 0.001
TukeyHSD(aov(Qubit_concentration~Time, data=mesocosm_qubit))
multcompLetters4(aov(Qubit_concentration~Time, data=mesocosm_qubit), TukeyHSD(aov(Qubit_concentration~Time, data=mesocosm_qubit)))
# 24h 12h  6h  1h 10m 
#"a" "b" "c" "c" "c"

summary(aov(Qubit_concentration~Extraction, data=mesocosm_qubit)) # p < 0.01
TukeyHSD(aov(Qubit_concentration~Extraction, data=mesocosm_qubit))
multcompLetters4(aov(Qubit_concentration~Extraction, data=mesocosm_qubit), TukeyHSD(aov(Qubit_concentration~Extraction, data=mesocosm_qubit)))
# BT  PS  PW  
#"a" "b" "b" 

summary(aov(Qubit_concentration~Time * Extraction, data=mesocosm_qubit)) # time: p < 0.001; extraction: p < 0.001； time:extraction: p < 0.001； 
TukeyHSD(aov(Qubit_concentration~Time * Extraction, data=mesocosm_qubit))
multcompLetters4(aov(Qubit_concentration~Time * Extraction, data=mesocosm_qubit), TukeyHSD(aov(Qubit_concentration~Time * Extraction, data=mesocosm_qubit)))
# 24h 12h  6h  1h 10m 
#"a" "b" "c" "d" "d"

# BT  PS  PW  
#"a" "b" "c"

#24h:BT 12h:BT 24h:PS 24h:PW 12h:PS 12h:PW  6h:BT  6h:PS  6h:PW  1h:BT  1h:PS 10m:BT  1h:PW 10m:PS 10m:PW 
#"a"    "b"   "bc"   "cd"   "de"   "ef"   "ef"  "efg"  "fgh"  "fgh"  "fgh"   "gh"    "h"    "h"    "h" 

# histogram-time and extraction
mesocosm_qubit_statistics = mesocosm_qubit %>%
  group_by(Time, Extraction) %>%
  summarise(
    Conc_Mean = mean(Qubit_concentration),
    Conc_SD = sd(Qubit_concentration),
    Conc_SEM = sd(Qubit_concentration) / sqrt(6),
    .groups = "drop"
  )

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
  geom_text(x="10m", y=30, label="h h gh", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="1h", y=30, label="h fgh fgh", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="6h", y=30, label="fhg efg ef",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="12h", y=30, label="ef de b",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="24h", y=30, label="cd bc a",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x=0.5, y=20, label = "(A)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qubit_figure.png", plot = mesocosm_qubit_figure, width = 8, height = 5, units = "cm", dpi = 600)

# histogram-time
mesocosm_qubit_statistics_time = mesocosm_qubit %>%
  group_by(Time) %>%
  summarise(
    Mean = mean(Qubit_concentration),
    SD = sd(Qubit_concentration),
    SEM = sd(Qubit_concentration) / sqrt(18),
    .groups = "drop"
  )

mesocosm_qubit_figure_time = ggplot(mesocosm_qubit_statistics_time, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", " 24 hours")) +
  #ggtitle("(A)") + 
  #xlab("Submersion time") + 
  ylab("Total eDNA conc.\n(ng/μL)") +
  geom_text(x="10m", y=20, label="d", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="1h", y=20, label="d", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="6h", y=20, label="c",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="12h", y=20, label="b",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="24h", y=20, label="a",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x=0.5, y=20, label = "(A)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qubit_figure_time.png", plot = mesocosm_qubit_figure_time, width = 8, height = 5, units = "cm", dpi = 600)

# histogram-extraction
mesocosm_qubit_statistics_extraction = mesocosm_qubit %>%
  group_by(Extraction) %>%
  summarise(
    Mean = mean(Qubit_concentration),
    SD = sd(Qubit_concentration),
    SEM = sd(Qubit_concentration) / sqrt(30),
    .groups = "drop"
  )

mesocosm_qubit_figure_extraction = ggplot(mesocosm_qubit_statistics_extraction, aes(x = Extraction, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_discrete(label = c("PowerWater", "PowerSoil", "Blood&Tissue")) +
  #ggtitle("(A)") + 
  #xlab("Submersion time") + 
  ylab("Total eDNA conc.\n(ng/μL)") +
  geom_text(x="PW", y=20, label="c", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="PS", y=20, label="b", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="BT", y=20, label="a",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x=0.5, y=20, label = "(D)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 2.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) 

ggsave("mesocosm_qubit_figure_extraction.png", plot = mesocosm_qubit_figure_extraction, width = 8, height = 5, units = "cm", dpi = 600) 





















### mesocosm qpcr 16s ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# import table
library(tidyverse)
mesocosm_qpcr_16s = read_csv("qPCR data/Chapter1_Mesocosm_16S.csv")

# clean data
mesocosm_qpcr_16s = as.data.frame(mesocosm_qpcr_16s) 
mesocosm_qpcr_16s$Time = factor (mesocosm_qpcr_16s$Time, levels = c("10m","1h","6h","12h","24h"))
mesocosm_qpcr_16s$Extraction = factor (mesocosm_qpcr_16s$Extraction, levels = c("PW","PS","BT"))
str(mesocosm_qpcr_16s)

# Shapiro-Wilk test (normality) 
shapiro.test(mesocosm_qpcr_16s$Pro_Concentration) 

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(mesocosm_qpcr_16s$Pro_Concentration ~ mesocosm_qpcr_16s$Extraction, data = mesocosm_qpcr_16s)
leveneTest(mesocosm_qpcr_16s$Pro_Concentration ~ mesocosm_qpcr_16s$Time, data = mesocosm_qpcr_16s) 

# ANOVA + Tukey HSD
library(multcompView)
summary(aov(Pro_Concentration~Time, data=mesocosm_qpcr_16s)) # p <0.001
TukeyHSD(aov(Pro_Concentration~Time, data=mesocosm_qpcr_16s))
multcompLetters4(aov(Pro_Concentration~Time, data=mesocosm_qpcr_16s), TukeyHSD(aov(Pro_Concentration~Time, data=mesocosm_qpcr_16s)))
# 24h 12h  6h  1h 10m 
#"a" "b" "c" "d" "

summary(aov(Pro_Concentration~Extraction, data=mesocosm_qpcr_16s)) # p <0.05
TukeyHSD(aov(Pro_Concentration~Extraction, data=mesocosm_qpcr_16s))
multcompLetters4(aov(Pro_Concentration~Extraction, data=mesocosm_qpcr_16s), TukeyHSD(aov(Pro_Concentration~Extraction, data=mesocosm_qpcr_16s)))
# BT  PS  PW  
#""a" "ab"  "b"

summary(aov(Pro_Concentration~Time * Extraction, data=mesocosm_qpcr_16s)) # time: p < 0.001; extraction: p < 0.001； time and extraction: p < 0.001； 
TukeyHSD(aov(Pro_Concentration~Time * Extraction, data=mesocosm_qpcr_16s))
multcompLetters4(aov(Pro_Concentration~Time * Extraction, data=mesocosm_qpcr_16s), TukeyHSD(aov(Pro_Concentration~Time * Extraction, data=mesocosm_qpcr_16s)))
# 24h 12h  6h  1h 10m 
#"a" "b" "c" "d" "e" 

# BT  PS  PW  
#""a" "b"  "b"

#24h:BT 12h:BT 24h:PS 24h:PW 12h:PW 12h:PS  6h:BT  6h:PW  1h:BT  6h:PS  1h:PS  1h:PW 10m:BT 10m:PS 10m:PW 
#"a"   "b"   "bc"   "cd"   "de"   "de"   "ef"   "fg"    "g"    "g"   "gh"   "gh"    "h"    "h"    "h"

# histogram-time and extraction
mesocosm_qpcr_16s_statistics = mesocosm_qpcr_16s %>% 
  group_by(Time, Extraction) %>%
  summarise(
    Mean = mean(Pro_Concentration),
    SD = sd(Pro_Concentration),
    SEM = sd(Pro_Concentration) / sqrt(6),
    .groups = "drop"
  )

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
  geom_text(x="10m", y=4, label="h h h", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="1h", y=4, label="gh gh g", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="6h", y=4, label="fg g ef ",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="12h", y=4, label="de de b",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="24h", y=4, label="cd bc a",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x=0.5, y=3, label = "(B)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qpcr_16s_figure.png", plot = mesocosm_qpcr_figure, width = 8, height = 5, units = "cm", dpi = 600)

# histogram-time
mesocosm_qpcr_16s_statistics_time = mesocosm_qpcr_16s %>%
  group_by(Time) %>%
  summarise(
    Mean = mean(Pro_Concentration),
    SD = sd(Pro_Concentration),
    SEM = sd(Pro_Concentration) / sqrt(18),
    .groups = "drop"
  )

mesocosm_qpcr_16s_figure_time = ggplot(mesocosm_qpcr_16s_statistics_time, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", " 24 hours")) +
  #ggtitle("(A)") + 
  #xlab("Submersion time") + 
  ylab("Prokaryotic eDNA yield\n(million copies/μL)") +
  geom_text(x="10m", y=4, label="e", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="1h", y=4, label="d", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="6h", y=4, label="c",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="12h", y=4, label="b",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="24h", y=4, label="a",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x=0.5, y=4, label = "(B)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qpcr_16s_figure_time.png", plot = mesocosm_qpcr_figure_time, width = 8, height = 5, units = "cm", dpi = 600)

# histogram-Extraction
mesocosm_qpcr_16s_statistics_extraction <- mesocosm_qpcr_16s %>%
  group_by(Extraction) %>%
  summarise(
    Mean = mean(Pro_Concentration),
    SD = sd(Pro_Concentration),
    SEM = sd(Pro_Concentration) / sqrt(30),
    .groups = "drop"
  )

mesocosm_qpcr_16s_figure_extraction = ggplot(mesocosm_qpcr_16s_statistics_extraction, aes(x = Extraction, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) +
  scale_x_discrete(label = c("PowerWater", "PowerSoil", "Blood&Tissue")) +
  #ggtitle("(A)") + 
  #xlab("Submersion time") + 
  ylab("Prokaryotic eDNA yield\n(million copies/μL)") +
  geom_text(x="PW", y=4, label="b", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="PS", y=4, label="b", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="BT", y=4, label="a",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x=0.5, y=4, label = "(E)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 2.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) 

ggsave("mesocosm_qpcr_16s_figure_time_extraction.png", plot = mesocosm_qpcr_figure_time_extraction, width = 8, height = 5, units = "cm", dpi = 600) 
























### mesocosm qpcr 18s ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# import table
library(tidyverse)
mesocosm_qpcr_18s = read_csv("qPCR data/Chapter1_Mesocosm_18S.csv")

# clean data
mesocosm_qpcr_18s = as.data.frame(mesocosm_qpcr_18s) 
mesocosm_qpcr_18s$Time = factor (mesocosm_qpcr_18s$Time, levels = c("10m","1h","6h","12h","24h"))
mesocosm_qpcr_18s$Extraction = factor (mesocosm_qpcr_18s$Extraction, levels = c("PW","PS","BT"))
str(mesocosm_qpcr_18s)

# Shapiro-Wilk test (normality) 
shapiro.test(mesocosm_qpcr_18s$Eu_Concentration) 

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(mesocosm_qpcr_18s$Eu_Concentration ~ mesocosm_qpcr_18s$Extraction, data = mesocosm_qpcr_18s)
leveneTest(mesocosm_qpcr_18s$Eu_Concentration ~ mesocosm_qpcr_18s$Time, data = mesocosm_qpcr_18s)

# ANOVA + Tukey HSD
library(multcompView)
summary(aov(Eu_Concentration~Time, data=mesocosm_qpcr_18s)) # p <0.001
TukeyHSD(aov(Eu_Concentration~Time, data=mesocosm_qpcr_18s))
multcompLetters4(aov(Eu_Concentration~Time, data=mesocosm_qpcr_18s), TukeyHSD(aov(Eu_Concentration~Time, data=mesocosm_qpcr_18s)))
# 24h 12h  6h  1h 10m 
#"a" "b" "b" "c" "d" 

summary(aov(Eu_Concentration~Extraction, data=mesocosm_qpcr_18s)) # p <0.001
TukeyHSD(aov(Eu_Concentration~Extraction, data=mesocosm_qpcr_18s))
multcompLetters4(aov(Eu_Concentration~Extraction, data=mesocosm_qpcr_18s), TukeyHSD(aov(Eu_Concentration~Extraction, data=mesocosm_qpcr_18s)))
# BT  PS  PW  
# "a" "b" "b" 

summary(aov(Eu_Concentration~Time * Extraction, data=mesocosm_qpcr_18s)) # time: p < 0.001; extraction: p < 0.001； time and extraction: p < 0.001； 
TukeyHSD(aov(Eu_Concentration~Time * Extraction, data=mesocosm_qpcr_18s))
multcompLetters4(aov(Eu_Concentration~Time * Extraction, data=mesocosm_qpcr_18s), TukeyHSD(aov(Eu_Concentration~Time * Extraction, data=mesocosm_qpcr_18s)))
# 24h 12h  6h  1h 10m 
# "a" "b" "c" "d" "e" 

#  BT  PS  PW 
#"a" "b" "c" 

# 24h:BT 12h:BT 24h:PS  6h:BT 24h:PW 12h:PS 12h:PW  1h:BT  6h:PS  6h:PW  1h:PS  1h:PW 10m:BT 10m:PW 10m:PS 
# "a"    "b"    "c"    "c"   "cd"   "de"   "ef"   "ef"   "ef"    "f"    "g"    "g"    "h"    "h"    "h

# histogram-time and extraction
library(tidyverse)
mesocosm_qpcr_18s_statistics = mesocosm_qpcr_18s %>% 
  group_by(Time, Extraction) %>%
  summarise(
    Mean = mean(Eu_Concentration),
    SD = sd(Eu_Concentration),
    SEM = sd(Eu_Concentration) / sqrt(6),
    .groups = "drop"
  ) 

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
  geom_text(x="10m", y=3, label="h h h", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="1h", y=3, label="g g ef", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="6h", y=3, label="f ef c",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="12h", y=3, label="ef de b",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="24h", y=3, label="cd c a",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x=0.5, y=2, label = "(C)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qpcr_18s_figure.png", plot = mesocosm_qpcr_18s_figure, width = 8, height = 5, units = "cm", dpi = 600)

# histogram-time
mesocosm_qpcr_18s_statistics_time <- mesocosm_qpcr_18s %>%
  group_by(Time) %>%
  summarise(
    Mean = mean(Eu_Concentration),
    SD = sd(Eu_Concentration),
    SEM = sd(Eu_Concentration) / sqrt(18),
    .groups = "drop"
  )

mesocosm_qpcr_18s_figure_time = ggplot(mesocosm_qpcr_18s_statistics_time, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", " 24 hours")) +
  #ggtitle("(A)") + 
  xlab("Submersion time") + 
  ylab("Microeukaryotic eDNA yield\n(million copies/μL)") +
  geom_text(x="10m", y=3, label="e", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="1h", y=3, label="d", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="6h", y=3, label="c",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="12h", y=3, label="b",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="24h", y=3, label="a",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x=0.5, y=3, label = "(C)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_qpcr_18s_figure_time.png", plot = mesocosm_qpcr_18s_figure_time, width = 8, height = 5, units = "cm", dpi = 600)

# histogram-extraction
mesocosm_qpcr_18s_statistics_extraction <- mesocosm_qpcr_18s %>%
  group_by(Extraction) %>%
  summarise(
    Mean = mean(Eu_Concentration),
    SD = sd(Eu_Concentration),
    SEM = sd(Eu_Concentration) / sqrt(30),
    .groups = "drop"
  )

mesocosm_qpcr_18s_figure_extraction =  ggplot(mesocosm_qpcr_18s_statistics_extraction, aes(x = Extraction, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "black", color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 3), breaks = c(0, 1, 2,3)) +
  scale_x_discrete(label = c("PowerWater", "PowerSoil", "Blood&Tissue")) +
  #ggtitle("(A)") + 
  xlab("Extraction method") + 
  ylab("Microeukaryotic eDNA yield\n(million copies/μL)") +
  geom_text(x="PW", y=3, label="c", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="PS", y=3, label="b", size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="BT", y=3, label="a",size =2.5,color = "black",family = "arial", fontface = "plain") +
  geom_text(x=0.5, y=3, label = "(F)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 2.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.text.y =element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) 

ggsave("mesocosm_qpcr_18s_figure_extraction.png", plot = mesocosm_qpcr_18s_figure_extraction, width = 8, height = 5, units = "cm", dpi = 600) 





















### combined figure mesocosm conc--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure_mesocosm_conc_legend = ggplot(mesocosm_qpcr_18s_statistics, aes(x = Time, y = Mean, fill = Extraction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  scale_color_manual(values = c("BT" = "#F27970","PS" = "#BB9727","PW"= "#05B9E2"), labels = c("PowerWater", "PowerSoil", "Blood & Tissue")) +
  scale_fill_manual(values = c("BT" = "#F27970","PS" = "#BB9727","PW"= "#05B9E2"), labels = c("PowerWater", "PowerSoil", "Blood & Tissue")) +
  scale_x_discrete(label = c("1/6", "1", "6", "12", "24")) +
  ggtitle("") + xlab("") + ylab("") +
  theme_bw() + 
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "arial", face = "plain"))

ggsave("figure_mesocosm_conc_legend.png", plot = figure_mesocosm_conc_legend, width = 10, height = 10, units = "cm", dpi = 600)

library(patchwork)
figure_mesocosm_conc = mesocosm_qubit_figure + mesocosm_qpcr_16s_figure + mesocosm_qpcr_18s_figure + plot_layout(ncol = 1, nrow = 3, byrow = T)

ggsave("figure_mesocosm_conc.png", plot = figure_mesocosm_conc, width = 8, height = 15, units = "cm", dpi = 600)

### combined figure mesocosm conc time extraction -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)
figure_mesocosm_conc_time_extraction = mesocosm_qubit_figure_time + mesocosm_qubit_figure_extraction + mesocosm_qpcr_16s_figure_time + mesocosm_qpcr_16s_figure_extraction + mesocosm_qpcr_18s_figure_time + mesocosm_qpcr_18s_figure_extraction + plot_layout(ncol = 2, nrow = 3, byrow = T)

ggsave("figure_mesocosm_conc_time_extraction.png", plot = figure_mesocosm_conc_time_extraction, width = 16, height = 15, units = "cm", dpi = 600)



















### 16s tidy up asv and tax tables ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
packageVersion("tidyverse")
library(vegan)
packageVersion("vegan")

# import tables
asv_16s = read_tsv("Sequence data/asv_16s.tsv")
dim(asv_16s)

tax_16s = read_tsv("Sequence data/tax_16s.tsv")
dim(tax_16s)

# merge tables
asv_tax_16s = inner_join(tax_16s,asv_16s, by = "Feature ID") 
dim(asv_tax_16s) # 31028 ASVs

# remove Chloroplast, Mitochondria and unassigned ASVs
asv_tax_16s_filtered = asv_tax_16s_filtered[!str_detect(asv_tax_16s_filtered$Taxon, "Chloroplast"),] 
dim(asv_tax_16s_filtered) # 380 Chloroplast
asv_tax_16s_filtered = asv_tax_16s_filtered[!str_detect(asv_tax_16s_filtered$Taxon, "Mitochondria"),] 
dim(asv_tax_16s_filtered) # 807 Mitochondria
asv_tax_16s_filtered = asv_tax_16s_filtered[!str_detect(asv_tax_16s_filtered$Taxon, "Unassigned"),] 
dim(asv_tax_16s_filtered) # 70 Unassigned
dim(asv_tax_16s_filtered) # 29771 ASVs left

# decontamination 
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

# remove NC, ID, confidence columns
asv_tax_16s_decontamed = asv_tax_16s_decontamed %>% select(-contains("-NC"))
dim(asv_tax_16s_decontamed)

asv_tax_16s_decontamed = as.data.frame(asv_tax_16s_decontamed) 
rownames(asv_tax_16s_decontamed) = asv_tax_16s_decontamed$`Feature ID`
asv_tax_16s_decontamed = asv_tax_16s_decontamed[,-1]
asv_tax_16s_decontamed = asv_tax_16s_decontamed[,-2]
dim(asv_tax_16s_decontamed)

# produce rarefaction figure
library(vegan)
rarecurve(t(asv_tax_16s_decontamed[,-1]), step=100, lwd=2, ylab="No. of prokaryptic ASVs", xlab="No. of sequences",label =F)

# Separate to mesocosm and field table, remove no longer existed ASVs
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

# rarefy ASV tables, remove no longer existed ASVs
asv_16s_clean_mesocosm_rarefy = rrarefy(asv_16s_clean_mesocosm, min(rowSums(asv_16s_clean_mesocosm)))
dim(asv_16s_clean_mesocosm_rarefy) 
asv_16s_clean_mesocosm_rarefy = asv_16s_clean_mesocosm_rarefy[,which(colSums(asv_16s_clean_mesocosm_rarefy)>0)] 
dim(asv_16s_clean_mesocosm_rarefy) 

asv_16s_clean_field_rarefy = rrarefy(asv_16s_clean_field, min(rowSums(asv_16s_clean_field)))
dim(asv_16s_clean_field_rarefy) 
asv_16s_clean_field_rarefy = asv_16s_clean_field_rarefy[,which(colSums(asv_16s_clean_field_rarefy)>0)] 
dim(asv_16s_clean_field_rarefy) 


















### 18s tidy up asv and tax tables------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(vegan)

# import tables
asv_18s = read_tsv("Sequence data/asv_18s.tsv")
dim(asv_18s)

tax_18s = read_tsv("Sequence data/tax_18s.tsv")
dim(tax_18s)

# merge tables
asv_tax_18s = inner_join(tax_18s,asv_18s, by = "Feature ID") 
dim(asv_tax_18s) # 13146 ASVs

# remove metazoans,land plants,and unassigned ASVs
asv_tax_18s_filtered = asv_tax_18s[!str_detect(asv_tax_18s$Taxon, "Metazoa"),] 
dim(asv_tax_18s_filtered) # 2019 metazoans
asv_tax_18s_filtered = asv_tax_18s_filtered[!str_detect(asv_tax_18s_filtered$Taxon, "Streptophyta"),] 
dim(asv_tax_18s_filtered) # 76 Streptophyta
asv_tax_18s_filtered = asv_tax_18s_filtered[!str_detect(asv_tax_18s_filtered$Taxon, "Unassigned"),] 
dim(asv_tax_18s_filtered) # 51 Unassigned
dim(asv_tax_18s_filtered) # 11000 ASVs left

# decontamination 
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

# remove NC, ID, confidence columns
asv_tax_18s_decontamed = asv_tax_18s_decontamed %>% select(-contains("-NC"))
dim(asv_tax_18s_decontamed)

asv_tax_18s_decontamed = as.data.frame(asv_tax_18s_decontamed) 
rownames(asv_tax_18s_decontamed) = asv_tax_18s_decontamed$`Feature ID`
asv_tax_18s_decontamed = asv_tax_18s_decontamed[,-1]
asv_tax_18s_decontamed = asv_tax_18s_decontamed[,-2]
dim(asv_tax_18s_decontamed)

# produce rarefaction figure
library(vegan)
rarecurve(t(asv_tax_18s_decontamed[,-1]), step=100, lwd=2, ylab="No. of microeukaryptic ASVs", xlab="No. of sequences",label =F)

# Separate to mesocosm and field table, remove no longer existed ASVs
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

# rarefy ASV tables, remove no longer existed ASVs
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
  select(Sample, Type, Water, Site, Repetition)
field_metadata$Type = factor(field_metadata$Type)
field_metadata$Water = factor(field_metadata$Water)
field_metadata$Site = factor(field_metadata$Site)
str(field_metadata)


















### mesocosm 16s richness ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(vegan)
mesocosm_16s_richness = estimateR(asv_16s_clean_mesocosm_rarefy)[1,]
mesocosm_16s_richness = data.frame(Sample = names(mesocosm_16s_richness), Richness = mesocosm_16s_richness)
mesocosm_16s_richness = left_join(mesocosm_metadata, mesocosm_16s_richness, by = "Sample")
str(mesocosm_16s_richness)
mesocosm_16s_richness = as.data.frame(mesocosm_16s_richness) 
str(mesocosm_16s_richness)

# Shapiro-Wilk test (normality) 
shapiro.test(mesocosm_16s_richness$Richness) 

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(subset(mesocosm_16s_richness)$Richness ~ subset(mesocosm_16s_richness)$Time, data = mesocosm_16s_richness)

# ANOVA + Tukey HSD
summary(aov(Richness~Time, data=mesocosm_16s_richness)) # p <0.001
TukeyHSD(aov(Richness~Time, data=mesocosm_16s_richness))

library(multcompView)
multcompLetters4(aov(Richness~Time, data=mesocosm_16s_richness), TukeyHSD(aov(Richness~Time, data=mesocosm_16s_richness)))
# 24h 12h  6h  1h 10m 
#"a"  "a" "ab"  "b"  "c" 

# gam
library(tidyverse)
mesocosm_16s_richness_gam_table = mesocosm_16s_richness %>%
  mutate(Time = case_match(as.character(Time), "10m" ~ 1/6, "1h" ~ 1, "6h" ~ 6,"12h" ~ 12, "24h" ~ 24))
str(mesocosm_16s_richness_gam_table)

library(mgcv)
summary(gam(Richness ~ s(Time,bs="cs", k =3), data=mesocosm_16s_richness_gam_table)) # p < 0.001; R2 (adj) = 0.58

# histogram
library(tidyverse)
mesocosm_16s_richness_statistics = mesocosm_16s_richness %>% group_by(Time) %>%
  summarise(
    Mean = mean(Richness),
    SD = sd(Richness),
    SEM = sd(Richness) / sqrt(6),
    .groups = "drop"
  ) 

mesocosm_16s_richness_figure_histogram = ggplot(mesocosm_16s_richness_statistics, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6,fill = "black",color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),color = "black", position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 2000), breaks = c(0, 500, 1000, 1500, 2000)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", "24 hours")) +
  ggtitle("") + 
  xlab("Submersion time") + 
  ylab("Prokaryotic richness") +
  geom_text(x="10m", y=2000, label="c",size =3,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="1h", y=2000, label="b",size =3,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="6h", y=2000, label="ab",size =3,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="12h", y=2000, label="a",size =3,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="24h", y=2000, label="a",size =3,color = "black",family = "arial", fontface = "plain") +
  #geom_text(x=0.5, y=2.9, label = "(C)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_16s_richness_figure_histogram.png", plot = mesocosm_16s_richness_figure_histogram, width = 8, height = 5, units = "cm", dpi = 600)



















### mesocosm 18s richness ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(vegan)
mesocosm_18s_richness = estimateR(asv_18s_clean_mesocosm_rarefy)[1,]
mesocosm_18s_richness = data.frame(Sample = names(mesocosm_18s_richness), Richness = mesocosm_18s_richness)
mesocosm_18s_richness = left_join(mesocosm_metadata, mesocosm_18s_richness, by = "Sample")
str(mesocosm_18s_richness)
mesocosm_18s_richness = as.data.frame(mesocosm_18s_richness) 
str(mesocosm_18s_richness)

# Shapiro-Wilk test (normality) 
shapiro.test(mesocosm_18s_richness$Richness) 

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(subset(mesocosm_18s_richness)$Richness ~ subset(mesocosm_18s_richness)$Time, data = mesocosm_18s_richness)

# ANOVA + Tukey HSD
summary(aov(Richness~Time, data=mesocosm_18s_richness)) # p < 0.001
TukeyHSD(aov(Richness~Time, data=mesocosm_18s_richness))

library(multcompView)
multcompLetters4(aov(Richness~Time, data=mesocosm_18s_richness), TukeyHSD(aov(Richness~Time, data=mesocosm_18s_richness)))
# 24h 12h  6h  1h 10m 
#"a"  "ab" "ab"  "b"  "c" 

# gam
mesocosm_18s_richness_gam_table = mesocosm_18s_richness %>%
  mutate(Time = case_match(as.character(Time), "10m" ~ 1/6, "1h" ~ 1, "6h" ~ 6,"12h" ~ 12, "24h" ~ 24))
str(mesocosm_18s_richness_gam_table)

library(mgcv)
summary(gam(Richness ~ s(Time,bs="cs", k =3), data=mesocosm_18s_richness_gam_table)) # p < 0.001; R2 (adj) = 0.55

# histogram
library(tidyverse)
mesocosm_18s_richness_statistics = mesocosm_18s_richness %>% group_by(Time) %>%
  summarise(
    Mean = mean(Richness),
    SD = sd(Richness),
    SEM = sd(Richness) / sqrt(6),
    .groups = "drop"
  ) 

mesocosm_18s_richness_figure_histogram = ggplot(mesocosm_18s_richness_statistics, aes(x = Time, y = Mean)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6,fill = "black",color = "black") + 
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),color = "black", position = position_dodge(width = 0.8), width = 0.2, size = 0.3) +
  scale_y_continuous(limits = c(0, 400), breaks = c(0, 100, 200, 300, 400)) +
  scale_x_discrete(label = c("10 mins", "1 hour", "6 hours", "12 hours", "24 hours")) +
  ggtitle("") + 
  xlab("Submersion time") + 
  ylab("Microeukaryotic richness") +
  geom_text(x="10m", y=400, label="c",size =3,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="1h", y=400, label="b",size =3,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="6h", y=400, label="ab",size =3,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="12h", y=400, label="ab",size =3,color = "black",family = "arial", fontface = "plain") +
  geom_text(x="24h", y=400, label="a",size =3,color = "black",family = "arial", fontface = "plain") +
  #geom_text(x=0.5, y=2.9, label = "(C)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  #geom_text(x=0.5, y=28, label = "N = 90", size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  #geom_text(x=0.5, y=27, label = expression(paste("",italic("p"), " < 0.001")), size =2, color = "black",family = "arial", fontface = "plain", hjust = 0) +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  #theme(legend.position = c(0.1, 0.9)) +
  #theme(legend.background = element_rect(fill = NA, color = NA, linewidth = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), linetype = "dotted", color = "gray", size = 0.4) + 
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "bold")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8,family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("mesocosm_18s_richness_figure_histogram.png", plot = mesocosm_18s_richness_figure_histogram, width = 8, height = 5, units = "cm", dpi = 600)


















### combined figure mesocosm richness gam ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# table manipulation 
mesocosm_16s_18s_richness_gam_table = mesocosm_16s_richness_gam_table %>% 
  inner_join(mesocosm_18s_richness_gam_table, by = "Sample") %>% 
  select(Sample, Time.x, Richness.x,Richness.y) 
str(mesocosm_16s_18s_richness_gam_table)

mesocosm_16s_18s_richness_gam_table = mesocosm_16s_18s_richness_gam_table %>%
  rename(Time = Time.x, Pro = Richness.x, Eu = Richness.y)

mesocosm_16s_18s_richness_gam_table = mesocosm_16s_18s_richness_gam_table %>%
  pivot_longer(
    cols = c(Pro, Eu),
    values_to = "Richness",
    names_to = "Group")

mesocosm_16s_18s_richness_gam_table$Group = factor (mesocosm_16s_18s_richness_gam_table$Group, levels = c("Pro","Eu"))
str(mesocosm_16s_18s_richness_gam_table)

mesocosm_16s_18s_richness_gam_table <- mesocosm_16s_18s_richness_gam_table %>%
  pivot_wider(names_from = Group, values_from = Richness) %>%
  mutate(
    Eu_scaled = Eu * 4
  )

# figure
figure_mesocosm_16s_18s_richness_gam = ggplot(mesocosm_16s_18s_richness_gam_table, aes(x = Time)) +
  #geom_point(aes(y = Total), color = "#FFBE7A", size = 2, show.legend = FALSE, alpha = 1) +
  geom_point(aes(y = Pro), color = "#9467BD", size = 2, show.legend = FALSE, alpha = 1) +
  geom_point(aes(y = Eu_scaled), color = "#2CA02C", size = 2, show.legend = FALSE, alpha = 1) +
  #stat_smooth(aes(y = Total), method = "gam", formula = y ~ s(x, bs = "cs", k = 5), color = "#FFBE7A", size = 1.5, se = TRUE, level = 0.95) +
  stat_smooth(aes(y = Pro), method = "gam", formula = y ~ s(x, bs = "cs", k = 3), color = "#9467BD", size = 1.5, se = TRUE, level = 0.95) +
  stat_smooth(aes(y = Eu_scaled), method = "gam", formula = y ~ s(x, bs = "cs", k = 3), color = "#2CA02C", size = 1.5, se = TRUE, level = 0.95) +
  scale_y_continuous(name = "Prokaryotic richness", limits = c(0, 2400), breaks = c(0, 600, 1200, 1800, 2400), sec.axis = sec_axis(~ . /3, breaks = c(0, 200, 400, 600, 800), name = "Microeukaryotic richness")) +
  #scale_y_continuous(name = "Prokaryotic richness", limits = c(0, 2400), breaks = c(0, 600, 1200, 1800, 2400), sec.axis = sec_axis(~ . /2, breaks = c(0, 300, 600, 900, 1200), name = "Eukaryotic richness")) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 6, 12, 18, 24)) +
  #geom_text(x=0.2, y=2300, label = "(B)", size =3.5, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  ggtitle("") + xlab("Submersion time (hour)") + ylab("") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 10,family = "arial", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 10,family = "arial", face = "plain")) +
  theme(axis.text.x =element_text(color="black",size=10,family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=10,family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5))  +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("figure_mesocosm_16s_18s_richness_gam.png", plot = figure_mesocosm_16s_18s_richness_gam, width = 8, height = 8, units = "cm", dpi = 600)

### combined figure mesocosm richness histogram------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)
figure_mesocosm_richness_histogram= mesocosm_16s_richness_figure_histogram + mesocosm_18s_richness_figure_histogram+ plot_layout(ncol = 1, nrow = 2, byrow = T)
ggsave("figure_mesocosm_richness_histogram.png", plot = figure_mesocosm_richness_histogram, width = 8, height = 10, units = "cm", dpi = 600)






### field qubit ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# import table
library(tidyverse)
field_qubit = read_csv("Qubit data/Chapter1_Field_Qubit.csv")

# clean data
field_qubit = as.data.frame(field_qubit) 
field_qubit$Water = factor (field_qubit$Water, levels = c("Clear","Turbid"))
field_qubit$Site = factor (field_qubit$Site, levels = c("1","2","3"))
field_qubit$Sampling = factor (field_qubit$Sampling, levels = c("Active","Passive"))
str(field_qubit)

# calculate mean, SD, and SEM
library(tidyverse)
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

# Shapiro-Wilk test (normality) 
shapiro.test(field_qubit_statistics$Mean) 

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(field_qubit_statistics$Mean ~ field_qubit_statistics$Sampling, data = field_qubit_statistics)

# paired t-test
t.test(subset(field_qubit_statistics, Sampling == "Active")$Mean, subset(field_qubit_statistics, Sampling == "Passive")$Mean, paired = TRUE) # p < 0.01

# boxplot
library(tidyverse)
library(ggsignif)

field_qubit_figure =ggplot(field_qubit_statistics, aes(x=Sampling, y=Mean)) +
  #stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.3, linetype = "solid", size = 0.3, aes(color = factor(Sampling))) +
  geom_point(size = 3, aes(color = factor(Sampling), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color =factor(Sampling)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("Active" = "#1F77B4", "Passive" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  #geom_text(x=0.5, y=45, label = "(A)", size =3, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  geom_signif(comparisons=list(c("Active", "Passive")), annotations="*", y_position = 45,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 50), breaks = c(0, 10, 20, 30, 40, 50)) +
  ggtitle("(A)") + xlab("") + ylab("Total eDNA conc.\n(ng/μL)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain", angle = 45,vjust = 0.9, hjust= 0.8)) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_qubit_figure.png", plot = field_qubit_figure, width = 5, height = 5, units = "cm", dpi = 600) 




















### field qpcr 16S ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# import table
library(tidyverse)
field_qpcr_16s = read_csv("qPCR data/Chapter1_Field_16S.csv")

# clean data
field_qpcr_16s = as.data.frame(field_qpcr_16s) 
field_qpcr_16s$Water = factor (field_qpcr_16s$Water, levels = c("Clear","Turbid"))
field_qpcr_16s$Site = factor (field_qpcr_16s$Site, levels = c("1","2","3"))
field_qpcr_16s$Sampling = factor (field_qpcr_16s$Sampling, levels = c("Active","Passive"))
str(field_qpcr_16s)

# calculate mean, SD, and SEM
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

# Shapiro-Wilk test (normality) 
shapiro.test(field_qpcr_16s_statistics$Mean) 

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(field_qpcr_16s_statistics$Mean ~ field_qpcr_16s_statistics$Sampling, data = field_qpcr_16s_statistics)

# paired t-test
t.test(subset(field_qpcr_16s_statistics, Sampling == "Active")$Mean, subset(field_qpcr_16s_statistics, Sampling == "Passive")$Mean, paired = TRUE) # p < 0.05

# boxplot
library(tidyverse)
library(ggsignif)

field_qpcr_16s_figure = ggplot(field_qpcr_16s_statistics, aes(x=Sampling, y=Mean)) +
  #stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.3, linetype = "solid", size = 0.3, aes(color = factor(Sampling))) +
  geom_point(size = 3, aes(color = factor(Sampling), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color= factor(Sampling)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("Active" = "#1F77B4", "Passive" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("Active", "Passive")), annotations="*", y_position = 83,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 90), breaks = c(0, 30, 60, 90)) +
  ggtitle("(B)") + xlab("") + ylab("Prokaryotic eDNA conc.\n(million copies/μL)") +
  #geom_text(x=0.5, y=85, label = "(B)", size =3, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain", angle = 45,vjust = 0.9, hjust= 0.8)) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_qpcr_16s_figure.png", plot = field_qpcr_16s_figure, width = 5, height = 5, units = "cm", dpi = 600)
















### field qpcr 18S ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# import table
library(tidyverse)
field_qpcr_18s = read_csv("qPCR data/Chapter1_Field_18S.csv")

# clean data
field_qpcr_18s = as.data.frame(field_qpcr_18s) 
field_qpcr_18s$Water = factor (field_qpcr_18s$Water, levels = c("Clear","Turbid"))
field_qpcr_18s$Site = factor (field_qpcr_18s$Site, levels = c("1","2","3"))
field_qpcr_18s$Sampling = factor (field_qpcr_18s$Sampling, levels = c("Active","Passive"))
str(field_qpcr_18s)

# calculate mean, SD, and SEM
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

# Shapiro-Wilk test (normality) 
shapiro.test(field_qpcr_18s_statistics$Mean) 

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(field_qpcr_18s_statistics$Mean ~ field_qpcr_18s_statistics$Sampling, data = field_qpcr_18s_statistics)

# paired t-test
t.test(subset(field_qpcr_18s_statistics, Sampling == "Active")$Mean, subset(field_qpcr_18s_statistics, Sampling == "Passive")$Mean, paired = TRUE) # p < 0.05

# boxplot
library(tidyverse)
library(ggsignif)

field_qpcr_18s_figure = ggplot(field_qpcr_18s_statistics, aes(x=Sampling, y=Mean)) +
  #stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.3, linetype = "solid", size = 0.3, aes(color = factor(Sampling))) +
  geom_point(size = 3, aes(color = factor(Sampling), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color= factor(Sampling)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("Active" = "#1F77B4", "Passive" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("Active", "Passive")), annotations="*", y_position = 5.5,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6)) +
  ggtitle("(C)") + xlab("") + ylab("Microeukaryotic eDNA conc.\n(million copies/μL)") +
  #geom_text(x=0.5, y=5, label = "(C)", size =3, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain", angle = 45,vjust = 0.9, hjust= 0.8)) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_qpcr_18s_figure.png", plot = field_qpcr_18s_figure, width = 5, height = 5, units = "cm", dpi = 600)














### field 16s richness ------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(vegan)
field_16s_richness = estimateR(asv_16s_clean_field_rarefy)[1,]
field_16s_richness = data.frame(Sample = names(field_16s_richness), Richness = field_16s_richness)
field_16s_richness = left_join(field_metadata, field_16s_richness, by = "Sample")
str(field_16s_richness)
field_16s_richness = as.data.frame(field_16s_richness) 
str(field_16s_richness)

# calculate mean, SD, and SEM
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

field_16s_richness_statistics$Water_Site = factor(field_16s_richness_statistics$Water_Site)
str(field_16s_richness_statistics)

# Shapiro-Wilk test (normality) 
shapiro.test(field_16s_richness_statistics$Mean)

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(field_16s_richness_statistics$Mean ~ field_16s_richness_statistics$Type, data = field_16s_richness_statistics)

# Paired t-test
t.test(subset(field_16s_richness_statistics, Type == "F")$Mean, subset(field_16s_richness_statistics, Type == "S")$Mean) # p < 0.05

# plot
library(ggsignif)
field_16s_richness_figure = ggplot(field_16s_richness_statistics, aes(x=Type, y=Mean)) +
  #stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.3, linetype = "solid", size = 0.3, aes(color = factor(Type))) +
  geom_point(size = 3, aes(color = factor(Type), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color = factor(Type)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("F", "S")), annotations="*", y_position = 2300,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 2500), breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
  ggtitle("(D)") + xlab("") + ylab("Prokaryotic richness") +
  #geom_text(x=0.5, y=2000, label = "(D)", size =3, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain", angle = 45,vjust = 0.9, hjust= 0.8)) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_richness_figure.png", plot = field_16s_richness_figure, width = 5, height = 5, units = "cm", dpi = 600)


















### field 18s richness ------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(vegan)
field_18s_richness = estimateR(asv_18s_clean_field)[1,]
field_18s_richness = data.frame(Sample = names(field_18s_richness), Richness = field_18s_richness)
field_18s_richness = left_join(field_metadata, field_18s_richness, by = "Sample")
str(field_18s_richness)
field_18s_richness = as.data.frame(field_18s_richness) 
str(field_18s_richness)

# calculate mean, SD, and SEM
library(tidyverse)
field_18s_richness_statistics = field_18s_richness %>% group_by(Water, Site, Type) %>%
  summarise(
    Mean = mean(Richness),
    SD = sd(Richness),
    SEM = sd(Richness) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 
str(field_18s_richness_statistics)

field_18s_richness_statistics$Water_Site = factor(field_18s_richness_statistics$Water_Site)
str(field_18s_richness_statistics)

# Shapiro-Wilk test (normality) 
shapiro.test(field_18s_richness_statistics$Mean)

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(field_18s_richness_statistics$Mean ~ field_18s_richness_statistics$Type, data = field_18s_richness_statistics)

# Paired t-test
t.test(subset(field_18s_richness_statistics, Type == "F")$Mean, subset(field_18s_richness_statistics, Type == "S")$Mean) # p < 0.05

#plot
library(ggsignif)
field_18s_richness_figure = ggplot(field_18s_richness_statistics, aes(x=Type, y=Mean)) +
  #stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.3, linetype = "solid", size = 0.3, aes(color = factor(Type))) +
  geom_point(size = 3, aes(color = factor(Type), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color = factor(Type)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("F", "S")), annotations="*", y_position = 2300,  tip_length = 0, size = 0.4, vjust=0.4, textsize = 4) +
  scale_y_continuous(limits = c(0, 2500), breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
  ggtitle("(E)") + xlab("") + ylab("Microeukaryotic richness") +
  #geom_text(x=0.5, y=2000, label = "(D)", size =3, color = "black",family = "arial", fontface = "bold", hjust = 0) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain", angle = 45,vjust = 0.9, hjust= 0.8)) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_richness_figure.png", plot = field_18s_richness_figure, width = 5, height = 5, units = "cm", dpi = 600)















### combined figure filed conc and richness ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)
figure_field_conc_richness = field_qubit_figure + field_qpcr_16s_figure + field_qpcr_18s_figure + field_16s_richness_figure + field_18s_richness_figure + plot_layout(ncol =5, nrow = 1,byrow = T)

ggsave("figure_field_conc_richness.png", plot = figure_field_conc_richness, width = 17, height =6, units = "cm", dpi = 600) 




















### field 16s community dissimilarity ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Bray-Curtis distance
library(vegan)
field_16s_bray = vegdist(asv_16s_clean_field_rarefy, method = "bray")
print(field_16s_bray)

# Bray-Curtis distance between AF and PS at each site
# C1
as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray"))))]))
# mean = 0.9163
sd(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray"))))])) / sqrt(9)
# sem = 0.001539892

# C2
as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray"))))]))
# mean = 0.8231
sd(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray"))))])) / sqrt(9)
# sem = 0.002219362

# C3
as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray"))))]))
# mean = 0.9466
sd(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray"))))])) / sqrt(9)
# sem = 0.0006532839

# T1
as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray"))))]))
# mean = 0.9425
sd(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray"))))])) / sqrt(9)
# sem = 0.001962767

# T2
as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray"))))]))
# mean = 0.9737
sd(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray"))))])) / sqrt(9)
# sem = 0.0009288749

# T3
as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray"))))]))
# mean = 0.9504
sd(c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray"))))])) / sqrt(9)
# sem = 0.001175558

# Bray-Curtis distance between AF and PS all sites
field_16s_bray_between_af_ps = c(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C1"),], method = "bray"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C2"),], method = "bray"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"C3"),], method = "bray"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T1"),], method = "bray"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T2"),], method = "bray"))))],
                                 as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"T3"),], method = "bray"))))])

field_16s_bray_between_af_ps = data.frame(value = field_16s_bray_between_af_ps*100, group = "pro")
summary(field_16s_bray_between_af_ps)

# plot 
field_16s_bray_between_af_ps_figure = ggplot(field_16s_bray_between_af_ps, aes(x = group, y = value)) +
  geom_boxplot(width=0.8, fill = "white", outlier.shape = NA, size = 0.5, color = "#9467BD") +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.8, linetype = "dotted", linewidth = 0.3, color = "#9467BD") +
  geom_point(position=position_jitter(width = 0.5, height = 0.5), size = 1, color = "#9467BD", alpha = 0.3, show.legend = F) +
  scale_color_manual(values = "#9467BD") +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("") + xlab("") + ylab("") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.6, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 10,family = "arial", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 5,family = "arial", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=5, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.3)) 

ggsave("field_16s_bray_between_af_ps_figure.png", plot = field_16s_bray_between_af_ps_figure, width = 3, height = 5, units = "cm", dpi = 600)

# PERMANOVA-Permutational Multivariate Analysis of Variance
library(vegan)
adonis2(field_16s_bray ~ field_metadata$Type, permutations = 999) # p < 0.001

# NMDS-Non-metric Multidimensional Scaling
library(car)
library(rgl)

field_16s_bray_nmds = metaMDS(field_16s_bray,k=3) 
field_16s_bray_nmds_table = data.frame(NMDS1 = field_16s_bray_nmds$points[,1], NMDS2 = field_16s_bray_nmds$points[,2],NMDS3 = field_16s_bray_nmds$points[,3])
field_16s_bray_nmds_table$Group = as.factor(field_metadata$Type)
levels(field_16s_bray_nmds_table$Group) <- c(".", "`")
str(field_16s_bray_nmds_table)

# plot
scatter3d(
  x = field_16s_bray_nmds_table$NMDS1, 
  y = field_16s_bray_nmds_table$NMDS2, 
  z = field_16s_bray_nmds_table$NMDS3, 
  groups = field_16s_bray_nmds_table$Group,
  bg.col = "white",
  surface = FALSE, 
  grid = FALSE,
  xlab = "", 
  ylab = "", 
  zlab = "",
  text.col = c("black", "black", "black"),
  axis.scales = FALSE,
  axis.col = c("black", "black", "black"),
  axis.ticks = FALSE,
  surface.col = c("#1F77B4", "#FF7F0E"),
  surface.alpha = 0.2,
  labels = FALSE, 
  ellipsoid = TRUE, 
  level = 0.5, 
  ellipsoid.alpha = 0.2, 
  id = FALSE,
  sphere.size = 1, 
  radius = 10, 
  threshold = 0.01, 
  speed = 10, 
  fov = 60
)


















### field 18s community dissimilarity ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Bray-Curtis distance
library(vegan)
field_18s_bray = vegdist(asv_18s_clean_field_rarefy, method = "bray")
str(field_18s_bray)

# Bray-Curtis distance between AF and PS at each site
# C1
as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray"))))]))
# mean = 0.8240
sd(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray"))))])) / sqrt(9)
# sem = 0.004503852

# C2
as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray"))))]))
# mean = 0.8122
sd(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray"))))])) / sqrt(9)
# sem = 0.00282997

# C3
as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray"))))]))
# mean = 0.8234
sd(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray"))))])) / sqrt(9)
# sem = 0.003105824

# T1
as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray"))))]))
# mean = 0.7155
sd(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray"))))])) / sqrt(9)
# sem = 0.002683083

# T2
as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray"))))]))
# mean = 0.8859
sd(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray"))))])) / sqrt(9)
# sem = 0.002248919

# T3
as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray"))))]

summary(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray"))))]))
# mean = 0.6608
sd(c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray"))))])) / sqrt(9)
# sem = 0.002925192


# Bray-Curtis distance between AF and PS all sites
field_18s_bray_between_af_ps = c(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C1"),], method = "bray"))))],
  as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C2"),], method = "bray"))))],
  as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"C3"),], method = "bray"))))],
  as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T1"),], method = "bray"))))],
  as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T2"),], method = "bray"))))],
  as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray"))[grepl("^F", rownames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray")))), grepl("^S", colnames(as.matrix(vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"T3"),], method = "bray"))))])

field_18s_bray_between_af_ps = data.frame(value = field_18s_bray_between_af_ps*100, group = "eu")
summary(field_18s_bray_between_af_ps)

# plot 
field_18s_bray_between_af_ps_figure = ggplot(field_18s_bray_between_af_ps, aes(x = group, y = value)) +
  geom_boxplot(width=0.8, fill = "white", outlier.shape = NA, size = 0.5, color = "#2CA02C") +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.8, linetype = "dotted", linewidth = 0.3, color = "#2CA02C") +
  geom_point(position=position_jitter(width = 0.5, height = 0.5), size = 1, color = "#2CA02C", alpha = 0.3, show.legend = F) +
  scale_color_manual(values = "#2CA02C") +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("") + xlab("") + ylab("") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.6, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10,family = "arial", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 10,family = "arial", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 5,family = "arial", face = "plain")) +
  theme(axis.text.x =element_blank()) +
  theme(axis.text.y =element_text(color="black",size=5, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.3)) 

ggsave("field_18s_bray_between_af_ps_figure.png", plot = field_18s_bray_between_af_ps_figure, width = 3, height = 5, units = "cm", dpi = 600)

# plot two to one
field_16s_18s_bray_between_af_ps = rbind(field_16s_bray_between_af_ps,field_18s_bray_between_af_ps)
str(field_16s_18s_bray_between_af_ps)
field_16s_18s_bray_between_af_ps$group = factor(field_16s_18s_bray_between_af_ps$group, levels = c("pro", "eu"))
str(field_16s_18s_bray_between_af_ps)

field_16s_18s_bray_between_af_ps_figure = ggplot(field_16s_18s_bray_between_af_ps, aes(x = group, y = value)) +
  geom_boxplot(width=0.8, fill = "white", outlier.shape = NA, size = 0.5) +
  stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.8, linetype = "dotted", size = 0.3) +
  geom_point(position=position_jitter(width = 0.3, height = 0.3), size = 1.2, aes(color = factor(group), fill = factor(group)), alpha = 0.5, shape =21, show.legend = F) +
  scale_color_manual(values = c("pro" = "#9467BD","eu" = "#2CA02C")) +
  scale_fill_manual(values = c("pro" = "#9467BD","eu" = "#2CA02C")) +
  scale_x_discrete(label = c("Prokaryotic\nCommunity","Microeukaryotic\nCommunity")) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ggtitle("") + xlab("") + ylab("Community dissimilarity (%)") +
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4,family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 5,family = "arial", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 8,family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5))+
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5)) 

ggsave("field_16s_18s_bray_between_af_ps_figure.png", plot = field_16s_18s_bray_between_af_ps_figure, width = 4, height = 6, units = "cm", dpi = 600)


# PERMANOVA-Permutational Multivariate Analysis of Variance
library(vegan)
adonis2(field_18s_bray ~ field_metadata$Type, permutations = 999) # p < 0.001

# NMDS-Non-metric Multidimensional Scaling
library(car)
library(rgl)

field_18s_bray_nmds = metaMDS(field_18s_bray,k=3) 
field_18s_bray_nmds_table = data.frame(NMDS1 = field_18s_bray_nmds$points[,1], NMDS2 = field_18s_bray_nmds$points[,2],NMDS3 = field_18s_bray_nmds$points[,3])
field_18s_bray_nmds_table$Group = as.factor(field_metadata$Type)
levels(field_18s_bray_nmds_table$Group) <- c(".", "`")  
str(field_18s_bray_nmds_table)

# plot
scatter3d(
  x = field_18s_bray_nmds_table$NMDS1, 
  y = field_18s_bray_nmds_table$NMDS2, 
  z = field_18s_bray_nmds_table$NMDS3, 
  groups = field_18s_bray_nmds_table$Group,
  bg.col = "white",
  surface = FALSE, 
  grid = FALSE,
  xlab = "", 
  ylab = "", 
  zlab = "",
  text.col = c("black", "black", "black"),
  axis.scales = FALSE,
  axis.col = c("black", "black", "black"),
  axis.ticks = FALSE,
  surface.col = c("#1F77B4", "#FF7F0E"),
  surface.alpha = 0.2,
  labels = FALSE, 
  ellipsoid = TRUE, 
  level = 0.5, 
  ellipsoid.alpha = 0.2, 
  id = FALSE,
  sphere.size = 1, 
  radius = 10, 
  threshold = 0.01, 
  speed = 10, 
  fov = 60
)


















### field 16s relative abundance  --------------------------------------------------------
library(tidyverse)
# calculate relative proportion out of 100
field_16s_rela_abun = asv_16s_clean_field_rarefy / rowSums(asv_16s_clean_field_rarefy) *100
field_16s_rela_abun  = as.data.frame(field_16s_rela_abun)

rowSums(field_16s_rela_abun) # 100

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
rowSums(field_16s_rela_abun_clean) # 100
colnames(field_16s_rela_abun_clean) 

# make top taxa
field_16s_rela_abun_top = field_16s_rela_abun_clean %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:10] %>% data.frame() 
sum(field_16s_rela_abun_top) # 3478.805 out of 3600
rownames(field_16s_rela_abun_top) 

# make below top taxa as others
field_16s_rela_abun_others = rowSums(field_16s_rela_abun_clean[, !(colnames(field_16s_rela_abun_clean) %in% rownames(field_16s_rela_abun_top))])
sum(field_16s_rela_abun_others) # 121.1949
sum(field_16s_rela_abun_top) + sum(field_16s_rela_abun_others) # 3600

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
  theme(legend.title = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "arial", face = "plain",angle = 90, vjust = 0.5)) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_rela_abun_figure_af.png", plot = field_16s_rela_abun_figure_af, width = 8, height = 6, units = "cm", dpi = 600) 

field_16s_rela_abun_figure_ps = ggplot(field_16s_rela_abun_table_ps_longer, aes(x = Sample, y = Proportion, fill = Group)) +
  geom_bar(stat="identity", position = "stack",  width = 0.8, color = "black", linewidth =0.2) +
  scale_fill_manual(values = c("#aec7e8", "#ffbb78", "#C8E6C9", "#ff9896", "#c5b0d5","#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5","#e0e0e080")) +
  ggtitle("")+ xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "arial", face = "plain",angle = 90, vjust = 0.5)) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_rela_abun_figure_ps.png", plot = field_16s_rela_abun_figure_ps, width = 8, height = 6, units = "cm", dpi = 600) 










### field 16s differential taxa --------------------------------------------------------------------------------------------------------------------------------
# average 3 replicates
field_16s_rela_abun_table_mean = field_16s_rela_abun_table %>%
  filter(str_detect(rownames(field_16s_rela_abun_table), "F-C1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-C1") %>%
  bind_rows(field_16s_rela_abun_table, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "F-C2")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-C2") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "F-C3")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-C3") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "F-T1")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-T1") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "F-T2")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-T2") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "F-T3")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-T3") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "S-C1")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-C1") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "S-C2")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-C2") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "S-C3")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-C3") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "S-T1")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-T1") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "S-T2")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-T2") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_16s_rela_abun_table_mean), "S-T3")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-T3") %>%
  bind_rows(field_16s_rela_abun_table_mean, .)

field_16s_rela_abun_table_mean = tail(field_16s_rela_abun_table_mean,12)
rownames(field_16s_rela_abun_table_mean) = field_16s_rela_abun_table_mean$Sample
field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  select(-Sample)

rowSums(field_16s_rela_abun_table_mean) # 100
sum(field_16s_rela_abun_table_mean) # 1200

field_16s_rela_abun_table_mean = field_16s_rela_abun_table_mean %>%
  mutate(Method = c(rep("AF", 6), rep("PS", 6)))

field_16s_rela_abun_table_mean$Method = factor(field_16s_rela_abun_table_mean$Method)
str(field_16s_rela_abun_table_mean)

# Wilcoxon rank sum test 
colnames(field_16s_rela_abun_table_mean)
wilcox.test(field_16s_rela_abun_table_mean$Pseudomonadota[1:6], field_16s_rela_abun_table_mean$Pseudomonadota[7:12], paired = TRUE) # p-value = 0.1563
wilcox.test(field_16s_rela_abun_table_mean$Bacteroidota[1:6], field_16s_rela_abun_table_mean$Bacteroidota[7:12], paired = TRUE) # p-value = 0.1563
wilcox.test(field_16s_rela_abun_table_mean$Cyanobacteriota[1:6], field_16s_rela_abun_table_mean$Cyanobacteriota[7:12], paired = TRUE) # p-value = 0.03125
wilcox.test(field_16s_rela_abun_table_mean$Actinomycetota[1:6], field_16s_rela_abun_table_mean$Actinomycetota[7:12], paired = TRUE) # p-value = 0.03125
wilcox.test(field_16s_rela_abun_table_mean$Thermoplasmatota[1:6], field_16s_rela_abun_table_mean$Thermoplasmatota[7:12], paired = TRUE) # p-value = 0.03125
wilcox.test(field_16s_rela_abun_table_mean$Thermoproteota[1:6], field_16s_rela_abun_table_mean$Thermoproteota[7:12], paired = TRUE) # p-value = 0.4375
wilcox.test(field_16s_rela_abun_table_mean$Planctomycetota[1:6], field_16s_rela_abun_table_mean$Planctomycetota[7:12], paired = TRUE) # p-value = 1
wilcox.test(field_16s_rela_abun_table_mean$Verrucomicrobiota[1:6], field_16s_rela_abun_table_mean$Verrucomicrobiota[7:12], paired = TRUE) # p-value = 0.2188
wilcox.test(field_16s_rela_abun_table_mean$Marinisomatota[1:6], field_16s_rela_abun_table_mean$Marinisomatota[7:12], paired = TRUE) # p-value = 0.03125
wilcox.test(field_16s_rela_abun_table_mean$Myxococcota_A_473307[1:6], field_16s_rela_abun_table_mean$Myxococcota_A_473307[7:12], paired = TRUE) #  p-value = 0.8438
# Cyanobacteriota Actinomycetota Thermoplasmatota Marinisomatota



# manipulate tables
field_16s_rela_abun_table_af_ps_longer = rbind (field_16s_rela_abun_table_af_longer, field_16s_rela_abun_table_ps_longer)
str(field_16s_rela_abun_table_af_ps_longer)

field_16s_rela_abun_table_af_ps_longer = field_16s_rela_abun_table_af_ps_longer %>%
  mutate(Sample_Clean = str_replace(Sample, "-\\d+$", ""))

field_16s_rela_abun_table_af_ps_longer$Method = factor(field_16s_rela_abun_table_af_ps_longer$Method)
str(field_16s_rela_abun_table_af_ps_longer)

field_16s_rela_abun_table_af_ps_longer$Sample_Clean = factor(field_16s_rela_abun_table_af_ps_longer$Sample_Clean)
str(field_16s_rela_abun_table_af_ps_longer)

# calculate mean, SD, and SEM
field_16s_rela_abun_table_af_ps_statistics = field_16s_rela_abun_table_af_ps_longer %>% group_by(Group, Sample_Clean, Method) %>%
  summarise(
    Mean = mean(Proportion),
    SD = sd(Proportion),
    SEM = sd(Proportion) / sqrt(3), 
    .groups = "drop") 

str(field_16s_rela_abun_table_af_ps_statistics)

field_16s_rela_abun_table_af_ps_statistics = field_16s_rela_abun_table_af_ps_statistics %>%
  mutate(Water = case_when(
    str_starts(Sample_Clean, "OS") ~ "OS",
    str_starts(Sample_Clean, "ES") ~ "ES"))

field_16s_rela_abun_table_af_ps_statistics$Water = factor(field_16s_rela_abun_table_af_ps_statistics$Water)
str(field_16s_rela_abun_table_af_ps_statistics)


# plot
library(tidyverse)
library(ggsignif)

# Cyanobacteriota
field_16s_differential_taxa_cyanobacteriota_figure = ggplot(field_16s_rela_abun_table_af_ps_statistics[field_16s_rela_abun_table_af_ps_statistics$Group == "Cyanobacteriota", ], aes(x=Method, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Method), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(color = factor(Method), ymin = Mean-SEM, ymax = Mean+SEM), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Sample_Clean), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_manual(values = c("AF" = "#1F77B4", "PS" = "#FF7F0E")) +
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("AF", "PS")), annotations="*", y_position = 45,  tip_length = 0, size = 0.3, vjust=0, textsize = 4) +
  scale_y_continuous(limits = c(0, 50), breaks = c(0,10, 20, 30, 40, 50)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "#C8E6C9", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_differential_taxa_cyanobacteriota_figure.png", plot = field_16s_differential_taxa_cyanobacteriota_figure, width = 4, height = 6, units = "cm", dpi = 600)

# Actinomycetota
field_16s_differential_taxa_actinomycetota_figure = ggplot(field_16s_rela_abun_table_af_ps_statistics[field_16s_rela_abun_table_af_ps_statistics$Group == "Actinomycetota", ], aes(x=Method, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Method), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(color = factor(Method), ymin = Mean-SEM, ymax = Mean+SEM), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Sample_Clean), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_manual(values = c("AF" = "#1F77B4", "PS" = "#FF7F0E")) +
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("AF", "PS")), annotations="*", y_position = 10,  tip_length = 0, size = 0.3, vjust=0, textsize = 4) +
  scale_y_continuous(limits = c(0, 12), breaks = c(0,4, 8, 12)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "#ff9896", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_differential_taxa_actinomycetota_figure.png", plot = field_16s_differential_taxa_actinomycetota_figure, width = 4, height = 6, units = "cm", dpi = 600)


# Thermoplasmatota
field_16s_differential_taxa_thermoplasmatota_figure = ggplot(field_16s_rela_abun_table_af_ps_statistics[field_16s_rela_abun_table_af_ps_statistics$Group == "Thermoplasmatota", ], aes(x=Method, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Method), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(color = factor(Method), ymin = Mean-SEM, ymax = Mean+SEM), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Sample_Clean), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_manual(values = c("AF" = "#1F77B4", "PS" = "#FF7F0E")) +
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("AF", "PS")), annotations="*", y_position = 16,  tip_length = 0, size = 0.3, vjust=0, textsize = 4) +
  scale_y_continuous(limits = c(0, 18), breaks = c(0, 9, 18)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "#c5b0d5", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_differential_taxa_thermoplasmatota_figure.png", plot = field_16s_differential_taxa_thermoplasmatota_figure, width = 4, height = 6, units = "cm", dpi = 600)

# Marinisomatota
field_16s_differential_taxa_marinisomatota_figure = ggplot(field_16s_rela_abun_table_af_ps_statistics[field_16s_rela_abun_table_af_ps_statistics$Group == "Marinisomatota", ], aes(x=Method, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Method), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(color = factor(Method), ymin = Mean-SEM, ymax = Mean+SEM), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Sample_Clean), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_manual(values = c("AF" = "#1F77B4", "PS" = "#FF7F0E")) +
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("AF", "PS")), annotations="*", y_position = 3.5,  tip_length = 0, size = 0.3, vjust=0, textsize = 4) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1, 2,3, 4)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "#dbdb8d", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_differential_taxa_marinisomatota_figure.png", plot = field_16s_differential_taxa_marinisomatota_figure, width = 4, height = 6, units = "cm", dpi = 600)














### combined figure field 16s relative abundance and differential taxa  --------------------------------------------------------
library(patchwork)

figure_field_16s_rela_abun = (field_16s_rela_abun_figure_af + field_16s_rela_abun_figure_ps) + plot_layout(ncol = 1, nrow = 2,byrow = T)
figure_field_16s_diff_taxa = (field_16s_differential_taxa_cyanobacteriota_figure + field_16s_differential_taxa_actinomycetota_figure + field_16s_differential_taxa_thermoplasmatota_figure + field_16s_differential_taxa_marinisomatota_figure) + plot_layout(ncol = 2, nrow = 2,byrow = T)
figure_field_16s_rela_abun_diff_taxa = (figure_field_16s_rela_abun | figure_field_16s_diff_taxa) + plot_layout(widths = c(1, 1))

ggsave("figure_field_16s_rela_abun_diff_taxa.png", plot = figure_field_16s_rela_abun_diff_taxa, width = 18, height = 12, units = "cm", dpi = 600)


















### field 18s relative abundance  --------------------------------------------------------
library(tidyverse)
# calculate relative proportion out of 100
field_18s_rela_abun = asv_18s_clean_field_rarefy / rowSums(asv_18s_clean_field_rarefy) *100
field_18s_rela_abun  = as.data.frame(field_18s_rela_abun)

rowSums(field_18s_rela_abun) # 100

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
rowSums(field_18s_rela_abun_clean) # 100
colnames(field_18s_rela_abun_clean) 

# make top taxa
field_18s_rela_abun_top = field_18s_rela_abun_clean %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:11] %>% data.frame() 
sum(field_18s_rela_abun_top) # 3467.537 out of 3600
rownames(field_18s_rela_abun_top) 

# make below top taxa as others
field_18s_rela_abun_others = rowSums(field_18s_rela_abun_clean[, !(colnames(field_18s_rela_abun_clean) %in% rownames(field_18s_rela_abun_top))])
sum(field_18s_rela_abun_others) # 132.4627
sum(field_18s_rela_abun_top) + sum(field_18s_rela_abun_others) # 3600

# add existed "Others" in top to Others
field_18s_rela_abun_others = field_18s_rela_abun_others + field_18s_rela_abun_clean[, "Others"]
sum(field_18s_rela_abun_others) # 384.946

# produce top+Others table
field_18s_rela_abun_table = cbind(field_18s_rela_abun_clean[,rownames(field_18s_rela_abun_top)] %>% select(-"Others"), Others = field_18s_rela_abun_others)
rowSums(field_18s_rela_abun_table) # 100
sum(field_18s_rela_abun_table) # 3600

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
  theme(legend.title = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "arial", face = "plain",angle = 90, vjust = 0.5)) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_rela_abun_figure_af.png", plot = field_18s_rela_abun_figure_af, width = 8, height = 6, units = "cm", dpi = 600) 

field_18s_rela_abun_figure_ps = ggplot(field_18s_rela_abun_table_ps_longer, aes(x = Sample, y = Proportion, fill = Group)) +
  geom_bar(stat="identity", position = "stack",  width = 0.8, color = "black", linewidth =0.2) +
  scale_fill_manual(values = c("#aec7e8", "#ffbb78", "#C8E6C9", "#ff9896", "#c5b0d5","#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5","#e0e0e080")) +
  ggtitle("")+ xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 10, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "arial", face = "plain",angle = 90, vjust = 0.5)) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_rela_abun_figure_ps.png", plot = field_18s_rela_abun_figure_ps, width = 8, height = 6, units = "cm", dpi = 600) 
















### field 18s differential taxa --------------------------------------------------------------------------------------------------------------------------------
# average 3 replicates
field_18s_rela_abun_table_mean = field_18s_rela_abun_table %>%
  filter(str_detect(rownames(field_18s_rela_abun_table), "F-C1")) %>%
  summarise(across(everything(), mean)) %>%
  mutate(Sample = "F-C1") %>%
  bind_rows(field_18s_rela_abun_table, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "F-C2")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-C2") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "F-C3")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-C3") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "F-T1")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-T1") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "F-T2")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-T2") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "F-T3")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "F-T3") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "S-C1")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-C1") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "S-C2")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-C2") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "S-C3")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-C3") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "S-T1")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-T1") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "S-T2")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-T2") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  filter(str_detect(rownames(field_18s_rela_abun_table_mean), "S-T3")) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(Sample = "S-T3") %>%
  bind_rows(field_18s_rela_abun_table_mean, .)

field_18s_rela_abun_table_mean = tail(field_18s_rela_abun_table_mean,12)
rownames(field_18s_rela_abun_table_mean) = field_18s_rela_abun_table_mean$Sample
field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  select(-Sample)

rowSums(field_18s_rela_abun_table_mean) # 100
sum(field_18s_rela_abun_table_mean) # 1200

field_18s_rela_abun_table_mean = field_18s_rela_abun_table_mean %>%
  mutate(Method = c(rep("AF", 6), rep("PS", 6)))

field_18s_rela_abun_table_mean$Method = factor(field_18s_rela_abun_table_mean$Method)
str(field_18s_rela_abun_table_mean)

# Wilcoxon rank sum test
colnames(field_18s_rela_abun_table_mean)
wilcox.test(field_18s_rela_abun_table$Gyrista[1:6], field_18s_rela_abun_table$Gyrista[7:12], paired = TRUE) # p-value = 0.6875
wilcox.test(field_18s_rela_abun_table$Ciliophora[1:6], field_18s_rela_abun_table$Ciliophora[7:12], paired = TRUE) # p-value = 0.03125
wilcox.test(field_18s_rela_abun_table$Dinoflagellata[1:6], field_18s_rela_abun_table$Dinoflagellata[7:12], paired = TRUE) # p-value = 0.4375
wilcox.test(field_18s_rela_abun_table$Bigyra[1:6], field_18s_rela_abun_table$Bigyra[7:12], paired = TRUE) # p-value =  0.0625
wilcox.test(field_18s_rela_abun_table$Chlorophyta_X[1:6], field_18s_rela_abun_table$Chlorophyta_X[7:12], paired = TRUE) # p-value = 1
wilcox.test(field_18s_rela_abun_table$Cryptophyta_X[1:6], field_18s_rela_abun_table$Cryptophyta_X[7:12], paired = TRUE) # p-value = 0.03125
wilcox.test(field_18s_rela_abun_table$Cercozoa[1:6], field_18s_rela_abun_table$Cercozoa[7:12], paired = TRUE) # p-value = 0.03125
wilcox.test(field_18s_rela_abun_table$Haptophyta_X[1:6], field_18s_rela_abun_table$Haptophyta_X[7:12], paired = TRUE) # p-value = 0.4375
wilcox.test(field_18s_rela_abun_table$Rhodophyta_X[1:6], field_18s_rela_abun_table$Rhodophyta_X[7:12], paired = TRUE) # p-value = 0.4004
wilcox.test(field_18s_rela_abun_table$Picozoa_X[1:6], field_18s_rela_abun_table$Picozoa_X[7:12], paired = TRUE) #  p-value = 0.03125
# Ciliophora Cryptophyta Cercozoa Picozoa

# manipulate tables
field_18s_rela_abun_table_af_ps_longer = rbind (field_18s_rela_abun_table_af_longer, field_18s_rela_abun_table_ps_longer)
str(field_18s_rela_abun_table_af_ps_longer)

field_18s_rela_abun_table_af_ps_longer = field_18s_rela_abun_table_af_ps_longer %>%
  mutate(Sample_Clean = str_replace(Sample, "-\\d+$", ""))

field_18s_rela_abun_table_af_ps_longer$Method = factor(field_18s_rela_abun_table_af_ps_longer$Method)
str(field_18s_rela_abun_table_af_ps_longer)

field_18s_rela_abun_table_af_ps_longer$Sample_Clean = factor(field_18s_rela_abun_table_af_ps_longer$Sample_Clean)
str(field_18s_rela_abun_table_af_ps_longer)

# calculate mean, SD, and SEM
field_18s_rela_abun_table_af_ps_statistics = field_18s_rela_abun_table_af_ps_longer %>% group_by(Group, Sample_Clean, Method) %>%
  summarise(
    Mean = mean(Proportion),
    SD = sd(Proportion),
    SEM = sd(Proportion) / sqrt(3), 
    .groups = "drop") 

str(field_18s_rela_abun_table_af_ps_statistics)

field_18s_rela_abun_table_af_ps_statistics = field_18s_rela_abun_table_af_ps_statistics %>%
  mutate(Water = case_when(
    str_starts(Sample_Clean, "OS") ~ "OS",
    str_starts(Sample_Clean, "ES") ~ "ES"))

field_18s_rela_abun_table_af_ps_statistics$Water = factor(field_18s_rela_abun_table_af_ps_statistics$Water)
str(field_18s_rela_abun_table_af_ps_statistics)


# plot
library(tidyverse)
library(ggsignif)

# Ciliophora
field_18s_differential_taxa_ciliophora_figure = ggplot(field_18s_rela_abun_table_af_ps_statistics[field_18s_rela_abun_table_af_ps_statistics$Group == "Ciliophora", ], aes(x=Method, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Method), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(color = factor(Method), ymin = Mean-SEM, ymax = Mean+SEM), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Sample_Clean), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_manual(values = c("AF" = "#1F77B4", "PS" = "#FF7F0E")) +
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("AF", "PS")), annotations="*", y_position = 67,  tip_length = 0, size = 0.3, vjust=0, textsize = 4) +
  scale_y_continuous(limits = c(0, 75), breaks = c(0,25, 50, 75)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_differential_taxa_ciliophora_figure.png", plot = field_18s_differential_taxa_ciliophora_figure, width = 4, height = 6, units = "cm", dpi = 600)

# Cryptophyta
field_18s_differential_taxa_cryptophyta_figure = ggplot(field_18s_rela_abun_table_af_ps_statistics[field_18s_rela_abun_table_af_ps_statistics$Group == "Cryptophyta_X", ], aes(x=Method, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Method), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(color = factor(Method), ymin = Mean-SEM, ymax = Mean+SEM), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Sample_Clean), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_manual(values = c("AF" = "#1F77B4", "PS" = "#FF7F0E")) +
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("AF", "PS")), annotations="*", y_position = 7.5,  tip_length = 0, size = 0.3, vjust=0, textsize = 4) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0,3, 6,9)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_differential_taxa_cryptophyta_figure.png", plot = field_18s_differential_taxa_cryptophyta_figure, width = 4, height = 6, units = "cm", dpi = 600)


# Cercozoa Picozoa
field_18s_differential_taxa_cercozoa_figure = ggplot(field_18s_rela_abun_table_af_ps_statistics[field_18s_rela_abun_table_af_ps_statistics$Group == "Cercozoa", ], aes(x=Method, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Method), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(color = factor(Method), ymin = Mean-SEM, ymax = Mean+SEM), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Sample_Clean), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_manual(values = c("AF" = "#1F77B4", "PS" = "#FF7F0E")) +
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("AF", "PS")), annotations="*", y_position = 9,  tip_length = 0, size = 0.3, vjust=0, textsize = 4) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0,5, 10)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_differential_taxa_cercozoa_figure.png", plot = field_18s_differential_taxa_cercozoa_figure, width = 4, height = 6, units = "cm", dpi = 600)

# Picozoa
field_18s_differential_taxa_picozoa_figure = ggplot(field_18s_rela_abun_table_af_ps_statistics[field_18s_rela_abun_table_af_ps_statistics$Group == "Picozoa_X", ], aes(x=Method, y=Mean)) +
  geom_point(size = 2, aes(color = factor(Method), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(color = factor(Method), ymin = Mean-SEM, ymax = Mean+SEM), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Sample_Clean), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_shape_manual(values = c(15, 16)) + 
  scale_color_manual(values = c("AF" = "#1F77B4", "PS" = "#FF7F0E")) +
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("AF", "PS")), annotations="*", y_position = 8,  tip_length = 0, size = 0.3, vjust=0, textsize = 4) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 3, 6, 9)) +
  ggtitle("") + xlab("") + ylab("Relative abundance (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_differential_taxa_picozoa_figure.png", plot = field_18s_differential_taxa_picozoa_figure, width = 4, height = 6, units = "cm", dpi = 600)












### combined figure field 18s relative abundance and differential taxa  --------------------------------------------------------
library(patchwork)

figure_field_18s_rela_abun = (field_18s_rela_abun_figure_af + field_18s_rela_abun_figure_ps) + plot_layout(ncol = 1, nrow = 2,byrow = T)
figure_field_18s_diff_taxa = (field_18s_differential_taxa_ciliophora_figure + field_18s_differential_taxa_cryptophyta_figure + field_18s_differential_taxa_cercozoa_figure + field_18s_differential_taxa_picozoa_figure) + plot_layout(ncol = 2, nrow = 2,byrow = T)
figure_field_18s_rela_abun_diff_taxa = (figure_field_18s_rela_abun | figure_field_18s_diff_taxa) + plot_layout(widths = c(1, 1))

ggsave("figure_field_18s_rela_abun_diff_taxa.png", plot = figure_field_18s_rela_abun_diff_taxa, width = 18, height = 12, units = "cm", dpi = 600)















### field 16s community similarity among 3 biological replicates ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# active filtration
field_16s_dissimi_f_among_3replicates = c(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-C1"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-C1"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-C2"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-C2"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-C3"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-C3"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-T1"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-T1"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-T2"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-T2"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-T3"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"F-T3"),], method = "bray")))])
print(field_16s_dissimi_f_among_3replicates)
summary (field_16s_dissimi_f_among_3replicates)

field_16s_simi_f_among_3replicates = (1-field_16s_dissimi_f_among_3replicates)*100
print(field_16s_simi_f_among_3replicates)
summary (field_16s_simi_f_among_3replicates)

field_16s_simi_f_among_3replicates_table = field_metadata[1:18, ]
field_16s_simi_f_among_3replicates_table$Similarity <- field_16s_simi_f_among_3replicates
str(field_16s_simi_f_among_3replicates_table)

# passive sampling
field_16s_dissimi_s_among_3replicates = c(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-C1"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-C1"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-C2"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-C2"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-C3"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-C3"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-T1"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-T1"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-T2"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-T2"),], method = "bray")))],
                                          as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-T3"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_16s_clean_field_rarefy[str_detect(rownames(asv_16s_clean_field_rarefy),"S-T3"),], method = "bray")))])
print(field_16s_dissimi_s_among_3replicates)
summary (field_16s_dissimi_s_among_3replicates)

field_16s_simi_s_among_3replicates = (1-field_16s_dissimi_s_among_3replicates)*100
print(field_16s_simi_s_among_3replicates)
summary (field_16s_simi_s_among_3replicates)

field_16s_simi_s_among_3replicates_table <- field_metadata[19:36, ]
field_16s_simi_s_among_3replicates_table$Similarity <- field_16s_simi_s_among_3replicates
str(field_16s_simi_s_among_3replicates_table)

#combine two tables
field_16s_simi_f_s_among_3replicates_table = rbind(field_16s_simi_f_among_3replicates_table, field_16s_simi_s_among_3replicates_table)
str(field_16s_simi_f_s_among_3replicates_table)

# mean SD SEM
field_16s_simi_f_s_among_3replicates_table_statistics = field_16s_simi_f_s_among_3replicates_table %>% group_by(Water, Site, Type) %>%
  summarise(
    Mean = mean(Similarity),
    SD = sd(Similarity),
    SEM = sd(Similarity) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 
str(field_16s_simi_f_s_among_3replicates_table_statistics)

field_16s_simi_f_s_among_3replicates_table_statistics$Water_Site = factor(field_16s_simi_f_s_among_3replicates_table_statistics$Water_Site)
str(field_16s_simi_f_s_among_3replicates_table_statistics)


# Shapiro-Wilk test (normality) 
shapiro.test(field_16s_simi_f_s_among_3replicates_table_statistics$Mean)

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(field_16s_simi_f_s_among_3replicates_table_statistics$Mean ~ field_16s_simi_f_s_among_3replicates_table_statistics$Type, data = field_16s_simi_f_s_among_3replicates_table_statistics) 

# Paired t-test
t.test(subset(field_16s_simi_f_s_among_3replicates_table_statistics, Type == "F")$Mean, subset(field_16s_simi_f_s_among_3replicates_table_statistics, Type == "S")$Mean, paired = TRUE) # p > 0.05

# plot
library(tidyverse)
library(ggsignif)

field_16s_simi_among_3replicates_figure =ggplot(field_16s_simi_f_s_among_3replicates_table_statistics, aes(x=Type, y=Mean)) +
  #stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.3, linetype = "solid", size = 0.3, aes(color = factor(Type))) +
  geom_point(size = 2, aes(color = factor(Type), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color = factor(Type)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("F", "S")), annotations="ns", y_position = 100,  tip_length = 0, size = 0.3, vjust=-0.2, textsize = 3) +
  scale_y_continuous(limits = c(75, 105), breaks = c(80, 85, 90, 95, 100)) +
  ggtitle("(A)") + xlab("") + ylab("Prokaryotic community similarity (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_16s_simi_among_3replicates_figure.png", plot = field_16s_simi_among_3replicates_figure, width = 5, height = 6, units = "cm", dpi = 600)



















### field 18s community similarity among 3 biological replicates ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# active filtration
field_18s_dissimi_f_among_3replicates = c(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-C1"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-C1"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-C2"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-C2"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-C3"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-C3"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-T1"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-T1"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-T2"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-T2"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-T3"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"F-T3"),], method = "bray")))])
print(field_18s_dissimi_f_among_3replicates)
summary (field_18s_dissimi_f_among_3replicates)

field_18s_simi_f_among_3replicates = (1-field_18s_dissimi_f_among_3replicates)*100
print(field_18s_simi_f_among_3replicates)
summary (field_18s_simi_f_among_3replicates)

field_18s_simi_f_among_3replicates_table = field_metadata[1:18, ]
field_18s_simi_f_among_3replicates_table$Similarity <- field_18s_simi_f_among_3replicates
str(field_18s_simi_f_among_3replicates_table)

# passive sampling
field_18s_dissimi_s_among_3replicates = c(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-C1"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-C1"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-C2"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-C2"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-C3"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-C3"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-T1"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-T1"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-T2"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-T2"),], method = "bray")))],
                                          as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-T3"),], method = "bray"))[lower.tri(as.matrix (vegdist(asv_18s_clean_field_rarefy[str_detect(rownames(asv_18s_clean_field_rarefy),"S-T3"),], method = "bray")))])
print(field_18s_dissimi_s_among_3replicates)
summary (field_18s_dissimi_s_among_3replicates)

field_18s_simi_s_among_3replicates = (1-field_18s_dissimi_s_among_3replicates)*100
print(field_18s_simi_s_among_3replicates)
summary (field_18s_simi_s_among_3replicates)

field_18s_simi_s_among_3replicates_table <- field_metadata[19:36, ]
field_18s_simi_s_among_3replicates_table$Similarity <- field_18s_simi_s_among_3replicates
str(field_18s_simi_s_among_3replicates_table)

#combine two tables
field_18s_simi_f_s_among_3replicates_table = rbind(field_18s_simi_f_among_3replicates_table, field_18s_simi_s_among_3replicates_table)
str(field_18s_simi_f_s_among_3replicates_table)

# mean SD SEM
field_18s_simi_f_s_among_3replicates_table_statistics = field_18s_simi_f_s_among_3replicates_table %>% group_by(Water, Site, Type) %>%
  summarise(
    Mean = mean(Similarity),
    SD = sd(Similarity),
    SEM = sd(Similarity) / sqrt(3), 
    .groups = "drop"
  ) %>%
  mutate(Water_Site = paste(Water, Site, sep = "_")) 
str(field_18s_simi_f_s_among_3replicates_table_statistics)

field_18s_simi_f_s_among_3replicates_table_statistics$Water_Site = factor(field_18s_simi_f_s_among_3replicates_table_statistics$Water_Site)
str(field_18s_simi_f_s_among_3replicates_table_statistics)


# Shapiro-Wilk test (normality) 
shapiro.test(field_18s_simi_f_s_among_3replicates_table_statistics$Mean)

# Levene test (Homogeneity of Variance)
library(car)
leveneTest(field_18s_simi_f_s_among_3replicates_table_statistics$Mean ~ field_18s_simi_f_s_among_3replicates_table_statistics$Type, data = field_18s_simi_f_s_among_3replicates_table_statistics) 

# Paired t-test
t.test(subset(field_18s_simi_f_s_among_3replicates_table_statistics, Type == "F")$Mean, subset(field_18s_simi_f_s_among_3replicates_table_statistics, Type == "S")$Mean, paired = TRUE) # p > 0.05

# plot
library(tidyverse)
library(ggsignif)

field_18s_simi_among_3replicates_figure =ggplot(field_18s_simi_f_s_among_3replicates_table_statistics, aes(x=Type, y=Mean)) +
  #stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = 0.3, linetype = "solid", size = 0.3, aes(color = factor(Type))) +
  geom_point(size = 2, aes(color = factor(Type), shape = factor(Water)), show.legend = F, alpha = 1) +
  geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM, color = factor(Type)), width = 0.03, size = 0.2, alpha = 1) +
  geom_line(aes(group = Water_Site), color = "grey", linetype = "solid", size = 0.3, alpha = 1) +
  scale_color_manual(values = c("F" = "#1F77B4", "S" = "#FF7F0E")) +
  scale_shape_manual(values = c(16, 15)) + 
  scale_x_discrete(label = c("Active\nFiltration","Passive\nSampling")) +
  geom_signif(comparisons=list(c("F", "S")), annotations="ns", y_position = 100,  tip_length = 0, size = 0.3, vjust=-0.2, textsize = 3) +
  scale_y_continuous(limits = c(75, 105), breaks = c(80, 85, 90, 95, 100)) +
  ggtitle("(B)") + xlab("") + ylab("Microeukaryotic community similarity (%)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid")) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 4, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 8, family = "arial", face = "bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(color = "black", size = 8, family = "arial", face = "plain")) +
  theme(axis.text.x = element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.5)) +
  theme(axis.ticks.y = element_line(color = "black", linewidth = 0.5))

ggsave("field_18s_simi_among_3replicates_figure.png", plot = field_18s_simi_among_3replicates_figure, width = 5, height = 6, units = "cm", dpi = 600)















### combined figure community similarity among 3 biological replicates --------------------------------------------------------------------------------------------------------------
library(patchwork)

figure_field_dissimilarity <- field_16s_simi_among_3replicates_figure + field_18s_simi_among_3replicates_figure + plot_layout(ncol = 2, nrow = 1, byrow = T)
ggsave("figure_field_dissimilarity.png", plot = figure_field_dissimilarity, width = 12, height =8, units = "cm", dpi = 600)













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

# plot
mantal_combined_plot = mantal_combined %>% 
  mutate(rd = cut(Mantel_r, breaks = c(-Inf, 0.6, Inf),
                  labels = c("< 0.6", "≥ 0.6")),
         pd = cut(p_value, breaks = c(-Inf, 0.05, Inf),
                  labels = c("< 0.05","≥ 0.05")))

field_mantel_af_figure = qcorrplot(correlate(env_factor, method = "pearson"), type = "lower", diag = TRUE) +
  geom_square() + 
  geom_couple(aes(colour = pd, size = rd), data = mantal_combined_plot[1:12,], 
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


ggsave("field_mantel_af_figure.png", plot = field_mantel_af_figure, width = 8, height = 10, units = "cm", dpi = 600)

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

ggsave("field_mantel_ps_figure.png", plot = field_mantel_ps_figure, width = 8, height = 8, units = "cm", dpi = 600)














### field 16s lefse ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(phyloseq)
# prepare asv table 
class (asv_16s_clean_field_rarefy)
class (asv_16s_clean_field_rarefy)
# prepare tax table
library(tidyverse)
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


# lefser
#library(lefser)
#library(mia)
# active filtration
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

# plot

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
  theme(legend.title = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 0, family = "arial", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=0, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.3)) +
  theme(axis.ticks.y = element_blank()) 

ggsave("lefse_16s_af_figure.png", plot = lefse_16s_af_figure, width = 5, height = 5, units = "cm", dpi = 600)

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
  theme(legend.title = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 0, family = "arial", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=0, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.3)) +
  theme(axis.ticks.y = element_blank()) 

ggsave("lefse_16s_ps_figure.png", plot = lefse_16s_ps_figure, width = 5, height = 5, units = "cm", dpi = 600)






















### field 18s lefse ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(phyloseq)
# prepare asv table 
class (asv_18s_clean_field_rarefy)
class (asv_18s_clean_field_rarefy)
# prepare tax table
library(tidyverse)
tax_18s_seperate_taxonomy = tax_18s %>%
  separate(Taxon, 
           into = c("Domain", "Supergroup", "Division", "Subdivision", 
                    "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";", 
           remove = FALSE,  # Keep the original Taxon column
           fill = "right")  # Fill missing ranks with NA on the right

tax_18s_seperate_taxonomy = tax_18s_seperate_taxonomy %>%
  mutate(across(c(Domain, Supergroup, Division, Subdivision, 
                  Class, Order, Family, Genus, Species), 
                ~gsub("^[a-z]+__", "", .)))

class(tax_18s_seperate_taxonomy)
tax_18s_seperate_taxonomy = as.data.frame(tax_18s_seperate_taxonomy)
class(tax_18s_seperate_taxonomy)
rownames(tax_18s_seperate_taxonomy) = tax_18s_seperate_taxonomy[,"Feature ID"] 
tax_18s_clean_field_rarefy = tax_18s_seperate_taxonomy %>%
  filter(`Feature ID` %in% rownames(t(asv_18s_clean_field_rarefy))) # remove extra taxa

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
# active filtration
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

# passive sampling
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
  theme(legend.title = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 0, family = "arial", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=0, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.3)) +
  theme(axis.ticks.y = element_blank()) 

ggsave("lefse_18s_af_figure.png", plot = lefse_18s_af_figure, width = 4, height = 6, units = "cm", dpi = 600)

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
  theme(legend.title = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(legend.text = element_text(colour = "black", size = 10, family = "arial", face = "plain")) +
  theme(plot.title = element_text(color = "black", size = 0, family = "arial", face = "plain")) +
  theme(axis.title.x = element_text(color = "black", size = 8)) +
  theme(axis.title.y = element_text(color = "black", size = 8)) +
  theme(axis.text.x =element_text(color="black",size=8, family = "arial", face = "plain")) +
  theme(axis.text.y =element_text(color="black",size=0, family = "arial", face = "plain")) +
  theme(axis.ticks.x = element_line(color = "black", linewidth = 0.3)) +
  theme(axis.ticks.y = element_blank())          

ggsave("lefse_18s_ps_figure.png", plot = lefse_18s_ps_figure, width = 4, height = 6, units = "cm", dpi = 600)













### combined figure lefse 16s 18s --------------------------------------------------------------------------------------------------------------
library(patchwork)

figure_field_lefse = lefse_16s_af_figure + lefse_16s_ps_figure + lefse_18s_af_figure + lefse_18s_ps_figure + plot_layout(ncol = 2, nrow = 2, byrow = F)
ggsave("figure_field_lefse.png", plot = figure_field_lefse, width = 10, height =10, units = "cm", dpi = 600)
