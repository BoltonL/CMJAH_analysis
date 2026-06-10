##Copyright: Larisse Bolton, 2026
##This script conducts a WISCA analysis on early-onset bloodstream infection data from Charlotte Maxeke Johannesburg Hospital
##Data custodian: Vindana Chibabhai (NICD)
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


library(openxlsx)
library(tidyverse)
library(lubridate)
library(AMR)

##source scripts from repository: 
source("./CMJAH_analysis/Code/parameters_amr.R")
source("./CMJAH_analysis/Code/wisca_AMR.R")
source("./CMJAH_analysis/Code/plot_wisca_nicd.R")

#--------------------------------------------------------------------------------------------------------------------------------------
#Read in data and address variable types
#--------------------------------------------------------------------------------------------------------------------------------------
data_path_v2 <- readRDS("./CMJAH_analysis/Data/gram_stain_ref.RDS") #generate gram stain reference set

data_dir <- "./wisca_data/data"
data_name <- "Larisse_CMJAH Neonatal BC data_complete_clean_25082024.xlsx"
data_raw <- read.xlsx(xlsxFile = file.path(data_dir, data_name), na.strings = c("NULL","NA"), sheet = 1)
data_raw <- select(data_raw,-ORGANISM_COMMENTS, -GROWTH_COMMENTS)
data_raw <- data_raw %>%
  arrange(REGISTRATION_DATE) %>%
  mutate(PATIENT.NUMBER = unlist(str_extract_all(PATIENT.NUMBER,"[0-9]+")),
         DATE_OF_BIRTH = janitor::convert_to_date(DATE_OF_BIRTH),
         REGISTRATION_DATE = ymd(REGISTRATION_DATE),
         ORGANISM_NAME = ifelse(ORGANISM_NAME == "PICHIA KUDRIAVZEVII","CANDIDA KRUSEI",ORGANISM_NAME),
         ORGANISM_NAME = ifelse(ORGANISM_NAME == "NAKASEOMYCES GLABRATA", "CANDIDA GLABRATA", ORGANISM_NAME),
         ORGANISM_NAME = ifelse(ORGANISM_NAME == "CLAVISPORA LUSITANIAE", "CANDIDA LUSITANIAE", ORGANISM_NAME),
         PATIENT.NUMBER = as.factor(PATIENT.NUMBER),
         GENDER = as.factor(GENDER),
         EPISODE_NO = as.factor(EPISODE_NO)) %>%
  mutate(DAY_ONSET = REGISTRATION_DATE - DATE_OF_BIRTH,
         infection_type = ifelse(DAY_ONSET < 3, "EONS","HAI")) %>%
  mutate(infection_type = as.factor(infection_type))
names(data_raw)[which(names(data_raw) == "MERO.+.VANOC")] <- "MERO.+.VANCO"
data_raw$GENDER <- recode_factor(data_raw$GENDER,
                                 `F` = "FEMALE",
                                 `M` = "MALE")
data_raw <- rename(data_raw, 
                   "PIPTAZ.+.AMIK" = "TAZ.+.AMIK")

exclude_candida_jhb <- data_raw %>%
  select(ORGANISM_NAME) %>%
  filter(str_detect(ORGANISM_NAME, "^CANDIDA")) %>%
  mutate(ORGANISM_NAME = str_to_title(ORGANISM_NAME))
exclude_path <- c("(Unknown Name)","Unknown", exclude_candida_jhb$ORGANISM_NAME)

infection <- levels(data_raw$infection_type)


#------------------------------------------------------------------------------------------------------------------------------------------
# Preliminary analysis
#-----------------------------------------------------------------------------------------------------------------------------------------
data_raw_episode <- data_raw %>% # an episode is defined as a first isolate within 14 days and unique even if more than 1 organism isolated 
  filter(!duplicated(EPISODE_NO))

# 1. Syndrome proportion
data_infection_type <- data.frame(table(data_raw_episode$infection_type)) #episodes-level
names(data_infection_type) <- c("syndrome","n")
data_infection_type <- data_infection_type %>%
  mutate(perc.n = (n/sum(n))*100)

data_infection_type_path <- data.frame(table(data_raw$infection_type)) #pathogen-level
names(data_infection_type_path) <- c("syndrome","n")
data_infection_type_path <- data_infection_type_path %>%
  mutate(perc.n = (n/sum(n))*100)

# 2. Microbial proportion

polymicro <- function(y){
  if (group_size(y)>1){
    y %>% mutate(POLYMICRO = "POLYMICROBIAL")
  } else {
    y %>% mutate(POLYMICRO = "MONOMICROBIAL")
  }
}

data_raw_micro <- data_raw %>%
  filter(infection_type == "HAI") %>%
  group_by(PATIENT.NUMBER,EPISODE_NO) %>%
  do(polymicro(.)) %>%
  ungroup()
  
data_micro_type_path <- data.frame(table(data_raw_micro$POLYMICRO)) #pathogen-level
names(data_micro_type_path) <- c("poly/mono","n")
data_micro_type_path <- data_micro_type_path %>%
  mutate(perc.n = (n/sum(n))*100)

data_micro_type <- data_raw_micro %>%
  filter(!duplicated(EPISODE_NO))
data_micro_type_epi <-  data.frame(table(data_micro_type$POLYMICRO)) #episode-level
names(data_micro_type_epi) <- c("poly/mono","n")
data_micro_type_epi <- data_micro_type_epi %>%
  mutate(perc.n = (n/sum(n))*100)

# 1. Early-onset bloodstream infections
data_eons_jhb <- data_raw %>%
  filter(infection_type == "EONS") %>%
  mutate(HOSPITAL = as.factor("CMJAH ")) %>%
  select(PATIENT.NUMBER, DATE_OF_BIRTH, GENDER, REGISTRATION_DATE, ORGANISM_NAME, 
         `AMPI.+.GENTA`, `AMPI.+.AMIK`, `PIPTAZ.+.AMIK`, infection_type, HOSPITAL)
#saveRDS(data_eons_jhb,"./Data/data_eons_jhb.RDS")

# 2. Hospital-acquired infections
data_hai_jhb <- data_raw %>%
  filter(infection_type == "HAI") 

## Pathogen distribution
data_hai_jhb_pathogens <- data_hai_jhb %>%
  group_by(ORGANISM_NAME) %>%
  count() %>%
  ungroup() %>%
  arrange(desc(n)) 
  
data_hai_jhb_pathogens_other <- data_hai_jhb_pathogens %>%
  filter(n <= 40) %>%
  summarise(other_paths = sum(n)) %>%
  bind_cols(c("OTHER")) 
names(data_hai_jhb_pathogens_other) <- c("n","ORGANISM_NAME")
data_hai_jhb_pathogens_other <- select(data_hai_jhb_pathogens_other,"ORGANISM_NAME","n")

data_hai_jhb_pathogens_main <- data_hai_jhb_pathogens %>%
  filter(n > 40) %>%
  bind_rows(data_hai_jhb_pathogens_other) %>%
  mutate(perc.n = (n/sum(n))*100) %>%
  mutate(ORGANISM_NAME = str_to_sentence(ORGANISM_NAME)) %>%
  arrange(desc(perc.n))

pie_path <-ggplot(data_hai_jhb_pathogens_other_main, aes(x = "",y = perc.n, fill = fct_inorder(ORGANISM_NAME))) +
  geom_bar(stat="identity", width = 1, color = "white") +
  geom_text(
    aes(label = paste0(round(perc.n,0),"%")),
    position = position_stack(vjust = 0.5),
    fontface = "bold"
  ) +
  labs(fill = "") +
  coord_polar(theta = "y", start = 0) +
  scale_fill_brewer(palette = "Set1") +
  theme_void() + 
  theme(
    legend.text = element_text(face = "italic"),
    legend.title = element_text(face = "italic"),
    legend.position = "bottom",
    axis.text.y = element_text(size = 16))

data_hai_jhb_pathogens_other_out <- data_hai_jhb_pathogens %>%
  filter(n <= 10) %>%
  mutate(perc.n = (n/sum(n))*100) %>%
  summarise(other_paths = sum(n)) %>%
  bind_cols(c("Other"))
names(data_hai_jhb_pathogens_other_out) <- c("n","ORGANISM_NAME")
data_hai_jhb_pathogens_other_out <- select(data_hai_jhb_pathogens_other_out,"ORGANISM_NAME","n")

data_hai_jhb_pathogens_other_main <- data_hai_jhb_pathogens %>%
  filter(n > 10, n <= 40) %>%
  bind_rows(data_hai_jhb_pathogens_other_out) %>%
  mutate(perc.n = (n/sum(n))*100)  %>%
  mutate(ORGANISM_NAME = str_to_sentence(ORGANISM_NAME))

donut_path <-ggplot(data_hai_jhb_pathogens_other_main,
       aes(x = 2, y = perc.n, fill = fct_inorder(ORGANISM_NAME))) +
  geom_col(color = "white") +
  geom_text(
    aes(label = paste0(round(perc.n,0),"%")),
    position = position_stack(vjust = 0.5),
    fontface = "bold"
  ) +
  labs(fill = "") +
  coord_polar(theta = "y", start = 0) +
  xlim(c(0.5, 2.5)) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_void() +
  theme(
    legend.text = element_text(face = "italic"),
    legend.title = element_text(face = "italic"),
    legend.position = "top",
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

ggpubr::ggarrange(pie_path, donut_path, ncol = 2, nrow = 1)

#saveRDS(data_hai_jhb,"./Data/data_hai_jhb.RDS")


antibiotic_eons <- c("AMP+GENT","AMP+AMIK","PIPTAZ+AMIK")
antibiotic_hai <- c("PIPTAZ+AMIK","MEM","MEM+VANCO","MEM+COL","MEM+AMIK")

include_pathogens <- 10
include_pathogens_hai <- 9

#------------------------------------------------------------------------------------------------------------------------------
# WISCA 
#------------------------------------------------------------------------------------------------------------------------------
deduped <- "yes"
analysis_level <- "none"
susceptible_I_in <- "S"
syndrome_name_eons <- "Early-onset bloodstream infections"
syndrome_name_hai <- "Hospital-associated bloodstream infections"

#Generate wisca input parameters -EONS
wisca_parameters_eons <- wisca_params(x = data_eons_jhb, 
                                      antibiotic_in = antibiotic_eons, 
                                      pathogen_in = include_pathogens,
                                      analysis = analysis_level,
                                      exclude = exclude_path, 
                                      isolate_first = deduped,
                                      susceptible_I = susceptible_I_in,
                                      infection_in = "EONS",
                                      infection_full = syndrome_name_eons)

plot_bugs_eons <- wisca_parameters_eons %>%
  #  group_by(hospital, mo) %>%
  mutate(perc.n = prop.n*100) %>%
  distinct(mo, .keep_all = TRUE) 
# %>%
# ungroup()

#Pathogen incidence plot
ggplot(plot_bugs_eons, aes(x = fct_inorder(fullname))) + geom_col(aes(y = perc.n)) +
  labs(title = paste0("Percentage contribution of bacterial pathogens to neonatal early-onset bloodstream infections"), x = "Pathogen classification",
       y = paste0("Percentage of top ",include_pathogens," pathogens")) + coord_flip() +
  theme_bw() + theme(axis.text.y = element_text(face = "italic"), text = element_text(size = 14)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(plot_bugs_eons$perc.n))))

#Generate wisca input parameters - HAI
wisca_parameters_hai <- wisca_params(x = data_hai_jhb, 
                                     antibiotic_in = antibiotic_hai, 
                                     pathogen_in = include_pathogens_hai,
                                     analysis = analysis_level,
                                     exclude = exclude_path, 
                                     isolate_first = deduped,
                                     susceptible_I = susceptible_I_in,
                                     infection_in = "HAI",
                                     syndrome_name_hai <- "Hospital-associated bloodstream infections")

#Pathogen incidence plot
plot_bugs_hai <- wisca_parameters_hai %>%
  mutate(mo = as.mo(mo)) %>%
  mutate(perc.n = prop.n*100) %>%
  distinct(mo, .keep_all = TRUE) %>%
  left_join(data_path_v2[,c("mo","GRAM1")], by = "mo") %>%
  distinct(mo, .keep_all = TRUE) %>%
  arrange(prop.n) %>%
  mutate(fullname = str_to_sentence(fullname)) %>%
  mutate(perc.n = prop.n*100) %>%
  mutate(GRAM1 = if_else(GRAM1 == "POS","Positive",GRAM1)) %>%
  mutate(GRAM1 = if_else(GRAM1 == "NEG","Negative",GRAM1)) %>%
  mutate(GRAM1 = if_else(GRAM1 == "FUNG","Fungi",GRAM1)) %>%
  mutate(GRAM1 = factor(GRAM1, levels = c("Positive","Negative","Fungi")))

# %>%
# ungroup()

ggplot(plot_bugs_hai, aes(x = fct_inorder(fullname))) + geom_col(aes(y = perc.n, fill = GRAM1)) +
  labs(title = paste0("Percentage contribution of bacterial pathogens to neonatal hospital-associated bloodstream infections"), x = "Pathogen classification",
       y = paste0("Percentage of top ",include_pathogens_hai," pathogens"), fill = "Gram stain") + coord_flip() +
  theme_bw() + theme(axis.text.y = element_text(face = "italic"), text = element_text(size = 14)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(max(plot_bugs_hai$perc.n))))

#Coverage estimate calculations
wisca_cover_eons <- wisca_funct(wisca_in = wisca_parameters_eons, analysis = analysis_level)
wisca_cover_hai <- wisca_funct(wisca_in = wisca_parameters_hai, analysis = analysis_level)

abo_first_level <- "MEM+VAN"
wisca_cover_hai_upd <- wisca_cover_hai %>%
  mutate(Regimen_full = if_else(Regimen_full == "Piperacillin/tazobactam + Amikacin","Piperacillin-tazobactam + Amikacin",Regimen_full)) %>%
  filter(Regimen != abo_first_level) %>%
  arrange(Coverage)
wisca_cover_plot_in <- bind_rows(subset(wisca_cover_hai, Regimen == abo_first_level),
                                 wisca_cover_hai_upd)
wisca_cover_plot_in$Regimen_full <- fct_inorder(wisca_cover_plot_in$Regimen_full)

#Output plots and tables
wisca_plot_out_eons <- wisca_plot(params_plot = wisca_parameters_eons,
                                  cover_plot = wisca_cover_eons, 
                                  pathogen_in = include_pathogens,
                                  analysis = analysis_level,
                                  infection_full = syndrome_name_eons)

wisca_plot_out_hai <- wisca_plot(params_plot = wisca_parameters_hai,
                                 cover_plot = wisca_cover_plot_in, 
                                 pathogen_in = include_pathogens_hai,
                                 analysis = analysis_level,
                                 infection_full = syndrome_name_hai)

#publication table output
wisca_eons_table <- wisca_plot_out_eons %>%
  select(`Antibiotic regimens`,`Microorganism fullname`:`Percentage pathogen sensitive`)
write.xlsx(wisca_eons_table,"./analysis/consult/Vindana/wisca_eons_table_out.xlsx")

wisca_hai_table <- wisca_plot_out_hai %>%
  select(`Antibiotic regimens`,`Microorganism fullname`:`Percentage pathogen sensitive`)
write.xlsx(wisca_hai_table,"./analysis/consult/Vindana/wisca_hai_table_out.xlsx")


