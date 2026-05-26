## Plots for Evolving Strategies for Virus Discovery (Microbial Genomics, 2026)
# Amanda Araújo Serrão de Andrade, Andrea Silverj, Theodore Josephs, Ann C. Gregory*
# Corresponding author: ann.gregory@ucalgary.ca

library(tidyverse)
library(readxl) 
library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)

options(scipen = 999)
setwd("/Users/amanda.araujoserraod/Downloads/Review_Viral_Discovery_Microbial_Genomics/")

# Define Color palletes

my_colors <- c(
  "All Viruses" = "#7b9e87",
  "RNA Viruses" = "#d8b365",
  "ssDNA Viruses" = "#5f4b8b",
  "dsDNA" = "#e78ac3",
  "Hybrid" = "#2E7D32",
  "Ai-Based" = "#F28E2B",
  "Similarity-Based" = "#4E79A7",
  "Structure-Based" = "#E15759"
)

env_order <- c(
  "Extreme",
  "Freshwater",
  "Marine",
  "Soil",
  "Other Environments",
  "Host"
)

env_colors <- c(
  "Extreme" = "#5f4b8b",          # Light green (near #7b9e87)
  "Freshwater" = "#A8D5BA",       # Light yellow (near #d8b365)  
  "Marine" = "#4F6F8F",           # Light purple (near #5f4b8b)
  "Soil" = "#DFAF5F",             # Light pink (near #e78ac3)
  "Other Environments" = "gray", # Muted green-teal
  "Host" = "#E15759"              # Light peach/orange
)

host_colors <- c(
  "Human" = "#B2182B",              # Muted red
  "Rodents" = "#8C510A",            # Muted blue-green
  "Birds/Poultry" = "#F4A582",      # Muted orange
  "Pigs/Swine" = "#4F6F8F",         # Muted green
  "Bats" = "#9E9AC8",               # Muted purple
  "Insects/Arthropods" = "#7b9e87", # Muted yellow
  "Ruminants" = "#D8B365",          
  "Fish/Aquatic" = "#92C5DE",       # Muted teal
  "Other Mammals" = "#F1B6DA",      # Muted pink
  "Plants" = "#B3DE69",             # Muted mint
  "Other Hosts" = "#BDBDBD"         # Muted gray
)

# Inputs

tools <- fread("./SupplementaryTable3_Andrade_Silverj_Josephs_2026.csv", header=T)
glimpse(tools)
table(tools$`Based on:`)

similarity <- tools[tools$`Based on:` == "Similarity-Based", ]
table(similarity$`Which virus`)
table(similarity$Taxonomy)
table(similarity$Annotation)

ai <- tools[tools$`Based on:` == "AI-Based", ]
table(ai$`Which virus`)
table(ai$Taxonomy)
table(ai$Annotation)

hybrid <- tools[tools$`Based on:` == "Hybrid", ]
table(hybrid$`Which virus`)
table(hybrid$Taxonomy)
table(hybrid$Annotation)

# -------------------------------

database_table <- fread("./SupplementaryTable2_Andrade_Silverj_Josephs_2026.csv", header=T)

database_tools <- database_table %>%
  select(Authors, Year, `Viral_id_tools`) %>%
  filter(!is.na(`Viral_id_tools`)) %>%
  separate_rows(`Viral_id_tools`, sep = ";\\s*") %>%
  mutate(tool_clean = str_trim(tolower(`Viral_id_tools`)))

tools_df_clean <- tools %>%
  mutate(tool_clean = str_trim(tolower(Tool)))

database_tools_methods <- database_tools %>%
  left_join(tools_df_clean %>% select(tool_clean, Method = `Based on:`),
            by = "tool_clean", relationship = "many-to-many")

method_counts <- database_tools_methods %>%
  filter(!is.na(Method)) %>%
  group_by(Method) %>%
  summarise(n_databases = n(), .groups = "drop")

# Plot 1. Number of Citations (Fig3 D in the manuscript)

data_plot <- tools %>%
  mutate(
    Citations = as.numeric(`Number of Citations (Google Scholar)`),
    Category = `Based on:`
  ) %>%
  slice_max(Citations, n = 20) %>%   # 
  arrange(Citations)

plot <- ggplot(data_plot, aes(x = Citations, y = reorder(Tool, Citations), color = Category)) +
  geom_segment(aes(x = 1, xend = Citations, y = Tool, yend = Tool),
               linewidth = 0.8, color = "gray60") +
  geom_point(size = 3) +
  scale_color_manual(values = my_colors) +
  scale_x_log10() +
  theme_pubr() +
  labs(
    x = "Number of Citations (log scale)",
    y = NULL,
    color = NULL
  )

plot

ggsave("citations.svg",
       plot = plot,
       width = 4.5,
       height = 4.5,
       dpi = 400,
       bg = "transparent")


# Plot 2. Density plot with publication years (Fig 3C in the manuscript)

database_table_include <- database_table %>%
  mutate(vOTU_count = as.numeric(gsub("[^0-9]", "", `votu_count`))) %>%
  filter(!is.na(vOTU_count), !is.na(Year))

database_tools <- database_table_include %>%
  select(Authors, Year, votu_count,
         tools_raw = Viral_id_tools) %>%
  filter(!is.na(tools_raw)) %>%
  separate_rows(tools_raw, sep = ";\\s*") %>%
  mutate(tool_clean = str_trim(tolower(tools_raw)))

tools_df_clean <- tools %>%
  mutate(tool_clean = str_trim(tolower(Tool)))

database_tools_methods <- database_tools %>%
  left_join(tools_df_clean %>% select(tool_clean, Method = `Based on:`),
            by = "tool_clean", relationship = "many-to-many") %>%
  filter(!is.na(Method))

# Standardize method names
database_tools_methods <- database_tools_methods %>%
  mutate(Method = str_trim(Method),
         Method = str_to_title(Method))  

p <- ggplot(database_tools_methods,
            aes(x = Year, weight = votu_count, fill = Method)) +
  geom_density(alpha = 0.7, adjust = 1.5) +
  scale_fill_manual(values = my_colors) +
  theme_pubr() + 
  scale_x_continuous(
    breaks = seq(2012, 2026, by = 2),
    limits = c(2012, 2027),
    expand = expansion(mult = c(0.01, 0.02))
    ) +
  labs(
    x = "Publishing year",
    y = "Density of vOTUs",
    fill = NULL,
    title = NULL) + theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p
ggsave("density_tool_votus.svg",
       plot = p,
       width = 8,
       height = 2,
       dpi = 600,
       bg = "transparent")


## Violin plot (Fig 3B in the manuscript)

plot_table <- database_table_include %>%
  filter(!is.na(Sequencing_type)) %>%
  filter(!is.na(Data_Source)) %>%    # Remove any remaining NA seq_type
  mutate(
    seq_type = factor(Sequencing_type, levels = c("Virome-Enriched", "Bulk Sequencing", "Both Strategies")),
    Data_Source = factor(Data_Source, levels = c("Primary Studies", "Meta-Analysis"))
  )

plot <- ggplot(plot_table, aes(x = seq_type, y = votu_count, fill = Data_Source)) +
  geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, 
               position = position_dodge(width = 0.9)) +
  geom_jitter(alpha = 0.6, size = 1.5, position = position_jitterdodge(dodge.width = 0.9)) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(
    values = c("Primary Studies" = "#E15759", 
               "Meta-Analysis" = "#D9D9D9"), 
    name = NULL
  ) +
  labs(x = NULL, y = "vOTU Count (log scale)", title = NULL) +
  theme_pubr(base_size = 12) +
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
  )

plot

ggsave("violin_plot.svg",
       plot = plot,
       width = 4,
       height = 2.5,
       dpi = 600,
       bg = "transparent")

# Heatmap (Fig 3A in the manuscript)

custom_palette_original <- c(
  "#3F5F7F",     # 1. Muted navy
  "#4F6F8F",
  "#5f8f86",
  "#7b9e87",
  "#6B5E8C",  # Soft purple
  "#8C6B8C",
  "#c86a6a",
  "#F28E2B",
  "#DFAF5F",  # Soft orange-yellow
  "#F0D18A",   # Pale yellow
  "#F8E8AF"   # Pale yellow
)


env_df <- database_table_include %>%
  mutate(
    votu_count = as.numeric(gsub("[^0-9]", "", votu_count)),
    nucleotide_type = toupper(`Nucleotide type`)
  ) %>%
  separate_rows(nucleotide_type, sep = ";\\s*") %>%
  mutate(
    nucleotide_type = ifelse(is.na(nucleotide_type) | nucleotide_type == "",
                             "UNKNOWN",
                             nucleotide_type)
  ) %>%
  filter(
    !is.na(votu_count),
    !is.na(Year),
    !is.na(Environment_Group),
    nucleotide_type != "UNKNOWN"
  )

studies_per_env_year <- database_table_include %>%
  filter(!is.na(Year), !is.na(Environment_Group)) %>%
  count(Year, Environment_Group, name = "n_studies")

env_heatmap <- env_df %>%
  filter(Year >= 2000, Year <= 2026) %>%
  
  group_by(Year, Environment_Group, nucleotide_type) %>%
  summarise(
    total_votus = sum(votu_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  left_join(studies_per_env_year, by = c("Year", "Environment_Group")) %>%
  
  mutate(
    votus_norm = total_votus / n_studies,
    log_votus = log10(votus_norm + 1)
  ) %>%
  
  complete(
    Year, Environment_Group, nucleotide_type,
    fill = list(total_votus = 0, votus_norm = 0, log_votus = 0, n_studies = 0)
  )

plot <- ggplot(env_heatmap,
               aes(
                 x = Environment_Group,
                 y = factor(Year),
                 fill = log_votus
               )) +
  
  geom_tile(
    color = "white",
    linewidth = 0.1,
    width = 1,
    height = 1
  ) +
  
  scale_fill_gradientn(
    name = "vOTU abundance\n(normalized by studies per environment-year, log10)",
    colors = custom_palette_original,
    limits = c(0, max(env_heatmap$log_votus, na.rm = TRUE))
  ) +
  
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  
  facet_grid(~ nucleotide_type) +
  
  coord_equal(clip = "off") +
  
  theme_pubr(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "top"
  ) +
  
  labs(
    x = NULL,
    y = NULL
  )

plot

ggsave("heatmap_votu.svg",
       plot = plot,
       width = 4,
       height = 5,
       dpi = 600,
       bg = "transparent")


### Area plot for the environments (Fig 2 in the manuscript)

env_area <- database_table_include %>%
  filter(!is.na(Environment_Group), !is.na(Year)) %>%
  mutate(Environment_Group = fct_relevel(Environment_Group, env_order)) %>%
  count(Environment_Group, Year, name = "paper_count")

env_area <- env_area %>%
  group_by(Year, Environment_Group) %>%
  summarise(paper_count = sum(paper_count), .groups = "drop")

env_area_complete <- env_area %>%
  complete(
    Year = seq(2000, 2026, by = 1),
    Environment_Group,
    fill = list(paper_count = 0)
  )

fig2c <- ggplot(env_area,
                aes(x = Year, y = paper_count, fill = Environment_Group)) +
  geom_area(alpha = 0.8, position = "stack") +
  scale_fill_manual(values = env_colors) +
  scale_x_continuous(
    breaks = seq(2010, 2026, by = 2),
    limits = c(2010, 2026),
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Publishing year",
    y = "Number of studies",
    fill = "Environment"
  ) +
  theme(legend.position = "top")

fig2c

table(env_area$environment_group)
head(env_area)

ggsave("density_env_group.svg",
       plot = fig2c,
       width = 8,
       height = 3,
       dpi = 600,
       bg = "transparent")

host_data <- database_table %>% 
filter(Environment_Group == "Host") %>%
filter(!is.na(Year), !is.na(Environment)) %>%  
mutate(
    # Clean and group hosts into meaningful categories
    host_group = case_when(
      # HUMAN
      str_detect(tolower(Environment), "human|infant|centenarian|host|gut|dental|saliva|respiratory|vagina|children|cerebrospinal|oral|plasma") ~ "Human",
      
      # RODENTS (mouse, rat, rodent)
      str_detect(tolower(Environment), "mouse|mice|rat|rodent|shrew|marmot") ~ "Rodents",
      
      # BIRDS & POULTRY
      str_detect(tolower(Environment), "bird|birds|chicken|chickens") ~ "Birds/Poultry",
      
      # PIGS & SWINE
      str_detect(tolower(Environment), "pig|pigs|piglet|piglets|swine") ~ "Pigs/Swine",
      
      # BATS
      str_detect(tolower(Environment), "bat|bats") ~ "Bats",
      
      # INSECTS & ARTHROPODS
      str_detect(tolower(Environment), "insect|invertebrates|mosquito|mosquitos|tick|ticks|mites|bee|ant|ants|bee|crustacean|triatominae") ~ "Insects/Arthropods",
      
      # RUMINANTS (cow, rumen)
      str_detect(tolower(Environment), "rumen|cow|beef") ~ "Ruminants",
      
      # FISH & AQUATIC
      str_detect(tolower(Environment), "fish|salmon|eel|sharks|hydra|sponges") ~ "Fish/Aquatic",
      
      # OTHER MAMMALS
      str_detect(tolower(Environment), "lemurs|dog|cat|dog|rhesus|monkey|lemurs|gibbon|pangolin|rhesus|marsupial|mammal|herbivorous|gibbons") ~ "Other Mammals",
      
      # PLANTS
      str_detect(tolower(Environment), "plant|wheat|parsley|pepper|clover|birch leag|camelids|algae|plants") ~ "Plants",
      # GENERIC/VAGUE (gut, host, fecal, etc.)
      TRUE ~ "Other hosts"
    )
  ) %>%
count(host_group, Year, name = "study_count")
print("Host groups summary:")
print(host_data %>% count(host_group, wt = study_count, sort = TRUE))

fig2c <- ggplot(host_data, aes(x = Year, y = study_count, fill = host_group)) +
  geom_area(alpha = 0.8, position = "stack", color = "white", size = 0.2) +
  scale_fill_manual(values = host_colors) +
  scale_x_continuous(breaks = seq(2010, 2026, by = 1)) +
  theme_pubr(base_size = 12) +
  labs(x = "Publishing Year", y = "Number of Studies", 
       fill = NULL) +
  theme(legend.position = "none")

print(fig2c)

ggsave("Fig2C_hosts_area.svg", fig2c, 
       width = 8, height = 3, dpi = 600, bg = "white")

# Summarize data for Table 2

df <- database_table %>%
  mutate(year_bin = case_when(
    Year >= 2010 & Year <= 2013 ~ "2010-2013",
    Year >= 2014 & Year <= 2017 ~ "2014-2017",
    Year >= 2018 & Year <= 2021 ~ "2018-2021", 
    Year >= 2022 & Year <= 2025 ~ "2022-2025",
    Year == 2026 ~ "2026*",
    TRUE ~ NA_character_
  ))

bin1 <- df[df$year_bin == "2010-2013", ]
table(bin1$environment_group)
table(bin1$assembler)
table(bin1$viral_id_tools)
cat(sprintf("2022-2025: %d studies, Total vOTUs: %s, Mean±SD: %.1f±%.1f\n", 
            nrow(bin1), format(sum(bin1$votu_count, na.rm=T), big.mark=","), 
            mean(bin1$votu_count, na.rm=T), sd(bin1$votu_count, na.rm=T)))


bin2 <- df[df$year_bin == "2014-2017", ]
table(bin2$environment_group)
table(bin2$assembler)
table(bin2$viral_id_tools)
cat(sprintf("2022-2025: %d studies, Total vOTUs: %s, Mean±SD: %.1f±%.1f\n", 
            nrow(bin2), format(sum(bin2$votu_count, na.rm=T), big.mark=","), 
            mean(bin2$votu_count, na.rm=T), sd(bin2$votu_count, na.rm=T)))

bin3 <- df[df$year_bin == "2018-2021", ]
table(bin3$environment_group)
table(bin3$assembler)
table(bin3$viral_id_tools)
cat(sprintf("2022-2025: %d studies, Total vOTUs: %s, Mean±SD: %.1f±%.1f\n", 
            nrow(bin3), format(sum(bin3$votu_count, na.rm=T), big.mark=","), 
            mean(bin3$votu_count, na.rm=T), sd(bin3$votu_count, na.rm=T)))

bin4 <- df[df$year_bin == "2022-2025", ]
table(bin4$environment_group)
table(bin4$assembler)
table(bin4$viral_id_tools)
cat(sprintf("2022-2025: %d studies, Total vOTUs: %s, Mean±SD: %.1f±%.1f\n", 
            nrow(bin4), format(sum(bin4$votu_count, na.rm=T), big.mark=","), 
            mean(bin4$votu_count, na.rm=T), sd(bin4$votu_count, na.rm=T)))

bin5 <- df[df$year_bin == "2026*", ]
table(bin5$environment_group)
table(bin5$assembler)
table(bin5$viral_id_tools)
cat(sprintf("2022-2025: %d studies, Total vOTUs: %s, Mean±SD: %.1f±%.1f\n", 
            nrow(bin5), format(sum(bin5$votu_count, na.rm=T), big.mark=","), 
            mean(bin5$votu_count, na.rm=T), sd(bin5$votu_count, na.rm=T)))
