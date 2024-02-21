# Load libraries

library(tidyverse)
library(tidymodels)
library(ssdtools)

# Define colors

colors <- c("#5F7800", "#EB8200", "#00A0BE", "#AAB400", "#FFB400", "#82C8DC")

# Read files

es <- read_csv("es.csv")
envirotox <- read_csv("envirotox.csv")
trex <- read_csv("trex.csv")
erc <- read_csv("erc.csv")

casrn_pesticides <- read_lines("casrn_pesticides.txt")
casrn_accase <- read_lines("casrn_accase.txt")

# Aquatic risk assessment

## Fish acute chemical toxicity distribution

### Select data

fish_acute <- es %>%
  filter(study == "Fish acute") %>%
  arrange(ai) %>%
  select(ai, group, species, result) %>%
  group_by(ai, group) %>%
  mutate(Conc = min(result),
         max = max(result)) %>%
  ungroup() %>%
  spread(species, result)

### Fit distributions

fish_acute_fit <- fish_acute %>%
  ssd_fit_dists()

### Model averaging and prediction

fish_acute_pred <- predict(fish_acute_fit, ci = FALSE) %>%
  mutate(type = "Fish acute (LC50)")

### Estimate hazard concentration at 5th percentile (HC5)

fish_acute_hc5 <- fish_acute_pred %>%
  filter(percent == 5) %>%
  pull(est)

### Tidy data for plotting

fish_acute_plot <- fish_acute %>%
  bind_cols(data.frame(pr = ssd_ecd_data(.))) %>%
  mutate(type = "Fish acute (LC50)") %>%
  gather(`Freshwater (cold)`:Saltwater, key = "species", value = result) %>%
  mutate(species = recode(species, "Freshwater (cold)" = "Freshwater",
                          "Freshwater (warm)" = "Freshwater")) %>%
  select(ai, group, species, type, result, pr, max)

## Fish chronic chemical toxicity distribution

### Select data

fish_chronic <- es %>%
  filter(study == "Fish early life stage") %>%
  arrange(ai) %>%
  select(ai, group, species, result) %>%
  group_by(ai, group) %>%
  mutate(Conc = min(result),
         max = max(result)) %>%
  ungroup() %>%
  spread(species, result)

### Fit distributions

fish_chronic_fit <- fish_chronic %>%
  ssd_fit_dists()

### Model averaging and prediction

fish_chronic_pred <- predict(fish_chronic_fit, ci = FALSE) %>%
  mutate(type = "Fish chronic (NOAEC)")

### Estimate hazard concentration at 5th percentile (HC5)

fish_chronic_hc5 <- fish_chronic_pred %>%
  filter(percent == 5) %>%
  pull(est)

### Tidy data for plotting

fish_chronic_plot <- fish_chronic %>%
  bind_cols(data.frame(pr = ssd_ecd_data(.))) %>%
  mutate(type = "Fish chronic (NOAEC)") %>%
  gather(Freshwater:Saltwater, key = "species", value = result) %>%
  select(ai, group, species, type, result, pr, max)

## Designate estimated environmental concentrations (EECs)

eec <- data.frame(type = c("Fish acute (LC50)", "Fish chronic (NOAEC)"),
                  eec = c(0.001468, 0.001396),
                  label = c("4-d EEC = 0.001468", "60-d EEC = 0.001396"))

## Combine aquatic data into one plot

bind_rows(fish_acute_plot, fish_chronic_plot) %>%
  na.omit() %>%
  ggplot() +
  geom_line(data = bind_rows(fish_acute_pred, fish_chronic_pred),
            aes(est, percent/100)) +
  geom_text(data = . %>% select(ai, type, pr, max) %>% unique(),
            aes(x = max * 1.6, y = pr, label = ai), size = 3, alpha = 0.6,
            hjust = 0) +
  geom_point(aes(result, pr, color = group, shape = species), size = 2,
             alpha = 0.8) +
  facet_wrap(~type, ncol = 1) +
  geom_vline(data = eec, aes(xintercept = eec),
             linetype = "dashed") +
  geom_label(data = eec, aes(x = eec, y = 0.95, label = label),
            hjust = 0.3) +
  scale_x_log10(labels = scales::label_comma(drop0trailing = TRUE),
                breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                limits = c(0.0005, 3000)) +
  scale_y_continuous(labels = label_percent(), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(name = "Group", values = colors) +
  scale_shape_discrete(name = "Species") +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin()) +
  xlab("Concentration (mg/L)") +
  ylab("Percent rank") +
  guides(color = guide_legend(order = 1), 
         shape = guide_legend(order = 2))

ggsave("fig_1.png", width = 6.5, height = 8, units = "in", dpi = 300)

## Fish vs daphnia

### Select daphnia and rainbow trout data for pesticides

daphnia <- envirotox %>%
  filter(CAS %in% casrn_pesticides,
         `Effect is 5X above water solubility` == 0,
         `Latin name` == "Daphnia magna",
         `Duration (hours)` == 48,
         `Test statistic` == "EC50") %>%
  mutate(result_log10 = log10(`Effect value`)) %>%
  group_by(CAS, `Chemical name`) %>%
  summarize(daphnia = mean(result_log10)) %>%
  ungroup()

rbt <- envirotox %>%
  filter(CAS %in% casrn_pesticides,
         `Effect is 5X above water solubility` == 0,
         `Latin name` == "Oncorhynchus mykiss",
         `Duration (hours)` == 96,
         `Test statistic` == "LC50") %>%
  mutate(result_log10 = log10(`Effect value`)) %>%
  group_by(CAS, `Chemical name`) %>%
  summarize(rbt = mean(result_log10)) %>%
  ungroup()

rbt_daphnia <- daphnia %>%
  left_join(rbt, by = join_by(CAS, `Chemical name`)) %>%
  mutate(`ACCase inhibitor` = ifelse(CAS %in% casrn_accase, TRUE, FALSE)) %>%
  na.omit()

### Construct linear model

rbt_daphnia_lm <- rbt_daphnia %>%
  lm(rbt ~ daphnia, data = .)

rbt_daphnia_lm %>% summary()

### Generate prediction for new ACCase inhibitor

rbt_pred <- rbt_daphnia_lm %>%
  predict(., data.frame(daphnia = c(log10(9.6))), interval = "prediction") %>%
  as.double() %>%
  10^.

### Compare estimate with fish acute chemical toxicity distribution

fish_acute_pred %>%
  filter(percent == 69) %>%
  pull(est)

### Plot data

rbt_daphnia %>%
  mutate(across(daphnia:rbt, ~ 10^.)) %>%
  ggplot(aes(daphnia, rbt)) +
  geom_point(aes(color = `ACCase inhibitor`), size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  theme_bw() +
  coord_fixed() +
  scale_x_log10(limits = c(0.00001, 10000), breaks = c(0.001, 0.1, 10, 1000),
                labels = scales::label_comma(drop0trailing = TRUE)) +
  scale_y_log10(limits = c(0.00001, 10000), breaks = c(0.001, 0.1, 10, 1000),
                labels = scales::label_comma(drop0trailing = TRUE)) +
  scale_color_manual(name = "ACCase inhibitor", values = colors[c(3,2)]) +
  xlab("48-h Daphnia EC50 (mg/L)") +
  ylab("96-h Rainbow trout LC50 (mg/L)") +
  geom_segment(aes(x = 9.6, xend = 9.6, y = 0, yend = rbt_pred[1]),
               linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 9.6, y = rbt_pred[1], yend = rbt_pred[1]),
             linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 9.6, y = rbt_pred[1], yend = rbt_pred[1]),
               linetype = "dashed") +
  annotate("text", x = 0.00002, y = 2000, size = 4, hjust = 0,
           label = "log10(y) = 0.197 + 0.613 * log10(x)") +
  annotate("text", x = 0.00002, y = 500, size = 4, hjust = 0, parse = TRUE,
           label = "R^2 == 0.543")

ggsave("fig_2.png", width = 6.5, height = 5.5, units = "in", dpi = 300)

# Terrestrial risk assessment

## Avian chronic chemical toxicity distribution

### Select data

avian_chronic <- es %>%
  filter(study == "Avian chronic",
         sign == "=") %>%
  arrange(ai) %>%
  select(ai, group, species, result) %>%
  group_by(ai, group) %>%
  mutate(Conc = min(result),
         max = max(result)) %>%
  ungroup() %>%
  spread(species, result)

### Fit distributions

avian_chronic_fit <- avian_chronic %>%
  ssd_fit_dists()

### Model averaging and prediction

avian_chronic_pred <- predict(avian_chronic_fit, ci = FALSE) %>%
  mutate(type = "Avian chronic (NOAEC)")

### Estimate hazard concentration at 5th percentile (HC5)

avian_chronic_hc5 <- avian_chronic_pred %>%
  filter(percent == 5) %>%
  pull(est)

### Tidy data for plotting

avian_chronic_plot <- avian_chronic %>%
  bind_cols(data.frame(pr = ssd_ecd_data(.))) %>%
  mutate(type = "Avian chronic (NOAEC)") %>%
  gather(`Mallard duck`:`Northern bobwhite`, key = "species",
         value = result) %>%
  select(ai, group, species, type, result, pr, max)

## Mammalian chronic

### Select data

mammalian_chronic <- es %>%
  filter(study == "Mammalian chronic",
         sign == "=") %>%
  arrange(ai) %>%
  mutate(species = "Mammal") %>%
  select(ai, group, species, result) %>%
  group_by(ai, group) %>%
  mutate(Conc = min(result),
         max = max(result)) %>%
  ungroup() %>%
  spread(species, result)

### Fit distributions

mammalian_chronic_fit <- mammalian_chronic %>%
  ssd_fit_dists()

### Model averaging and predictions

mammalian_chronic_pred <- predict(mammalian_chronic_fit, ci = FALSE) %>%
  mutate(type = "Mammalian chronic (NOAEC)")

### Estimate hazard concentration at 5th percentile (HC5)

mammalian_chronic_hc5 <- mammalian_chronic_pred %>%
  filter(percent == 5) %>%
  pull(est)

### Tidy data for plotting

mammalian_chronic_plot <- mammalian_chronic %>%
  bind_cols(data.frame(pr = ssd_ecd_data(.))) %>%
  mutate(type = "Mammalian chronic (NOAEC)") %>%
  gather(Mammal, key = "species", value = result) %>%
  select(ai, group, species, type, result, pr, max)

## Combine terrestrial data into one plot

bind_rows(avian_chronic_plot, mammalian_chronic_plot) %>%
  na.omit() %>%
  mutate(species = factor(species, levels = c("Mallard duck",
                                              "Northern bobwhite",
                                             "Mammal"))) %>%
  ggplot() +
  geom_line(data = bind_rows(avian_chronic_pred, mammalian_chronic_pred),
            aes(est, percent/100)) +
  geom_text(data = . %>% select(ai, type, pr, max) %>% unique(),
            aes(x = max*1.2, y = pr, label = ai), size = 3, alpha = 0.6,
            hjust = 0) +
  geom_point(aes(result, pr, color = group, shape = species), size = 2,
             alpha = 0.8) +
  facet_wrap(~type, ncol = 1) +
  scale_x_log10(labels = scales::label_comma(drop0trailing = TRUE),
                limits = c(5, 3500),
                breaks = c(1, 10, 100, 1000)) +
  scale_y_continuous(labels = label_percent(), breaks = seq(0, 1, by = 0.2)) +
  scale_shape_discrete(name = "Species") +
  scale_color_manual(name = "Group", values = colors) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box="vertical",
        legend.margin=margin()) +
  xlab("Concentration (mg/kg-diet)") +
  ylab("Percent rank")

ggsave("fig_3.png", width = 6.5, height = 8, units = "in", dpi = 300)

## Select terrestrial exposure decline data

trex_decline <- trex %>%
  mutate(item = factor(item, levels = c("Short grass", "Broadleaf plants",
                                       "Tall grass", "Arthropods",
                                       "Fruits/pods"))) %>%
  gather(upper:mean, key = "type", value = "value") %>%
  mutate(type = recode(type, "upper" = "Upper Kenaga",
                       "mean" = "Mean Kenaga"),
         type = factor(type, levels = c("Upper Kenaga", "Mean Kenaga")))

## Provide annotations for plotting

trex_decline_ann <- data.frame(type = rep("Upper Kenaga", 2),
                               hc5 = c(avian_chronic_hc5,
                                       mammalian_chronic_hc5),
                               label = c("Avian HC5", "Mammalian HC5"))

## Plot terrestrial exposure declined data with HC5s

trex_decline %>%
  ggplot(aes(days, value)) +
  geom_line(aes(color = item)) +
  geom_hline(yintercept = mammalian_chronic_hc5, linetype = "dotted") +
  geom_hline(yintercept = avian_chronic_hc5, linetype = "dashed") +
  geom_label(data = trex_decline_ann, aes(x = 280, y  = hc5, label = label)) +
  facet_wrap(~type) +
  theme_bw() +
  scale_color_manual(values = colors, name = "Dietary item") +
  xlab("Days") +
  ylab("Concentration (mg/kg-diet)")

ggsave("fig_4.png", width = 6.5, height = 4, units = "in", dpi = 300,
       scale = 1.2)

## Prepare data for ecotoxicity risk calculator

### Generate terrestrial exposure distributions at 0 and 70 days based on
### upper and mean Kenaga exposure estimates

set.seed(100)

trex %>%
  filter(item == "Short grass",
         days %in% c(0, 70)) %>%
  mutate(across(c(upper, mean), log10),
         sd = (upper-mean)/1.96,
         eec = map2(mean, sd, rnorm, n = 10000)) %>%
  unnest(cols = c(eec)) %>%
  select(item, days, eec) %>%
  mutate(eec = 10^eec) %>%
  arrange(days, eec) %>%
  write_csv("erc_exposure.csv")

### Select mammalian chronic data

mammalian_chronic %>%
  select(ai, endpoint = Conc) %>%
  arrange(endpoint) %>%
  write_csv("erc_effects.csv")

### Data used to generate risk curves with ecotoxicity risk calculator
### Curves extracted and plotted below

erc %>%
  mutate(across(probability_effect:probability_exceedance, ~ .x / 100),
         type = factor(type, levels = c("Risk curve",
                                        "De minimis-low risk boundary",
                                        "Low-intermediate risk boundary",
                                        "Intermediate-high risk boundary"))) %>%
ggplot(aes(probability_effect, probability_exceedance)) +
  geom_line(aes(color = type, linetype = type, linewidth = type)) +
  facet_wrap(~time) +
  scale_color_manual(name = "Legend", values = c(colors[3], rep("Gray60", 3))) +
  scale_linetype_manual(name = "Legend", values = c("solid", "dotted",
                                                    "dashed", "longdash")) +
  scale_linewidth_manual(name = "Legend", values = c(1, 0.5, 0.5, 0.5)) +
  theme_bw() +
  coord_fixed() +
  scale_y_continuous(labels = label_percent(), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(labels = label_percent(), breaks = seq(0, 1, by = 0.2)) +
  xlab("Probability of effects") +
  ylab("Probability of exposure")

ggsave("fig_5.png", width = 6.5, height = 3, units = "in", dpi = 300,
       scale = 1.4)