# Load libraries

library(tidyverse)
library(tidymodels)
library(ssdtools)

# Define colors

colors <- c("#5F7800", "#EB8200", "#00A0BE", "#AAB400", "#FFB400", "#82C8DC")

# Read files

es <- read_csv("es.csv")
bcpc <- read_csv("bcpc.csv")
trex <- read_csv("trex.csv")
erc <- read_csv("erc.csv")

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

fish_acute_pred %>%
  filter(percent == 58) %>%
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

### Select data

rbt_daphnia_acute <- bcpc %>%
  mutate(result = log10(result)) %>%
  group_by(ai, group, species) %>%
  summarize(result = mean(result)) %>%
  mutate(result = 10^result) %>%
  ungroup() %>%
  spread(species, result) %>%
  na.omit()

### Fit linear model

rbt_daphnia_acute_fit <- rbt_daphnia_acute %>%
  mutate(across(rainbow_trout:daphnia, log10)) %>%
  lm(rainbow_trout ~ daphnia, data = .)

### Summarize model

rbt_daphnia_acute_fit %>% summary()

### Generate prediction for rainbow trout

rbt_daphnia_acute_fit %>%
  predict(., data.frame(daphnia = c(log10(9.6))), interval = "confidence") %>%
  as.double() %>%
  10^.

### Plot linear relationship

rbt_daphnia_acute %>%
  ggplot(aes(daphnia, rainbow_trout)) +
  geom_smooth(method = "lm", color = "black", alpha = 0.2, parse = TRUE) +
  geom_point(aes(color = group), size = 2, alpha = 0.6) +
  theme_bw() +
  coord_fixed() +
  scale_x_log10(limits = c(0.01, 1000), breaks = c(0.01, 0.1, 1, 10, 100, 1000),
                labels = scales::label_comma(drop0trailing = TRUE)) +
  scale_y_log10(limits = c(0.01, 1000), breaks = c(0.01, 0.1, 1, 10, 100, 1000),
                labels = scales::label_comma(drop0trailing = TRUE)) +
  scale_color_manual(name = "Group", values = colors) +
  xlab("48-h Daphnia EC50 (mg/L)") +
  ylab("96-h Rainbow trout LC50 (mg/L)") +
  geom_segment(aes(x = 9.6, xend = 9.6, y = 0, yend = 2.68),
               linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 9.6, y = 2.68, yend = 2.68),
             linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 9.6, y = 2.68, yend = 2.68),
               linetype = "dashed") +
  annotate("text", x = 0.02, y = 750, size = 4, hjust = 0,
           label = "log10(y) = -0.330 + 0.772 * log10(x)") +
  annotate("text", x = 0.02, y = 450, size = 4, hjust = 0, parse = TRUE,
           label = "R^2 == 0.686")

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

ggsave("fig_5.png", width = 6.5, height = 3, units = "in", dpi = 300, scale = 1.4)