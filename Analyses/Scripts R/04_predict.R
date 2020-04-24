library(tidyverse)
library(magrittr)
# library(gganimate)
# library(transformr)
library(rstan)
library(loo)
library(progress)
library(furrr)
library(patchwork)
plan(multiprocess)
options(mc.cores = ifelse(availableCores() > 4,
                          10,
                          availableCores()))
filter  <- dplyr::filter
lag     <- dplyr::lag
extract <- rstan::extract

theme_set(theme_bw())

source("00_functions.R")

## ---------------- Loading the models

print(str_glue("Loading the models"))

# Loading meta-data
mods <- readRDS("../STAN_mods.rds")
ids  <- readRDS("../FinalIDs.rds")


## ---------------- Computing predicted values

# Computing the predicted values
mods <-
    mods %>%
    mutate(FitFunc = future_pmap(list(ID, Dist, Model_Year),
                                 load_and_generate_fitfunc,
                                 progress = FALSE,
                                 .progress = TRUE))

## ---------------- Outputting some graphs

## Selecting the best models
loo  <- readRDS("../loo.rds")
best <- inner_join(mods,
                   loo %>% 
                       filter(Delta_LOOIC == 0) %>% 
                       select(-LOOIC, -Delta_LOOIC))

## Plotting the graphs for the best model
restart_pb(nrow(best))
best %$%
    pmap(list(Data, FitFunc, DataID, Species, Population),
         ~ {pb$tick();
            ggsave(plot = plot_fitfunc(..1, ..2),
                   filename = str_glue("../../Figures/Fitted_STAN_best/{..3}_{..4}_{..5}.pdf"),
                   device   = cairo_pdf(),
                   width    = 15,
                   height   = 15);
            dev.off()},
         .progress = TRUE)

## Rendering the animation of the best model
best <-
    best %>%
    mutate(Nyear = map_int(Data, ~ nlevels(.[["Year_fac"]])))
           
restart_pb(nrow(best))
best %$%
    pmap(list(Data, FitFunc, DataID, Species, Population, Nyear),
         ~ {pb$tick();
            anim_save(animation = animate(animate_fitfunc(..1, ..2),
                                          renderer = ffmpeg_renderer(),
                                          nframes = 20 * ..6),
                      filename = str_glue("../../Figures/Fitted_STAN_best/{..3}_{..4}_{..5}.mp4"))},
         .progress = TRUE)

## Plotting the AR1 (most complete) model
ar1 <-
    mods %>%
    filter(Model_Year == "AR1") %T>%
    {restart_pb(nrow(.))} %>%
    mutate(Curve_df = map2(Data, FitFunc, ~ construct_curve_df(.x, .y)),
           Plot = map2(Curve_df, FitFunc, ~ plot_fitfunc(.x, .y)))

ar1 %>%
    with(
        pmap(list(Plot, DataID, Species, Population),
             ~ {ggsave(plot = ..1,
                       filename = str_glue("../../Figures/Fitted_STAN_AR1/{..2}_{..3}_{..4}.pdf"),
                       device   = cairo_pdf(),
                       width    = 7,
                       height   = 7);
                dev.off()},
             .progress = TRUE)
    )

## Using the new IDs
ar1 <-
    ar1 %>%
    select(-ID) %>%
    left_join(ids) %>%
    select(ID, everything()) %>%
    arrange(ID)

## Removing the datasets for which NoSel is the best
ar1 <- filter(ar1, !(ID %in% c("Mel", "Oam")))

## Unnest the curve dfs and standardise the "year"
ar1_unst <-
    ar1 %>%
    unnest(Curve_df) %>%
    group_by(ID) %>%
    mutate(Year_fac   = factor(Year_fac, levels = sort(as.numeric(levels(Year_fac)))),
           Year_scale = scale_minmax(Year_index)) %>%
    ungroup()

## Generating the overall plot
p <-
    ggplot(ar1_unst) +
    geom_line(mapping   = aes(x      = Pheno,
                              y      = Y,
                              colour = Year_scale,
                              group  = Year_fac),
              size      = 1,
              alpha     = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = "Phenological trait", y = "Fitness trait") +
    facet_wrap(~ ID, nrow = 7, ncol = 6, scales = "free") +
    scale_size_area(name = "Counts") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_colour_gradient(low = "#0055ff", high = "#ffd500", guide = "none") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 30),
          strip.text   = element_text(family = "Linux Biolinum O", size = 16))
   
cairo_pdf("../../Figures/Fitted_STAN_AR1/All.pdf", height = 20, width = 15)
plot(p)
dev.off()
