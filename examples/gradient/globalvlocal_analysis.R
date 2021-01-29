#!/usr/bin/Rscript
library(picante) ## provides both library(vegan) and library(ape) ## ade4?
library(lme4) ## for (generalized) linear mixed effects models
library(lmerTest) ## p-values in summary for (generalized) linear mixed effects models
library(foreach)
library(tidyverse)
library(corrplot)
library(cowplot) ## arrange ggplots in a grid
library(FactoMineR) ## additional ordination methods
library(factoextra) ## additional ordination visualisation
library(xtable) ## exports tables in latex code
library(ggsci) ## scientific color scales
library(MuMIn)

### LOAD AND PREPARE DATA

dday = 750 ## time step to analyse
dispmode = "both" ## 'local' or 'global'
mytworesults = Sys.glob(paste0("data/*/*tsv"))

rawresults = tibble()
for (filepath in mytworesults) {
    rawresults = bind_rows(rawresults, read_tsv(filepath))
}

## only use replicates where all scenarios are finished
repstable = rawresults %>% filter(time==dday) %>% select(replicate, conf) %>% group_by(conf) %>% unique %>% table
doublereps = which(rowSums(repstable) == 4) %>% names %>% as.numeric
filteredresults = rawresults %>% filter(replicate %in% doublereps) %>% filter(time<=dday)

tworesults = filteredresults %>%
    select(-ngenesstd, -area, -contains("compat"), -contains("reptol"), -contains("adaption")) %>%
    mutate(linkage_degree=ngenesmean/nlnkgunitsmean,
           scenario=ifelse(grepl("constant", conf), "static", "variable"),
           disptype=ifelse(grepl("global", conf), "b-global", "a-local"),
           config=sub("constant", "static", sub("_", " ", sub(".conf", "", sub("examples/gradient/", "", conf))))) %>%
    select(-contains("lnkgunits"), -conf) %>%
    ##XXX This doesn't calculate the full range, but as the model doesn't output min/max values
    ## anymore, it's the best we can do.
    mutate(mintemprange=(tempoptmean-(2*tempoptstd))-(temptolmean+(2*temptolstd)),
           maxtemprange=(tempoptmean+(2*tempoptstd))+(temptolmean+(2*temptolstd)),
           minprecrange=(precoptmean-(2*precoptstd))-(prectolmean+(2*prectolstd)),
           maxprecrange=(precoptmean+(2*precoptstd))+(prectolmean+(2*prectolstd))) %>%
    select(-ends_with("min"), -ends_with("max"), -ends_with("sdstd")) %>% na.omit()
names(tworesults) = names(tworesults) %>% gsub("std", "_pop._var.", .) %>% gsub("sdmean", "_gen._var.", .)  %>%
    ## This next line is a nasty hack arising from a switch in the data from `med` values to `mean`
    gsub("dispmean", "dispme-an", .) %>% gsub("mean", "", .) %>% gsub("dispme-an", "dispmean", .)

tworesults = tworesults %>% mutate(species = paste0(scenario, ".", lineage)) %>%
    rename(mean_dispersal_distance = dispmean, 
           long_distance_dispersal = dispshape, 
           precipitation_optimum = precopt, 
           precipitation_tolerance = prectol, 
           seed_size = seedsize,
           adult_body_size = repsize, 
           temperature_optimum = tempopt, 
           temperature_tolerance = temptol, 
           number_of_genes = ngenes) %>%
    mutate(mean_dispersal_distance_intra_CV_median = dispmean_pop._var. / mean_dispersal_distance,
           mean_dispersal_distance_genetic_CV_median = dispmean_gen._var. / mean_dispersal_distance,
           long_distance_dispersal_intra_CV_median = dispshape_pop._var. / long_distance_dispersal,
           long_distance_dispersal_genetic_CV_median = dispshape_gen._var. / long_distance_dispersal,
           precipitation_tolerance_intra_CV_median = prectol_pop._var. / precipitation_tolerance,
           precipitation_tolerance_genetic_CV_median = prectol_gen._var. / precipitation_tolerance,
           seed_size_intra_CV_median = seedsize_pop._var. / seed_size,
           seed_size_genetic_CV_median = seedsize_gen._var. / seed_size,
           adult_body_size_intra_CV_median = repsize_pop._var. / adult_body_size,
           adult_body_size_genetic_CV_median = repsize_gen._var. / adult_body_size,
           temperature_tolerance_intra_CV_median = temptol_pop._var. / temperature_tolerance,
           temperature_tolerance_genetic_CV_median = temptol_gen._var. / temperature_tolerance)
genetic_variation = tworesults %>% select(ends_with("genetic_CV_median")) %>% as_tibble() %>% rowMeans
tworesults = bind_cols(tworesults, mean_genetic_variation=genetic_variation)

mybeta = tibble(time=numeric(), config=character(), replicate=numeric(), betadiv=numeric(), zetadiv=numeric(), zetasd=numeric())
for (ts in unique(tworesults$time)) {
    for (cf in unique(tworesults$config)) {
        for (r in unique(tworesults$replicate)) {
            mycom = tworesults %>% filter(time==ts, config==cf, replicate==r) %>% mutate(abundance=juveniles+adults) %>% group_by(x, y) %>%
                select(abundance, lineage) %>% spread(key=lineage, value=abundance, fill=0) %>% ungroup() %>% select(-x, -y)
            if(nrow(mycom > 0)) {
                betadiv = betadiver(mycom, "w")
            } else {
                betadiv = NA
            }
            mybeta = bind_rows(mybeta, list(time=ts, config=cf, replicate=r, beta_diversity=mean(betadiv)))
        }
    }
}


### POPULATIONS AND DIVERSITY OVER TIME

lclrich = tworesults %>% filter(time>=50) %>% group_by(time, x, y, config, replicate) %>%
    summarize(alpha_diversity = length(unique(lineage))) %>%
    ungroup %>% group_by(time, config, replicate) %>%
    summarize_at(vars(alpha_diversity), mean) %>%
    ggplot(aes(time, alpha_diversity, group=config)) +
    stat_summary(aes(color=config), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) +
    scale_color_viridis_d(name = "Environment") +
    ylab(expression(paste(alpha, "-diversity", sep = ""))) +
    xlab("Year") + theme_bw()
ggsave(paste0("localrichness_over_time_", dispmode, ".pdf"), lclrich, width=6, height=4)

beta =  mybeta %>% filter(time>=50) %>% ggplot(aes(time, beta_diversity, group=config)) + stat_summary(aes(color=config), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d() + theme_bw() + ylab(expression(paste(beta, "-diversity", sep = ""))) + xlab("Year")
ggsave(paste0("betadiv_over_time_", dispmode, ".pdf"), beta, width=6, height=4)

ttlrich = tworesults %>% filter(time>=50) %>% select(-x, -y) %>% group_by(time, config, replicate) %>% summarize(gamma_diversity = length(unique(lineage))) %>%
    ggplot(aes(time, gamma_diversity, group=config)) + stat_summary(aes(color=config), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d() + theme_bw() + ylab(expression(paste(gamma, "-diversity", sep = ""))) + xlab("Year")
ggsave(paste0("totalrichness_over_time_", dispmode, ".pdf"), ttlrich, width=6, height=4)

juvs = tworesults %>% filter(time>=50) %>% select(-x, -y) %>% group_by(time, config, replicate) %>%
    ggplot(aes(time, juveniles, group=config)) + stat_summary(aes(color=config), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d(name="Environment") + theme_bw() + ylab("Number of juveniles") + xlab("Year")
ggsave(paste0("juveniles_over_time_", dispmode, ".pdf"), juvs, width=6, height=4)

adlts = tworesults %>% filter(time>=50) %>% select(-x, -y) %>% group_by(time, config, replicate) %>%
    ggplot(aes(time, adults, group=config)) + stat_summary(aes(color=config), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d() + theme_bw() + ylab("Number of adults") + xlab("Year")
ggsave(paste0("adults_over_time_", dispmode, ".pdf"), adlts, width=6, height=4)

myenv = tworesults %>% group_by(time, config, replicate) %>% select(temp, prec) %>% unique() %>% ungroup()
myspecs = tworesults %>% group_by(time, config, replicate, lineage) %>% select(ends_with("range")) %>%
    summarize(minprecrange=min(minprecrange), maxprecrange=max(maxprecrange),
              mintemprange=min(mintemprange), maxtemprange=max(maxtemprange)) %>% mutate(rangefilling=0) %>% ungroup()
myspecs = myspecs %>% inner_join(myenv) %>% mutate(habitable = temp>=mintemprange & temp<=maxtemprange & prec>=minprecrange & prec<=maxprecrange) %>%
    group_by(time, config, replicate, lineage) %>% select(habitable) %>% summarise(rangefilling=sum(habitable)/length(habitable)) %>% ungroup()
range =  myspecs %>% filter(time>=50) %>% mutate(replicate=as.factor(replicate), config=as.factor(config)) %>% group_by(time, config, replicate) %>%
    ggplot(aes(time, rangefilling, group=config)) + stat_summary(aes(color=config), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d() + theme_bw() + ylab("Range-filling") + xlab("Year")
ggsave(paste0("rangefilling_over_time_", dispmode, ".pdf"), range, width=6, height=4)

ecogrid = plot_grid(lclrich + theme(legend.position="none"),
          beta + theme(legend.position="none"),
          ttlrich + theme(legend.position=c(.6, .75)),
          juvs + theme(legend.position="none"),
          adlts + theme(legend.position="none"),
          range + theme(legend.position="none"), labels="auto", ncol=3, align="vh")
pattsleg = plot_grid(ecogrid, ncol=1, rel_heights=c(1,.1)) # get_legend(juvs), 
ggsave(paste0("ecopatts_", dispmode, ".pdf"), pattsleg, width=7, height=5)


### TRAIT ANALYSIS (PCA AND TRAIT MEANS)

mainendtraits = tworesults %>% filter(time == dday) %>%
    rename(Scenario = config) %>%
    dplyr::select(Scenario, mean_dispersal_distance, number_of_genes, precipitation_tolerance,
           adult_body_size, temperature_tolerance, linkage_degree, mean_genetic_variation,
           long_distance_dispersal, seed_size) %>%
    mutate_at(vars(-Scenario), function(x) log(x + 1)) %>%
    rename(`Mean dispersal distance` = mean_dispersal_distance, `Number of genes` = number_of_genes,
           `Precipitation tolerance` = precipitation_tolerance, `Adult biomass (g)` = adult_body_size,
           `Temperature tolerance` = temperature_tolerance, `Genetic linkage` = linkage_degree,
           `Mean genetic variation` = mean_genetic_variation, `Long distance dispersal` = long_distance_dispersal,
           `Seed biomass (g)` = seed_size)

endpca = prcomp(mainendtraits[,-1], scale=T)
print(xtable(endpca, digits = 3, floating = FALSE, booktabs = TRUE, include.rownames=FALSE))

endpcaviz = fviz_pca_biplot(endpca, col.var=factor(c("ecological", "genetic", "ecological", "ecological", "ecological", "genetic", "genetic", "ecological", "ecological")),
                geom.ind="point", fill.ind=mainendtraits$Scenario, pointsize=1, pointshape=21, addEllipses = TRUE) + #, ellipse.alpha=0.1, ellipse.type = "convex") +
    theme_bw() + scale_fill_viridis_d("Scenario") + scale_color_brewer(palette="Set2", name="Trait")
ggsave(paste0("pca_t",dday,"_maintraits_", dispmode, ".pdf"), endpcaviz, width=4.5, height=4)

## prepare data for separate analysis (global vs local dispersal)
scen = "variable"
myendresults = tworesults %>% filter(time == dday) %>% filter(scenario==scen) %>%
    mutate_at(vars(contains("tolerance"), contains("size"), contains("gene"), contains("dispersal")), function(x) log(x + 1)) %>%
    mutate(disptype=as.factor(disptype)) %>%
    rename(Dispersal=disptype, log_adult_body_size=adult_body_size) %>%
    na.omit()

## Trait means:
traitnames = myendresults %>% dplyr::select(mean_dispersal_distance, long_distance_dispersal, number_of_genes, precipitation_tolerance, seed_size, 
                                            log_adult_body_size, temperature_tolerance, linkage_degree, mean_genetic_variation) %>%
    names() 

endtraits_lme = foreach(trait=traitnames) %do% {
    lmer(get(trait) ~ Dispersal + (1|replicate), data = myendresults)
}
names(endtraits_lme) = traitnames

endtraits_lme_summary = lapply(endtraits_lme, summary)

lme_table = bind_cols(names = names(endtraits_lme_summary), as_tibble(t(sapply(endtraits_lme_summary, function(x) unlist(as.tibble(x$coefficients)[2,])))))
lme_table$names = factor(c("Mean dispersal distance", "Long distance dispersal", "Number of genes",
                           "Precipitation tolerance", "Seed biomass (g)", "Adult biomass (g)",
                           "Temperature tolerance", "Genetic linkage", "Mean genetic variation"),
                         levels = rev(c("Mean dispersal distance", "Long distance dispersal",
                                        "Precipitation tolerance", "Seed biomass (g)", "Adult biomass (g)", "Temperature tolerance", 
                                        "Number of genes", "Genetic linkage", "Mean genetic variation")))
print(xtable(lme_table, digits = c(0, 0, 3, 3, 0, 3, 3)), floating = FALSE, booktabs = TRUE, include.rownames=FALSE)

lme_table %>%
    ggplot(aes(names, Estimate, fill = ifelse(Estimate > 0, "1", "-1"))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1) +
    geom_bar(stat = "identity", width = 0.5, position = "dodge") +
    geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`), position = position_dodge(.5), width = 0) +
    scale_y_continuous(limits = c(min(lme_table[,"Estimate"] - lme_table[,"Std. Error"]),
                                  max(lme_table[,"Estimate"] + lme_table[,"Std. Error"]) + 0.01)) +
    xlab("") + ylab("Mean trait values with global relative to local pollen movement") + coord_flip() +
    scale_fill_npg(guide = FALSE)
ggsave(paste0("diffs_means_", dispmode, "_", scen, "_separate.pdf"), width = 5, height = 5)


subtraitnames = myendresults %>% dplyr::select(contains("CV_median")) %>% names() 
endsubtraits_lme = foreach(trait=subtraitnames) %do% {
    lmer(get(trait) ~ Dispersal + (1|replicate), data = myendresults)
}
names(endsubtraits_lme) = subtraitnames
endsubtraits_lme_summary = lapply(endsubtraits_lme, summary)
endsubtraits_lme_table = bind_cols(names = names(endsubtraits_lme_summary),
                                   as.tibble(t(sapply(endsubtraits_lme_summary, function(x) unlist(as.tibble(x$coefficients)[2,])))))
endsubtraits_lme_table[,6] <= 0.05
endsubtraits_lme_table$level = factor(ifelse(grepl("genetic", endsubtraits_lme_table$names), "Genetic standing variation", "Phenotypic standing variation"),
                                       levels = c("Community means", "Phenotypic standing variation", "Genetic standing variation"))
endsubtraits_lme_table$names = factor(rep(c("Mean dispersal distance", "Long distance dispersal",
                                            "Precipitation tolerance", "Seed biomass (g)", "Adult biomass (g)",
                                            "Temperature tolerance"), each = 2),
                                      levels = rev(c("Mean dispersal distance", "Long distance dispersal",
                                                     "Precipitation tolerance", "Seed biomass (g)", "Adult biomass (g)", "Temperature tolerance", 
                                                     "Number of genes", "Mean genetic variation")))
print(xtable(endsubtraits_lme_table[,c(1,7,2:6)], digits = c(0, 0, 0, 3, 3, 0, 3, 3)), floating = FALSE, booktabs = TRUE, include.rownames=FALSE)
lme_table$level = factor("Community means", levels = c("Community means", "Phenotypic standing variation", "Genetic standing variation"))

combdiffs = bind_rows(lme_table, endsubtraits_lme_table) %>%
    ggplot(aes(names, Estimate, fill = ifelse(Estimate < 0, "-1", "1"))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1) +
    geom_bar(stat = "identity", width = 0.5, position = "dodge") +
    geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`), position = position_dodge(.5), width = 0) +
    xlab("") + ylab("Mean trait values with global relative to local pollen movement") +
    coord_flip() + scale_fill_npg(guide = FALSE) + facet_grid(.~level, scales = "free")
ggsave(paste0("diffs_means_", dispmode, "_", scen, "_all.pdf"), combdiffs, width = 7, height = 3)



## prepare data for combined analysis (global and local dispersal)
myendresults = tworesults %>% filter((time == 0 & scenario == "static") | time == dday) %>% 
    mutate_at(vars(contains("tolerance"), contains("size"), contains("gene"), contains("dispersal")), function(x) log(x + 1)) %>%
    mutate(scenario = ifelse(time == 0, "initial", scenario)) %>%
    mutate(scenario=as.factor(scenario)) %>%
    filter(scenario != "initial") %>%
    rename(Environment=scenario, log_adult_body_size=adult_body_size) %>%
    na.omit()

## Trait means:
traitnames = myendresults %>% dplyr::select(mean_dispersal_distance, long_distance_dispersal, number_of_genes, precipitation_tolerance, seed_size, 
                                            log_adult_body_size, temperature_tolerance, linkage_degree, mean_genetic_variation) %>%
    names() 

endtraits_lme = foreach(trait=traitnames) %do% {
    lmer(get(trait) ~ Environment + (1|replicate), data = myendresults)
}
names(endtraits_lme) = traitnames

endtraits_lme_summary = lapply(endtraits_lme, summary)

lme_table = bind_cols(names = names(endtraits_lme_summary), as_tibble(t(sapply(endtraits_lme_summary, function(x) unlist(as.tibble(x$coefficients)[2,])))))
lme_table$names = factor(c("Mean dispersal distance", "Long distance dispersal", "Number of genes",
                           "Precipitation tolerance", "Seed biomass (g)", "Adult biomass (g)",
                           "Temperature tolerance", "Genetic linkage", "Mean genetic variation"),
                         levels = rev(c("Mean dispersal distance", "Long distance dispersal",
                                        "Precipitation tolerance", "Seed biomass (g)", "Adult biomass (g)", "Temperature tolerance", 
                                        "Number of genes", "Genetic linkage", "Mean genetic variation")))
print(xtable(lme_table, digits = c(0, 0, 3, 3, 0, 3, 3)), floating = FALSE, booktabs = TRUE, include.rownames=FALSE)

lme_table %>%
    ggplot(aes(names, Estimate, fill = ifelse(Estimate > 0, "1", "-1"))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1) +
    geom_bar(stat = "identity", width = 0.5, position = "dodge") +
    geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`), position = position_dodge(.5), width = 0) +
    scale_y_continuous(limits = c(min(lme_table[,"Estimate"] - lme_table[,"Std. Error"]),
                                  max(lme_table[,"Estimate"] + lme_table[,"Std. Error"]) + 0.01)) +
    xlab("") + ylab("Difference in means between environments") + coord_flip() +
    scale_fill_npg(guide = FALSE)
ggsave(paste0("diffs_means_", dispmode, "_combined.pdf"), width = 5, height = 5)

