##install.packages("picante", dependencies = TRUE)
library(picante) ## provides both library(vegan) and library(ape) ## ade4?
##install.packages("lme4")
library(lme4) ## for (generalized) linear mixed effects models
##install.packages("lmerTest")
library(lmerTest) ## p-values in summary for (generalized) linear mixed effects models
##install.packages("foreach")
library(foreach)
##install.packages("tidyverse")
library(tidyverse)
##install.packages("corrplot")
library(corrplot)
##install.packages("cowplot")
library(cowplot) ## arrange ggplots in a grid
##install.packages("FactoMineR", dependencies = TRUE)
library(FactoMineR) ## additional ordination methods
##install.packages("factoextra", dependencies = TRUE)
library(factoextra) ## additional ordination visualisation
##install.packages("xtable")
library(xtable) ## exports tables in latex code
##install.packages("ggsci")
library(ggsci) ## scientific color scales
##install.packages("MuMIn") # not available for R 3.4.4
library(MuMIn)

## make sure you have run the shell code to produce the file `twothree` with all completed replicates
## tworuns = read_table("twothree", col_names=c("nruns", "seed"))
## allfiles = Sys.glob("2019-0*/stats*tsv")
## mytworesults = allfiles[grep(paste(tworuns$seed, collapse="|"), allfiles)]
## rm(allfiles, tworuns)

## Incomplete runs are filtered out later
dispmode = "local" ## 'local' or 'global'
mytworesults = Sys.glob(paste0("data/2020*", dispmode, "*/*tsv"))

 rawresults = tibble()
     for (filepath in mytworesults) {
         ##if(length(grep("_sg_", filepath)) == 1 | length(grep("_sgv_", filepath)) == 1) rawresults = bind_rows(rawresults, read_tsv(filepath))
         rawresults = bind_rows(rawresults, read_tsv(filepath))
     }

## XXX this was originally `filter(x==1, y==1, time==1000)`, but that didn't work for global
 repstable = rawresults %>% filter(time==1000) %>% select(replicate, conf) %>% group_by(conf) %>% unique %>% table
 doublereps = which(rowSums(repstable) == 2) %>% names %>% as.numeric
 filteredresults = rawresults %>% filter(replicate %in% doublereps)
 
 tworesults = filteredresults %>% #rawresults %>%
     select(-ngenesstd, -area, -contains("compat"), -contains("reptol"), -contains("adaption")) %>%
     mutate(linkage_degree=ngenesmean/nlnkgunitsmean,
            scenario=ifelse(grepl("constant", conf), "static", "variable")) %>%
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

tworesults = tworesults %>% mutate(species = paste0(scenario, ".", lineage)) %>% #select(-contains("seedsize")) %>%  # seedsize is similar between scenarios + correlated/somewhat redundant with repsize
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

## species stats

worldend = tworesults[which(tworesults$time==1000),]
speccon = unique(worldend[which(worldend$scenario=="static"),]$lineage)
specvar = unique(worldend[which(worldend$scenario=="variable"),]$lineage)
specconex = setdiff(speccon, specvar)
specvarex = setdiff(specvar, speccon)
specboth = intersect(specvar, speccon)

print(paste("Species surviving only in", dispmode, "static runs:", length(specconex)))
print(paste("Species surviving only in", dispmode, "variable runs:", length(specvarex)))
print(paste("Species surviving in both", dispmode, "runs:", length(specboth)))
print(paste("(Considering n =", length(unique(worldend$replicate)), "replicates.)"))

## static plot over all replicates:
all = tworesults %>% rename(Environment = scenario) %>% filter(time == 500) %>% group_by(x, y, Environment) %>%
    summarize(`Richness (n.spp.)` = length(unique(lineage)) / length(unique(replicate)), `Adult biomass (g)` = mean(adult_body_size), `Number of genes` = mean(number_of_genes)) %>%
    mutate(Environment = ifelse(Environment == "static", "Static environment", "Variable environment"))
allrich = all %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Richness (n.spp.)`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic() +
    xlab("x-coordinate (temperature)") + ylab("y-coordinate (precipitation)")
allmass = all %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Adult biomass (g)`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic() +
    xlab("x-coordinate (temperature)") + ylab("y-coordinate (precipitation)")
allgenes = all %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Number of genes`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic() +
    xlab("x-coordinate (temperature)") + ylab("y-coordinate (precipitation)")
spatgridall = plot_grid(allrich, allmass, allgenes, labels="auto", ncol=1, align="vh")
ggsave(paste0("mapplots_all_", dispmode, ".pdf"), spatgridall, width=5, height=6)

## static plot for replicate 18:
## XXX Why replicate 18?
## XXX replaced with 22 for global
n = 18
d = tworesults %>% rename(Environment = scenario) %>% filter(replicate == n, time == 500) %>% group_by(x, y, Environment) %>%
    summarize(`Richness (n.spp.)` = length(unique(lineage)), `Adult biomass (g)` = mean(adult_body_size), `Number of genes` = mean(number_of_genes)) %>%
    mutate(Environment = ifelse(Environment == "static", "Static environment", "Variable environment"))
rich = d %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Richness (n.spp.)`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic() +
    xlab(expression("x-coordinate (" %~~% "temperature)")) + ylab("y-coordinate (precipitation)")
mass = d %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Adult biomass (g)`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic() +
    xlab("x-coordinate (temperature)") + ylab("y-coordinate (precipitation)")
genes = d %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Number of genes`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic() +
    xlab("x-coordinate (temperature)") + ylab("y-coordinate (precipitation)")
spatgrid = plot_grid(rich, mass, genes, labels="auto", ncol=1, align="vh")
ggsave(paste0("mapplots_s", n, "_", dispmode, ".pdf"), spatgrid, width=5, height=6)

## DV works
spatgridboth = plot_grid(rich, allrich, mass, allmass, genes, allgenes, labels="auto", ncol=2, align="vh")
ggsave(paste0("mapplots_both_", dispmode, ".pdf"), spatgridboth, width=10, height=6)

## animated: ##FIXME DV doesn't work yet
da = tworesults %>% rename(Environment = scenario) %>% filter(replicate == 18, time > 0) %>%
    mutate(Environment = ifelse(Environment == "static", "Static environment", "Variable environment")) %>%
    group_by(x, y, Environment, time) %>%
    summarize(`Richness / n.spp.` = length(unique(lineage)), `Adult biomass (g)` = mean(repsize), `Number of genes` = mean(ngenes))
richs = da %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Richness / n.spp.`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic() + transition_time(time) + labs(title = "Year: {frame_time}")
anim_save(paste0("mapplots_rich_s", 18, "_", dispmode, ".gif"), richs, nframes = length(unique(da$time)), fps = 2)
masss = da %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Adult biomass (g)`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic() + transition_time(time) + labs(title = "Year: {frame_time}")
anim_save(paste0("mapplots_mass_s", 18, "_", dispmode, ".gif"), masss, nframes = length(unique(da$time)), fps = 2)
geness = da %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Number of genes`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic() + transition_time(time) + labs(title = "Year: {frame_time}")
anim_save(paste0("mapplots_genes_s", 18, "_", dispmode, ".gif"), geness, nframes = length(unique(da$time)), fps = 2)
richa = da %>% ggplot(aes(x, y)) + geom_tile(aes(fill = Richness)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic()
massa = da %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Adult biomass`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic()
genesa = da %>% ggplot(aes(x, y)) + geom_tile(aes(fill = `Number of genes`)) + coord_fixed() +
    facet_grid(. ~ Environment) + scale_fill_viridis_c(option="magma") + theme_classic()
spatgrida = plot_grid(richa, massa, genesa, labels="auto", ncol=1, align="vh") + transition_time(time)
anim_save(paste0("mapplots_s", 18, "_", dispmode, ".gif"), spatgrida, width=9, height=7)

## diversity over time

lclrich = tworesults %>% filter(time>=50) %>% group_by(time, x, y, scenario, replicate) %>% summarize(alpha_diversity = length(unique(lineage))) %>%
      ungroup %>% group_by(time, scenario, replicate) %>% summarize_at(vars(alpha_diversity), mean) %>%
    ggplot(aes(time, alpha_diversity, group=scenario)) + stat_summary(aes(color=scenario), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d(name = "Environment") + theme_bw() + ylab(expression(paste(alpha, "-diversity", sep = ""))) + xlab("Year")
ggsave(paste0("localrichness_over_time_", dispmode, ".pdf"), lclrich, width=6, height=4)

mybeta = tibble(time=numeric(), scenario=character(), replicate=numeric(), betadiv=numeric(), zetadiv=numeric(), zetasd=numeric())
for (ts in unique(tworesults$time)) {
    for (cf in unique(tworesults$scenario)) {
        for (r in unique(tworesults$replicate)) {
            mycom = tworesults %>% filter(time==ts, scenario==cf, replicate==r) %>% mutate(abundance=juveniles+adults) %>% group_by(x, y) %>%
                select(abundance, lineage) %>% spread(key=lineage, value=abundance, fill=0) %>% ungroup() %>% select(-x, -y)
            if(nrow(mycom > 0)) {
                betadiv = betadiver(mycom, "w")
            } else {
                betadiv = NA
            }
            mybeta = bind_rows(mybeta, list(time=ts, scenario=cf, replicate=r, beta_diversity=mean(betadiv)))
        }
    }
}

beta =  mybeta %>% filter(time>=50) %>% ggplot(aes(time, beta_diversity, group=scenario)) + stat_summary(aes(color=scenario), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d() + theme_bw() + ylab(expression(paste(beta, "-diversity", sep = ""))) + xlab("Year")
ggsave(paste0("betadiv_over_time_", dispmode, ".pdf"), beta, width=6, height=4)

ttlrich = tworesults %>% filter(time>=50) %>% select(-x, -y) %>% group_by(time, scenario, replicate) %>% summarize(gamma_diversity = length(unique(lineage))) %>%
    ggplot(aes(time, gamma_diversity, group=scenario)) + stat_summary(aes(color=scenario), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d() + theme_bw() + ylab(expression(paste(gamma, "-diversity", sep = ""))) + xlab("Year")
ggsave(paste0("totalrichness_over_time_", dispmode, ".pdf"), ttlrich, width=6, height=4)

## population development

juvs = tworesults %>% filter(time>=50) %>% select(-x, -y) %>% group_by(time, scenario, replicate) %>%
    ggplot(aes(time, juveniles, group=scenario)) + stat_summary(aes(color=scenario), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d(name="Environment") + theme_bw() + ylab("Number of juveniles") + xlab("Year")
ggsave(paste0("juveniles_over_time_", dispmode, ".pdf"), juvs, width=6, height=4)

adlts = tworesults %>% filter(time>=50) %>% select(-x, -y) %>% group_by(time, scenario, replicate) %>%
    ggplot(aes(time, adults, group=scenario)) + stat_summary(aes(color=scenario), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d() + theme_bw() + ylab("Number of adults") + xlab("Year")
ggsave(paste0("adults_over_time_", dispmode, ".pdf"), adlts, width=6, height=4)

## range filling

myenv = tworesults %>% group_by(time, scenario, replicate) %>% select(temp, prec) %>% unique() %>% ungroup()
myspecs = tworesults %>% group_by(time, scenario, replicate, lineage) %>% select(ends_with("range")) %>%
    summarize(minprecrange=min(minprecrange), maxprecrange=max(maxprecrange),
              mintemprange=min(mintemprange), maxtemprange=max(maxtemprange)) %>% mutate(rangefilling=0) %>% ungroup()
myspecs = myspecs %>% inner_join(myenv) %>% mutate(habitable = temp>=mintemprange & temp<=maxtemprange & prec>=minprecrange & prec<=maxprecrange) %>%
    group_by(time, scenario, replicate, lineage) %>% select(habitable) %>% summarise(rangefilling=sum(habitable)/length(habitable)) %>% ungroup()
range =  myspecs %>% filter(time>=50) %>% mutate(replicate=as.factor(replicate), scenario=as.factor(scenario)) %>% group_by(time, scenario, replicate) %>%
    ggplot(aes(time, rangefilling, group=scenario)) + stat_summary(aes(color=scenario), fun.y = mean, geom="line", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d() + theme_bw() + ylab("Range-filling") + xlab("Year")
ggsave(paste0("rangefilling_over_time_", dispmode, ".pdf"), range, width=6, height=4)

## complete fig. 2

ecogrid = plot_grid(lclrich + theme(legend.position=c(.6, .75)),
          beta + theme(legend.position="none"),
          ttlrich + theme(legend.position="none"),
          juvs + theme(legend.position="none"),
          adlts + theme(legend.position="none"),
          range + theme(legend.position="none"), labels="auto", ncol=3, align="vh")
pattsleg = plot_grid(ecogrid, ncol=1, rel_heights=c(1,.1)) # get_legend(juvs), 
ggsave(paste0("ecopatts_", dispmode, ".pdf"), pattsleg, width=7, height=5)

##XXX Not sure what this next section does, or whether we still need it
tworesults %>% filter(time>=50) %>% select(-x, -y) %>% group_by(time, scenario, replicate) %>%
    ggplot(aes(time, mean_genetic_variation, group=scenario)) + stat_summary(aes(color=scenario), fun.y = mean, geom="smooth", size=1) +
    stat_summary(fun.data=mean_cl_boot, geom="ribbon", alpha=0.1) + scale_color_viridis_d() + theme_bw() + ylab("Mean genetic variation") + xlab("Year")

lineagevec = tworesults %>% filter(time==500) %>% select(scenario, lineage)
 sharedspecies = intersect(lineagevec[lineagevec$scenario=="static",]$lineage, lineagevec[lineagevec$scenario=="variable",]$lineage)
 staticspecies = setdiff(lineagevec[lineagevec$scenario=="static",]$lineage, lineagevec[lineagevec$scenario=="variable",]$lineage)
 variablespecies = setdiff(lineagevec[lineagevec$scenario=="variable",]$lineage, lineagevec[lineagevec$scenario=="static",]$lineage)
 uniquespecies = c(staticspecies, variablespecies)
 c(length(lineagevec[lineagevec$scenario=="static",]$lineage), length(lineagevec[lineagevec$scenario=="variable",]$lineage))
 summary(is.element(unique(lineagevec[lineagevec$scenario=="static",]$lineage), unique(lineagevec[lineagevec$scenario=="variable",]$lineage)))
 summary(is.element(unique(lineagevec[lineagevec$scenario=="variable",]$lineage), unique(lineagevec[lineagevec$scenario=="static",]$lineage)))

 pdf(paste0("venn_sp_500_", dispmode, ".pdf"), width=4, height=3)
 draw.pairwise.venn(length(unique(lineagevec[lineagevec$scenario=="static",]$lineage)),
                    length(unique(lineagevec[lineagevec$scenario=="variable",]$lineage)),
                    sum(is.element(unique(lineagevec[lineagevec$scenario=="static",]$lineage), unique(lineagevec[lineagevec$scenario=="variable",]$lineage))),
                    c("static", "variable"), fill=viridisLite::viridis(2),
                    cat.dist=0.05, cat.pos=c(-45,45), margin=0.04, fontfamily="sans",
                    cat.fontfamily="sans")
 dev.off()

 summary(lineagevec$lineage %in% uniquespecies)

## PCA
mainendtraits = tworesults %>% filter(time == 500) %>%
    rename(Environment = scenario) %>%
    dplyr::select(Environment, mean_dispersal_distance, number_of_genes, precipitation_tolerance,
           adult_body_size, temperature_tolerance, linkage_degree, mean_genetic_variation,
           long_distance_dispersal, seed_size) %>%
    mutate_at(vars(-Environment), function(x) log(x + 1)) %>%
    rename(`Mean dispersal distance` = mean_dispersal_distance, `Number of genes` = number_of_genes,
           `Precipitation tolerance` = precipitation_tolerance, `Adult biomass (g)` = adult_body_size,
           `Temperature tolerance` = temperature_tolerance, `Genetic linkage` = linkage_degree,
           `Mean genetic variation` = mean_genetic_variation, `Long distance dispersal` = long_distance_dispersal,
           `Seed biomass (g)` = seed_size)

endpca = prcomp(mainendtraits[,-1], scale=T)
print(xtable(endpca, digits = 3, floating = FALSE, booktabs = TRUE, include.rownames=FALSE))


endpcaviz = fviz_pca_biplot(endpca, col.var=factor(c("ecological", "genetic", "ecological", "ecological", "ecological", "genetic", "genetic", "ecological", "ecological")),
                geom.ind="point", fill.ind=mainendtraits$Environment, pointsize=1, pointshape=21, addEllipses = TRUE) + #, ellipse.alpha=0.1, ellipse.type = "convex") +
    theme_bw() + scale_fill_viridis_d("Environment") + scale_color_brewer(palette="Set2", name="Trait")
ggsave(paste0("pca_t500_maintraits_", dispmode, ".pdf"), endpcaviz, width=4.5, height=4)
endpcascree = fviz_eig(endpca) + theme_bw() # => all traits similarly important for characterisation of trait space
pca_grid = plot_grid(endpcaviz, endpcascree, ncol=2, rel_widths=c(1, 0.5), labels="auto")
ggsave(paste0("pca_t500_maintraits_scree_", dispmode, ".pdf"), pca_grid, width=7, height=4)

## prepare data:
  myendresults = tworesults %>% filter((time == 0 & scenario == "static") | time == 500) %>% #mutate(shared=as.factor(ifelse(lineage %in% sharedspecies, "shared", "unique"))) %>%
      mutate_at(vars(contains("tolerance"), contains("size"), contains("gene"), contains("dispersal")), function(x) log(x + 1)) %>%
      mutate(scenario = ifelse(time == 0, "initial", scenario)) %>%
      mutate(scenario=as.factor(scenario)) %>%
      filter(scenario == "static" | scenario == "variable") %>%
      rename(Environment=scenario, log_adult_body_size=adult_body_size) %>%
      na.omit()

  ##DV works
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
ggsave(paste0("differences_traits_environments_replicate_means_", dispmode, ".pdf"),
       width = 5, height = 5)
  ggsave(paste0("diffs_means_", dispmode, ".pdf"), width = 5, height = 5)

  ## Trait variances (/CV):
  subtraitnames = myendresults %>% dplyr::select(contains("CV_median")) %>% names() 

  endsubtraits_lme = foreach(trait=subtraitnames) %do% {
      lmer(get(trait) ~ Environment + (1|replicate), data = myendresults)
  }
  names(endsubtraits_lme) = subtraitnames

  endsubtraits_lme_summary = lapply(endsubtraits_lme, summary)

  endsubtraits_lme_table = bind_cols(names = names(endsubtraits_lme_summary), as.tibble(t(sapply(endsubtraits_lme_summary, function(x) unlist(as.tibble(x$coefficients)[2,])))))
  endsubtraits_lme_table[,6] <= 0.05
  endsubtraits_lme_table$level =  factor(ifelse(grepl("genetic", endsubtraits_lme_table$names), "genetic", "phenotypic"),
                                         levels = c("Community means", "phenotypic", "genetic"))
  endsubtraits_lme_table$names = factor(rep(c("Mean dispersal distance", "Long distance dispersal",
     "Precipitation tolerance", "Seed biomass (g)", "Adult biomass (g)",
     "Temperature tolerance"), each = 2),
             levels = rev(c("Mean dispersal distance", "Long distance dispersal",
     "Precipitation tolerance", "Seed biomass (g)", "Adult biomass (g)", "Temperature tolerance", 
     "Number of genes", "Genetic linkage", "Mean genetic variation")))
  print(xtable(endsubtraits_lme_table[,c(1,7,2:6)], digits = c(0, 0, 0, 3, 3, 0, 3, 3)), floating = FALSE, booktabs = TRUE, include.rownames=FALSE)

##DV I have no idea what this is supposed to show...
endsubtraits_lme_table %>%
      ggplot(aes(names, Estimate, fill = ifelse(Estimate < 0, "1", "-1"))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1) +
      geom_bar(stat = "identity", width = 0.5, position = "dodge") +
      geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`), position = position_dodge(.5), width = 0) +
      xlab("") + ylab("Difference in variances between environments") +
      scale_y_continuous(limits = c(min(endsubtraits_lme_table[,"Estimate"] - endsubtraits_lme_table[,"Std. Error"]),
                                    max(endsubtraits_lme_table[,"Estimate"] + endsubtraits_lme_table[,"Std. Error"]) + 0.005)) +
      coord_flip() + scale_fill_npg(guide = FALSE) + facet_grid(.~level)
ggsave(paste0("differences_traits_environments_replicate_variances_", dispmode, ".pdf"),
       width = 10, height = 5)
  ggsave(paste0("diffs_variances_", dispmode, ".pdf"), width = 10, height = 5)

  lme_table$level = factor("Community means", levels = c("Community means", "Intraspecific variation", "Genetic variation"))

  combdiffs = bind_rows(lme_table, endsubtraits_lme_table) %>%
      ggplot(aes(names, Estimate, fill = ifelse(Estimate < 0, "-1", "1"))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1) +
      geom_bar(stat = "identity", width = 0.5, position = "dodge") +
      geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`), position = position_dodge(.5), width = 0) +
      xlab("") + ylab("Differences of variable compared to static environments") +
      coord_flip() + scale_fill_npg(guide = FALSE) + facet_grid(.~level, scales = "free")
  ggsave(paste0("all_diffs_variances_", dispmode, ".pdf"), combdiffs, width = 7, height = 3)

## PCA standard deviation plot (appendix fig. A1)
##FIXME not quite right yet :-/

reps = tworesults$replicate %>% unique() %>% sample()
tsteps = seq(800,1000,50)
anavar = tibble()
for (nreps in seq(10, length(table(tworesults$replicate)), 10)) {
    for (ntsteps in 1:length(tsteps)) {
        temp.res = tworesults %>% filter(time %in% tsteps[1:ntsteps], replicate %in% reps[1:nreps]) %>%
            filter(scenario=="variable")
        mypca = temp.res %>% select(-ends_with("sdmin")) %>%
            ##XXX deleted -maxage and -maxsize below:
            select(-(x:prec), -adults, -juveniles, -time, -replicate) %>%
            select_if(is.numeric) %>% select_if(function(x){!any(is.na(x))}) %>%
            select_if(function(x){var(x)!=0}) %>% prcomp(scale=T)
        anavar = bind_rows(anavar, c(Number_of_replicates=nreps, Number_of_timesteps = ntsteps, PC=mypca$sdev))
    }
}

pcasds = anavar %>% gather(contains("PC"), key=component, value=Standard_deviation, factor_key=T) %>%
    mutate(Number_of_replicates=as.factor(Number_of_replicates)) %>%
    ggplot(aes(Number_of_timesteps, Standard_deviation)) + geom_line(aes(color=Number_of_replicates)) +
    facet_wrap(.~component, scales="free_y") + scale_color_viridis_d() + theme_bw()
ggsave(paste0("pcasd_t800_timesteps_replicates_", dispmode, ".pdf"), pcasds, width=12, height=8)
