# ---
# title: "Lygus morphometrics"
# author: "Namyatova et al."
# ---

## Packages required
#install.packages("vegan")
#install.packages("ggplot2")
#install.packages("grid")
#install.packages("ggrepel")
#install.packages("cowplot")
#install.packages("car")

library(vegan)
library(ggplot2)
library(grid)
library(ggrepel)
library(cowplot)
library(scales)
library(car)

ar <- arrow(length = unit(0.2, "cm")) #arrow for loading plot from grid package
theme_set(theme_bw(base_size = 15))

#Centering function
center <- function(x){
  x - mean(x, na.rm = TRUE)} #NB! for double-centering it's enough to use this function to center only by column as centering by row is already incorporated into rda function we're using for PCA

## Data -----
lygus_meas <- read.table(file="Measurements_lygus.csv", sep=",", dec=".", header=TRUE, stringsAsFactors = TRUE)
lygus_meas$Species_Gender <- as.factor(paste(lygus_meas$Species, lygus_meas$Gender))
lygus_meas <- lygus_meas[, c(1:3, 21, 4:20)]
str(lygus_meas)
colSums(is.na(lygus_meas))

# Palette to tie the exact color to each species
palette_S = setNames(object = scales::hue_pal()(5), nm = unique(lygus_meas$Species))
palette_SG = setNames(object = brewer_pal('qual', 'Paired')(10), nm = unique(lygus_meas$Species_Gender))

## 1. PCA without double-centering all data-----
lygus_pca <- rda(lygus_meas[, -c(1:4)], scale = TRUE)
summary(lygus_pca)
eig <- eigenvals(lygus_pca)[1:17]
eig*100/sum(eig)

# Estimating the number of significant PCs using broken-stick model
screeplot(lygus_pca, type = "lines", bstick = TRUE)

### Loadings and loading plot ----
scores(lygus_pca, display = "species", choices = c(1, 2, 3), scaling = 0)
df_load <- as.data.frame(scores(lygus_pca, display = "species", choices = c(1, 2, 3), scaling = "species"))
p_load <- ggplot(df_load) +
  geom_text_repel(aes(x = PC1, y = PC2, label = rownames(df_load)), max.overlaps = 50,
                  size = 2, direction ='both', nudge_x = 0.2, seed = 256) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               color = "coral2", arrow = ar, alpha = 0.7) +
  coord_equal(xlim = c(-0.5, 2.5), ylim = c(-1.5, 1.5)) + 
  ggtitle('Loading plot, all,\nnot double-centered')
  

# Data for score plots
df_scores <- data.frame(lygus_meas[, 1:4],
                        scores(lygus_pca, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_scores)

## 1a. Score plot, grouped by sex ----
p_scores_G <- ggplot(df_scores, aes(x = PC1, y = PC2, color = Gender)) + 
  geom_point(size = 3.5, alpha = 0.6) + 
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle('Score plot, all,\nnot double-centered, by sex')+
  theme(legend.text = element_text(face = "italic"))

# Loading and score plots
both_plots_G <- plot_grid(p_load, p_scores_G, align = "h", rel_widths = c(0.45, 0.57))
#ggsave('no_dbc_all_gender_loadings_ordination.tiff', plot = both_plots_G, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

## 1b. Score plot, grouped by sex and species ----
p_scores_SG <- ggplot(df_scores, aes(x = PC1, y = PC2, color = Species_Gender)) + 
  geom_point(size = 3.5, alpha = 0.6) + 
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle('Score plot, all, not double-centered,\n by sex and species') +
  scale_color_manual(values=palette_SG) +
  theme(legend.text = element_text(face = "italic"))

# Loading and score plot
both_plots_SG <- plot_grid(p_load, p_scores_SG, align = "h", rel_widths = c(0.38, 0.62))
#ggsave('no_dbc_all_speciesgender_loadings_ordination.tiff', plot = both_plots_SG, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

### Analysis of variance of PC1 ----
df_SG <- data.frame(species_gender = lygus_meas$Species_Gender,
                 scores(lygus_pca, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_SG)
mod1SG <- lm(PC1 ~ species_gender, data = df_SG)
anova(mod1SG)

# Testing for conditions
mod_diag1SG <- fortify(mod1SG)
res_p1SG <- ggplot(data = mod_diag1SG, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val1SG <- mean(mod_diag1SG$.stdresid)
sd_val1SG <- sd(mod_diag1SG$.stdresid)
norm_p1SG <- ggplot(mod_diag1SG, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val1SG, slope = sd_val1SG)
plot_grid(res_p1SG, norm_p1SG, ncol = 2, rel_widths = c(0.55, 0.45))

# Point range graph for PC1
df_SG$species_gender <- reorder(df$species_gender, df_SG$PC1, FUN=mean)

pc1_plot_SG<-ggplot(df_SG, aes(x = species_gender, y = PC1, color = species_gender)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9)) + 
  ggtitle('Point range graph for PC1, all,\nnot double-centered, by sex and species') +
  scale_color_manual(values = palette_SG) +
  theme(legend.position = 'none', axis.text.x = element_text(face = "italic"))

# Tukey's test of PC1
TukeyHSD(aov(mod1SG))

### Analysis of variance of PC2 ----
mod2SG <- lm(PC2 ~ species_gender, data = df_SG)
anova(mod2SG)

# Testing for conditions
mod_diag2SG <- fortify(mod2SG)
res_p2SG <- ggplot(data = mod_diag2SG, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val2SG <- mean(mod_diag2SG$.stdresid)
sd_val2SG <- sd(mod_diag2SG$.stdresid)
norm_p2SG <- ggplot(mod_diag2SG, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val2SG, slope = sd_val2SG)
plot_grid(res_p2SG, norm_p2SG, ncol = 2, rel_widths = c(0.55, 0.45))

# Point range graph for PC2
df_SG$species_gender <- reorder(df_SG$species_gender, df$PC2, FUN=mean)

pc2_plot_SG<-ggplot(df_SG, aes(x = species_gender, y = PC2, color = species_gender)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9)) + 
  ggtitle('Point range graph for PC2, all,\nnot double-centered, by sex and species') +
  scale_color_manual(values = palette_SG) +
  theme(legend.text = element_text(face = "italic"), axis.text.x = element_text(face = "italic"))

# Tukey's test of PC2
TukeyHSD(aov(mod2SG))

both_plots_pc_SG <- plot_grid(pc1_plot_SG, pc2_plot_SG, rel_widths = c(0.40, 0.63))
#ggsave('no_dbc_all_speciesgender_PC.tiff', plot = both_plots_pc_SG, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

## 2. PCA after double-centering all data ----

dbc_lygus <- t(apply(lygus_meas[, 5:21], 1, center))
lygus_pca_dbc <- rda(dbc_lygus, scale = TRUE)
summary(lygus_pca_dbc)
eig <- eigenvals(lygus_pca_dbc)[1:16]
eig*100/sum(eig)

# Estimating the number of significant PCs using broken-stick model
screeplot(lygus_pca_dbc, type = "lines", bstick = TRUE)

### Loadings and loading plot ----
scores(lygus_pca_dbc, display = "species", choices = c(1, 2, 3), scaling = 0)
df_load_dbc <- as.data.frame(scores(lygus_pca_dbc, display = "species", choices = c(1, 2, 3), scaling = "species"))
p_load_dbc <- ggplot(df_load_dbc) +
  geom_text_repel(aes(x = PC1, y = PC2, label = rownames(df_load_dbc)),
                  size = 2, max.overlaps = 30, force_pull=0, seed = 230) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               color = "coral2", arrow = ar, alpha = 0.7) +
  coord_equal(xlim = c(-2, 2), ylim = c(-2, 2)) + 
  ggtitle('Loading plot, all,\ndouble-centered')

# Data for score plots
df_scores_dbc <- data.frame(lygus_meas[, 1:4],
                        scores(lygus_pca_dbc, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_scores_dbc)

## 2a. Score plot, grouped by sex ----
p_scores_G_dbc <- ggplot(df_scores_dbc, aes(x = PC1, y = PC2, color = Gender)) + 
  geom_point(size = 3.5, alpha = 0.6) + 
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle('Score plot, all,\ndouble-centered, by sex') +
  theme(legend.text = element_text(face = "italic"))

# Loading and score plots
both_plots_G_dbc <- plot_grid(p_load_dbc, p_scores_G_dbc, align = "h", rel_widths = c(0.45, 0.57)) 
#ggsave('dbc_all_gender_loadings_ordination.tiff', plot = both_plots_G_dbc, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

## 2b. Score plot, grouped by sex and species ----
p_scores_SG_dbc <- ggplot(df_scores_dbc, aes(x = PC1, y = PC2, color = Species_Gender)) + 
  geom_point(size = 3.5, alpha = 0.6) + 
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle('Score plot, all, double-centered,\nby sex and species')  +
  scale_color_manual(values=palette_SG) +
  theme(legend.text = element_text(face = "italic"))

# Loading and score plots
both_plots_SG_dbc <- plot_grid(p_load_dbc, p_scores_SG_dbc, align = "h", rel_widths = c(0.38, 0.62))
#ggsave('dbc_all_speciesgender_loadings_ordination.tiff', plot = both_plots_SG_dbc, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

### Analysis of variance of PC1 ----
df_SG_dbc <- data.frame(species_gender = lygus_meas$Species_Gender,
                    scores(lygus_pca_dbc, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_SG_dbc)
mod1SG_dbc <- lm(PC1 ~ species_gender, data = df_SG_dbc)
anova(mod1SG_dbc)

# Testing for conditions
mod_diag1SG_dbc <- fortify(mod1SG_dbc)
res_p1SG_dbc <- ggplot(data = mod_diag1SG_dbc, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val1SG_dbc <- mean(mod_diag1SG_dbc$.stdresid)
sd_val1SG_dbc <- sd(mod_diag1SG_dbc$.stdresid)
norm_p1SG_dbc <- ggplot(mod_diag1SG_dbc, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val1SG_dbc, slope = sd_val1SG_dbc)
plot_grid(res_p1SG_dbc, norm_p1SG_dbc, ncol = 2, rel_widths = c(0.55, 0.45))

# Point range graph for PC1
#df_SG$species_gender <- reorder(df$species_gender, df_SG$PC1, FUN=mean)
pc1_plot_SG_dbc<-ggplot(df_SG_dbc, aes(x = species_gender, y = PC1, color = species_gender)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9)) + 
  ggtitle('Point range graph for PC1, all,\ndouble-centered, by sex and species') +
  scale_color_manual(values=palette_SG) + 
  theme(legend.position="none", axis.text.x = element_text(face = "italic"))

# Tukey's test of PC1
TukeyHSD(aov(mod1SG_dbc))

### Analysis of variance of PC2 ----
mod2SG_dbc <- lm(PC2 ~ species_gender, data = df_SG_dbc)
anova(mod2SG_dbc)

# Testing for conditions
mod_diag2SG_dbc <- fortify(mod2SG_dbc)
res_p2SG_dbc <- ggplot(data = mod_diag2SG_dbc, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val2SG_dbc <- mean(mod_diag2SG_dbc$.stdresid)
sd_val2SG_dbc <- sd(mod_diag2SG_dbc$.stdresid)
norm_p2SG_dbc <- ggplot(mod_diag2SG_dbc, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val2SG_dbc, slope = sd_val2SG_dbc)
plot_grid(res_p2SG_dbc, norm_p2SG_dbc, ncol = 2, rel_widths = c(0.55, 0.45))

# Point range graph for PC2
#df_SG$species_gender <- reorder(df_SG$species_gender, df$PC2, FUN=mean)
pc2_plot_SG_dbc<-ggplot(df_SG_dbc, aes(x = species_gender, y = PC2, color = species_gender)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9)) + 
  ggtitle('Point range graph for PC2, all,\ndouble-centered, by sex and species') +
  scale_color_manual(values=palette_SG) +
  theme(legend.text = element_text(face = "italic"), axis.text.x = element_text(face = "italic"))

# Tukey's test of PC2
TukeyHSD(aov(mod2SG_dbc))

both_plots_pc_SG_dbc <- plot_grid(pc1_plot_SG_dbc, pc2_plot_SG_dbc, rel_widths = c(0.40, 0.63))
#ggsave('dbc_all_speciesgender_PC.tiff', plot = both_plots_pc_SG_dbc, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

## 3a. PCA without double-centering data, only females -----
lygus_meas_F <- lygus_meas[lygus_meas$Gender=='female', ]

lygus_pca_F <- rda(lygus_meas_F[, -c(1:4)], scale = TRUE)
summary(lygus_pca_F)
eig <- eigenvals(lygus_pca_F)[1:17]
eig*100/sum(eig)

# Estimating the number of significant PCs using broken-stick model
screeplot(lygus_pca_F, type = "lines", bstick = TRUE)

### Loadings and loading plot ----
scores(lygus_pca_F, display = "species", choices = c(1, 2, 3), scaling = 0)
df_load_F <- as.data.frame(scores(lygus_pca_F, display = "species", choices = c(1, 2, 3), scaling = "species"))

p_load_F <- ggplot(df_load_F) +
  geom_text_repel(aes(x = PC1, y = PC2, label = rownames(df_load_F)),
                  size = 2, max.overlaps = 50, force_pull=2, nudge_x = 0.2, seed = 256) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               color = "coral2", arrow = ar, alpha = 0.7) +
  coord_equal(xlim = c(-0.5, 2.5), ylim = c(-1.5, 1.5)) + 
  ggtitle('Loading plot, females,\nnot double-centered')

# Data for score plots
df_scores_F <- data.frame(lygus_meas_F[, 1:4],
                        scores(lygus_pca_F, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_scores_F)

### Score plot ----
p_scores_F <- ggplot(df_scores_F, aes(x = PC1, y = PC2, color = Species)) + 
  geom_point(size = 3.5, alpha = 0.6) + 
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle('Score plot, females,\nnot double-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic"))

# Loading and score plots
both_plots_F<-plot_grid(p_load_F, p_scores_F, align = "h", rel_widths = c(0.45, 0.67))
#ggsave('no_dbc_female_loadings_ordination.tiff', plot = both_plots_F, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

### Analysis of variance of PC1  ----
df_F <- data.frame(species = lygus_meas_F$Species,
                 scores(lygus_pca_F, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_F)
mod1_F <- lm(PC1 ~ species, data = df_F)
anova(mod1_F)

# Testing for conditions
mod_diag1_F <- fortify(mod1_F)
res_p1_F <- ggplot(data = mod_diag1_F, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val1_F <- mean(mod_diag1_F$.stdresid)
sd_val1_F <- sd(mod_diag1_F$.stdresid)
norm_p1_F <- ggplot(mod_diag1_F, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val1_F, slope = sd_val1_F)
plot_grid(res_p1_F, norm_p1_F, ncol = 2, rel_widths = c(0.55, 0.45))

leveneTest(PC1 ~ species, data = df_F)
ggplot(mod_diag1_F, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag1_F, aes(x = species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for PC1
df_F$species <- reorder(df_F$species, df_F$PC1, FUN=mean)

pc1_plot_F <- ggplot(df_F, aes(x = species, y = PC1, color = species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph for PC1, females,\nnot double-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.position = 'none', axis.text.x = element_text(face = "italic"))

# Tukey's test of PC1
TukeyHSD(aov(mod1_F))

### Analysis of variance of PC2 ----
mod2_F <- lm(PC2 ~ species, data = df_F)
anova(mod2_F)

# Testing for conditions
mod_diag2_F <- fortify(mod2_F)
res_p2_F <- ggplot(data = mod_diag2_F, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val2_F <- mean(mod_diag2_F$.stdresid)
sd_val2_F <- sd(mod_diag2_F$.stdresid)
norm_p2_F <- ggplot(mod_diag2_F, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val2_F, slope = sd_val2_F)
plot_grid(res_p2_F, norm_p2_F, ncol = 2, rel_widths = c(0.55, 0.45))

leveneTest(PC2 ~ species, data = df_F)
ggplot(mod_diag2_F, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag2_F, aes(x = species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for PC2
df_F$species <- reorder(df_F$species, df_F$PC2, FUN=mean)

pc2_plot_F <- ggplot(df_F, aes(x = species, y = PC2, color = species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph for PC2, females,\nnot double-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic"), axis.text.x = element_text(face = "italic"))

# Tukey's test of PC2
TukeyHSD(aov(mod2_F))

both_plots_pc_F<- plot_grid(pc1_plot_F, pc2_plot_F, rel_widths = c(0.40, 0.58))
#ggsave('no_dbc_female_PC.tiff', plot = both_plots_pc_F, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

## 3b.  PCA after double-centering data, only females ----
dbc_lygus_F <- t(apply(lygus_meas_F[, 5:21], 1, center))

lygus_pca_dbc_F <- rda(dbc_lygus_F, scale = TRUE)
summary(lygus_pca_dbc_F)
eig <- eigenvals(lygus_pca_dbc_F)[1:16]
eig*100/sum(eig)

# Estimating the number of significant PCs using broken-stick model
screeplot(lygus_pca_dbc_F, type = "lines", bstick = TRUE)

### Loadings and loading plot ----
scores(lygus_pca_dbc_F, display = "species", choices = c(1, 2, 3), scaling = 0)
df_load_dbc_F <- as.data.frame(scores(lygus_pca_dbc_F, display = "species", choices = c(1, 2, 3), scaling = "species"))

p_load_dbc_F <- ggplot(df_load_dbc_F) +
  geom_text_repel(aes(x = PC1, y = PC2, label = rownames(df_load_dbc_F)),
                  size = 2, max.overlaps = 50, force_pull=0, seed = 100000) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               color = "coral2", arrow = ar, alpha = 0.7) +
  coord_equal(xlim = c(-2, 2), ylim = c(-2, 2)) + 
  ggtitle('Loading plot, females,\ndouble-centered')

# Datafor score plots
df_scores_dbc_F <- data.frame(lygus_meas_F[, 1:4],
                            scores(lygus_pca_dbc_F, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_scores_dbc_F)

### Score plot ----
p_scores_dbc_F <- ggplot(df_scores_dbc_F, aes(x = PC1, y = PC2, color = Species)) + 
  geom_point(size = 3.5, alpha = 0.6) + 
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle('Score plot, females,\ndouble-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic"))

# Loading and score plots
both_plots_F_dbc <- plot_grid(p_load_dbc_F, p_scores_dbc_F, align = "h", rel_widths = c(0.45, 0.67))
#ggsave('dbc_female_loadings_ordination.tiff', plot = both_plots_F_dbc, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

### Analysis of variance of PC1 ----
df_dbc_F <- data.frame(species = lygus_meas_F$Species,
                     scores(lygus_pca_dbc_F, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_dbc_F)
mod1_dbc_F <- lm(PC1 ~ species, data = df_dbc_F)
anova(mod1_dbc_F)

# Testing for conditions
mod_diag1_dbc_F <- fortify(mod1_dbc_F)
res_p1_dbc_F <- ggplot(data = mod_diag1_dbc_F, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val1_dbc_F <- mean(mod_diag1_dbc_F$.stdresid)
sd_val1_dbc_F <- sd(mod_diag1_dbc_F$.stdresid)
norm_p1_dbc_F <- ggplot(mod_diag1_dbc_F, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val1_dbc_F, slope = sd_val1_dbc_F)
plot_grid(res_p1_dbc_F, norm_p1_dbc_F, ncol = 2, rel_widths = c(0.55, 0.45))

leveneTest(PC1 ~ species, data = df_dbc_F)
ggplot(mod_diag1_dbc_F, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag1_dbc_F, aes(x = species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for PC1
df_dbc_F$species <- reorder(df_dbc_F$species, df_dbc_F$PC1, FUN=mean)

pc1_plot_F_dbc <- ggplot(df_dbc_F, aes(x = species, y = PC1, color = species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph for PC1, females,\ndouble-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.position = 'none', axis.text.x = element_text(face = "italic"))

# Tukey's test of PC1
TukeyHSD(aov(mod1_dbc_F))

### Analysis of variance of PC2 ----
mod2_dbc_F <- lm(PC2 ~ species, data = df_dbc_F)
anova(mod2_dbc_F)

# Testing for conditions
mod_diag2_dbc_F <- fortify(mod2_dbc_F)
res_p2_dbc_F <- ggplot(data = mod_diag2_dbc_F, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val2_dbc_F <- mean(mod_diag2_dbc_F$.stdresid)
sd_val2_dbc_F <- sd(mod_diag2_dbc_F$.stdresid)
norm_p2_dbc_F <- ggplot(mod_diag2_dbc_F, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val2_dbc_F, slope = sd_val2_dbc_F)
plot_grid(res_p2_dbc_F, norm_p2_dbc_F, ncol = 2, rel_widths = c(0.55, 0.45))

leveneTest(PC2 ~ species, data = df_dbc_F)
ggplot(mod_diag2_dbc_F, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag2_dbc_F, aes(x = species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for PC2
df_dbc_F$species <- reorder(df_dbc_F$species, df_dbc_F$PC2, FUN=mean)

pc2_plot_F_dbc <- ggplot(df_dbc_F, aes(x = species, y = PC2, color = species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph for PC2, females,\ndouble-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic"), axis.text.x = element_text(face = "italic"))

# Tukey's test of PC2
TukeyHSD(aov(mod2_dbc_F))

both_plots_pc_F_dbc <- plot_grid(pc1_plot_F_dbc, pc2_plot_F_dbc, rel_widths = c(0.41, 0.58))
#ggsave('dbc_female_PC.tiff', plot = both_plots_pc_F_dbc, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

## 4a.  PCA without double-centering data, only males -----
lygus_meas_M <- lygus_meas[lygus_meas$Gender=='male', ]

lygus_pca_M <- rda(lygus_meas_M[, -c(1:4)], scale = TRUE)
summary(lygus_pca_M)
eig <- eigenvals(lygus_pca_M)[1:17]
eig*100/sum(eig)

## Estimating the number of significant PCs using broken-stick model
screeplot(lygus_pca_M, type = "lines", bstick = TRUE)

### Loadings and loading plot for data on males without double-centering ----
scores(lygus_pca_M, display = "species", choices = c(1, 2, 3), scaling = 0)

df_load_M <- as.data.frame(scores(lygus_pca_M, display = "species", choices = c(1, 2, 3), scaling = "species"))

p_load_M <- ggplot(df_load_M) +
  geom_text_repel(aes(x = PC1, y = PC2, label = rownames(df_load_M)),
                  size = 2, max.overlaps = 50, force_pull=2, nudge_x = 0.2, seed = 256) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               color = "coral2", arrow = ar, alpha = 0.7) +
  coord_equal(xlim = c(-0.5, 2.5), ylim = c(-1.5, 1.5)) + 
  ggtitle('Loading plot, males,\nnot double-centered')

# Data for score plots
df_scores_M <- data.frame(lygus_meas_M[, 1:4],
                          scores(lygus_pca_M, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_scores_M)

### Score plot ----
p_scores_M <- ggplot(df_scores_M, aes(x = PC1, y = PC2, color = Species)) + 
  geom_point(size = 3.5, alpha = 0.6) + 
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle('Score plot for data on males\nwithout double-centering') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic"))

# Loading and score plots
both_plots_M <- plot_grid(p_load_M, p_scores_M, align = "h", rel_widths = c(0.45, 0.67))
#ggsave('no_dbc_male_loadings_ordination.tiff', plot = both_plots_M, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

### Analysis of variance of PC1 ----
df_M <- data.frame(species = lygus_meas_M$Species,
                   scores(lygus_pca_M, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_M)
mod1_M <- lm(PC1 ~ species, data = df_M)
anova(mod1_M)

# Testing for conditions
mod_diag1_M <- fortify(mod1_M)
res_p1_M <- ggplot(data = mod_diag1_M, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val1_M <- mean(mod_diag1_M$.stdresid)
sd_val1_M <- sd(mod_diag1_M$.stdresid)
norm_p1_M <- ggplot(mod_diag1_M, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val1_M, slope = sd_val1_M)
plot_grid(res_p1_M, norm_p1_M, ncol = 2, rel_widths = c(0.55, 0.45))

leveneTest(PC1 ~ species, data = df_M)
ggplot(mod_diag1_M, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag1_M, aes(x = species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for PC1
df_M$species <- reorder(df_M$species, df_M$PC1, FUN=mean)

pc1_plot_M <- ggplot(df_M, aes(x = species, y = PC1, color = species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph for PC1, males,\nnot double-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.position = 'none', axis.text.x = element_text(face = "italic"))

# Tukey's test of PC1
TukeyHSD(aov(mod1_M))

### Analysis of variance of PC2 ----
mod2_M <- lm(PC2 ~ species, data = df_M)
anova(mod2_M)

# Testing for conditions
mod_diag2_M <- fortify(mod2_M)
res_p2_M <- ggplot(data = mod_diag2_M, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val2_M <- mean(mod_diag2_M$.stdresid)
sd_val2_M <- sd(mod_diag2_M$.stdresid)
norm_p2_M <- ggplot(mod_diag2_M, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val2_M, slope = sd_val2_M)
plot_grid(res_p2_M, norm_p2_M, ncol = 2, rel_widths = c(0.55, 0.45))

leveneTest(PC2 ~ species, data = df_M)
ggplot(mod_diag2_M, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag2_M, aes(x = species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for PC2
df_M$species <- reorder(df_M$species, df_M$PC2, FUN=mean)

pc2_plot_M <- ggplot(df_M, aes(x = species, y = PC2, color = species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph for PC2, males,\nnot double-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic"), axis.text.x = element_text(face = "italic"))

# Tukey's test of PC2
TukeyHSD(aov(mod2_M))

both_plots_pc_M<- plot_grid(pc1_plot_M, pc2_plot_M, rel_widths = c(0.40, 0.58))
#ggsave('no_dbc_male_PC.tiff', plot = both_plots_pc_M, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

## 4b. PCA after double-centering data, only males ----

dbc_lygus_M <- t(apply(lygus_meas_M[, 5:21], 1, center))

lygus_pca_dbc_M <- rda(dbc_lygus_M, scale = TRUE)
summary(lygus_pca_dbc_M)
eig <- eigenvals(lygus_pca_dbc_M)[1:16]
eig*100/sum(eig)

# Estimating the number of significant PCs using broken-stick model
screeplot(lygus_pca_dbc_M, type = "lines", bstick = TRUE)

### Loadings and loading plot ----
scores(lygus_pca_dbc_M, display = "species", choices = c(1, 2, 3), scaling = 0)

df_load_dbc_M <- as.data.frame(scores(lygus_pca_dbc_M, display = "species", choices = c(1, 2, 3), scaling = "species"))

p_load_dbc_M <- ggplot(df_load_dbc_M) +
  geom_text_repel(aes(x = PC1, y = PC2, label = rownames(df_load_dbc_M)),
                  size = 2, max.overlaps = 50, force_pull=0, seed = 10000) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               color = "coral2", arrow = ar, alpha = 0.7) +
  coord_equal(xlim = c(-2, 2), ylim = c(-2, 2)) + 
  ggtitle('Loading plot, males,\ndouble-centered')

# Data for score plots
df_scores_dbc_M <- data.frame(lygus_meas_M[, 1:4],
                              scores(lygus_pca_dbc_M, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_scores_dbc_M)

### Score plot ----
p_scores_dbc_M <- ggplot(df_scores_dbc_M, aes(x = PC1, y = PC2, color = Species)) + 
  geom_point(size = 3.5, alpha = 0.6) + 
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle('Score plot, males,\ndouble-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic"))

# Loading and score plots
both_plots_M_dbc <- plot_grid(p_load_dbc_M, p_scores_dbc_M, align = "h", rel_widths = c(0.45, 0.67))
#ggsave('dbc_male_loadings_ordination.tiff', plot = both_plots_M_dbc, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

### Analysis of variance of PC1 ----
df_dbc_M <- data.frame(species = lygus_meas_M$Species,
                       scores(lygus_pca_dbc_M, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
str(df_dbc_M)
mod1_dbc_M <- lm(PC1 ~ species, data = df_dbc_M)
anova(mod1_dbc_M)

# Testing for conditions
mod_diag1_dbc_M <- fortify(mod1_dbc_M)
res_p1_dbc_M <- ggplot(data = mod_diag1_dbc_M, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val1_dbc_M <- mean(mod_diag1_dbc_M$.stdresid)
sd_val1_dbc_M <- sd(mod_diag1_dbc_M$.stdresid)
norm_p1_dbc_M <- ggplot(mod_diag1_dbc_M, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val1_dbc_M, slope = sd_val1_dbc_M)
plot_grid(res_p1_dbc_M, norm_p1_dbc_M, ncol = 2, rel_widths = c(0.55, 0.45))

leveneTest(PC1 ~ species, data = df_dbc_M)
ggplot(mod_diag1_dbc_M, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag1_dbc_M, aes(x = species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for PC1
df_dbc_M$species <- reorder(df_dbc_M$species, df_dbc_M$PC1, FUN=mean)

pc1_plot_M_dbc <- ggplot(df_dbc_M, aes(x = species, y = PC1, color = species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph for PC1, males,\ndouble-centered') +
  theme(legend.position = 'none', axis.text.x = element_text(face = "italic")) +
  scale_color_manual(values=palette_S)

# Tukey's test of PC1
TukeyHSD(aov(mod1_dbc_M))

### Analysis of variance of PC2 ----
mod2_dbc_M <- lm(PC2 ~ species, data = df_dbc_M)
anova(mod2_dbc_M)

# Testing for conditions
mod_diag2_dbc_M <- fortify(mod2_dbc_M)
res_p2_dbc_M <- ggplot(data = mod_diag2_dbc_M, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val2_dbc_M <- mean(mod_diag2_dbc_M$.stdresid)
sd_val2_dbc_M <- sd(mod_diag2_dbc_M$.stdresid)
norm_p2_dbc_M <- ggplot(mod_diag2_dbc_M, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val2_dbc_M, slope = sd_val2_dbc_M)
plot_grid(res_p2_dbc_M, norm_p2_dbc_M, ncol = 2, rel_widths = c(0.55, 0.45))

leveneTest(PC2 ~ species, data = df_dbc_M)
ggplot(mod_diag2_dbc_M, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag2_dbc_M, aes(x = species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for PC2
df_dbc_M$species <- reorder(df_dbc_M$species, df_dbc_M$PC2, FUN=mean)

pc2_plot_M_dbc <- ggplot(df_dbc_M, aes(x = species, y = PC2, color = species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph for PC2, males,\ndouble-centered') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic"), axis.text.x = element_text(face = "italic"))

# Tukey's test of PC2
TukeyHSD(aov(mod2_dbc_M))

both_plots_pc_M_dbc <- plot_grid(pc1_plot_M_dbc, pc2_plot_M_dbc, rel_widths = c(0.40, 0.58))
#ggsave('dbc_male_PC.tiff', plot = both_plots_pc_M_dbc, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

## 5. Functions from Baur & Leuenberger (2011) ----

## 5a. PCA in isometry free shape space ---------
shapePCA <- function(x, npc=3, rpc=4){
  X <- log(x)
  X <- scale(X, center=TRUE, scale=FALSE)
  p <- dim(X)[2]
  I <- diag(1,p,p)
  a0 <- as.vector(rep(1,p)/p)
  P <- I-p*(a0%*%t(a0))
  Y <- X%*%P
  colnames(Y) <- colnames(X)
  isosize <- X%*%a0; colnames(isosize) <- "isosize"
  shapePCA <- prcomp(Y,center=FALSE, scale=FALSE)
  loadings <- shapePCA$rotation[,1:p]
  PCmatrix <- shapePCA$x
  components.a <- shapePCA$sdev^2; components.b <- components.a/sum(components.a); components.c <- cumsum(components.b); components <- rbind(components.a, components.b, components.c)
  colnames(loadings) <- paste("shape.PC", 1:p, sep="")
  colnames(PCmatrix) <- colnames(loadings)
  colnames(components) <- colnames(loadings); rownames(components) <- c("Variance", "Proportion of Variance", "Cumulative Proportion")
  list(shapePCA=shapePCA, PCmatrix=PCmatrix[,1:npc], isosize=isosize, loadings=round(loadings[,1:npc], rpc), components=round(components[,1:npc], rpc))
}

pca_M <- shapePCA(lygus_meas_M[, -c(1:4)], npc=3)
pca_F <- shapePCA(lygus_meas_F[, -c(1:4)], npc=3)
pca_M$loadings
pca_F$loadings
pca_M$components
pca_F$components
rot_M <- pca_M$PCmatrix
rot_F <- pca_F$PCmatrix

ifs_plot_1_M <- ggplot(lygus_meas_M, aes(x = pca_M$isosize, y = rot_M[,1], color = Species)) +
  geom_point(size = 3.5, alpha = 0.8) + 
  ggtitle('Isometric size versus first shape component, males') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic")) +
  xlab("isometric size") + ylab("shape PC1")
#ggsave('ifs_plot_1_M.tiff', plot = ifs_plot_1_M, path = './figs', scale = 1.5, width = 20, device='tiff', dpi=300,  units = c("cm"))

ifs_plot_1_F <- ggplot(lygus_meas_F, aes(x = pca_F$isosize, y = rot_F[,1], color = Species)) +
  geom_point(size = 3.5, alpha = 0.8) + 
  ggtitle('Isometric size versus first shape component, females') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic")) +
  xlab("isometric size") + ylab("shape PC1")
#ggsave('ifs_plot_1_F', plot = ifs_plot_1_F, path = './figs', scale = 1.5, width = 20, device='tiff', dpi=300,  units = c("cm"))

ifs_plot_2_M <- ggplot(lygus_meas_M, aes(x = pca_M$isosize, y = rot_M[,2], color = Species)) +
  geom_point(size = 3.5, alpha = 0.8) + 
  ggtitle('Isometric size versus second shape component, males') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic")) +
  xlab("isometric size") + ylab("shape PC2")
#ggsave('ifs_plot_2_M', plot = ifs_plot_2_M, path = './figs', scale = 1.5, width = 20, device='tiff', dpi=300,  units = c("cm"))

ifs_plot_2_F <- ggplot(lygus_meas_F, aes(x = pca_F$isosize, y = rot_F[,2], color = Species)) +
  geom_point(size = 3.5, alpha = 0.8) + 
  ggtitle('Isometric size versus second shape component, females') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic")) +
  xlab("isometric size") + ylab("shape PC2")
#ggsave('ifs_plot_2_F', plot = ifs_plot_2_F, path = './figs', scale = 1.5, width = 20, device='tiff', dpi=300,  units = c("cm"))

ifs_plot_12_M <- ggplot(lygus_meas_M, aes(x = rot_M[,1], y = rot_M[,2], color = Species)) +
  geom_point(size = 3.5, alpha = 0.8) + 
  ggtitle('First versus second shape component, males') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic")) +
  xlab("shape PC1") + ylab("shape PC2")
#ggsave('ifs_plot_12_M', plot = ifs_plot_12_M, path = './figs', scale = 1.5, width = 20, device='tiff', dpi=300,  units = c("cm"))

ifs_plot_12_F <- ggplot(lygus_meas_F, aes(x = rot_F[,1], y = rot_F[,2], color = Species)) +
  geom_point(size = 3.5, alpha = 0.8) + 
  ggtitle('First versus second shape component, females') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic")) +
  xlab("shape PC1") + ylab("shape PC2")
#ggsave('ifs_plot_12_F', plot = ifs_plot_12_F, path = './figs', scale = 1.5, width = 20, device='tiff', dpi=300,  units = c("cm"))

### Analysis of variance of isometric size of data on males----
mod_size_M <- lm(pca_M$isosize ~ Species, data = lygus_meas_M)
anova(mod_size_M)

# Testing for conditions
mod_diag_size_M <- fortify(mod_size_M)
res_p_size_M <- ggplot(data = mod_diag_size_M, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val_size_M <- mean(mod_diag_size_M$.stdresid)
sd_val_size_M <- sd(mod_diag_size_M$.stdresid)
norm_p_size_M <- ggplot(mod_diag_size_M, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val_size_M, slope = sd_val_size_M)
plot_grid(res_p_size_M, norm_p_size_M, ncol = 2, rel_widths = c(0.55, 0.45))

ggplot(mod_diag_size_M, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag_size_M, aes(x = Species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for isometric size
lygus_meas_M$Species <- reorder(lygus_meas_M$Species, pca_M$isosize, FUN=mean)

size_plot_M <- ggplot(lygus_meas_M, aes(x = Species, y = pca_M$isosize, color = Species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph\nfor isometric size, males') +
  theme(legend.position = 'none', axis.text.x = element_text(face = "italic")) +
  scale_color_manual(values=palette_S)

# Tukey's test
TukeyHSD(aov(mod_size_M))

### Analysis of variance of isometric size of data on females ----
mod_size_F <- lm(pca_F$isosize ~ Species, data = lygus_meas_F)
anova(mod_size_F)

# Testing for conditions
mod_diag_size_F <- fortify(mod_size_F)
res_p_size_F <- ggplot(data = mod_diag_size_F, aes(x = .fitted, y = .stdresid)) + 
  geom_point(aes(size = .cooksd)) + 
  geom_hline(yintercept = 0) + 
  geom_smooth(method="loess", se=FALSE)
mean_val_size_F <- mean(mod_diag_size_F$.stdresid)
sd_val_size_F <- sd(mod_diag_size_F$.stdresid)
norm_p_size_F <- ggplot(mod_diag_size_F, aes(sample = .stdresid)) + 
  geom_point(stat = "qq") + 
  geom_abline(intercept = mean_val_size_F, slope = sd_val_size_F)
plot_grid(res_p_size_F, norm_p_size_F, ncol = 2, rel_widths = c(0.55, 0.45))

ggplot(mod_diag_size_F, aes(x = .fitted, y = .stdresid)) + geom_point() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)
ggplot(mod_diag_size_F, aes(x = Species, y = .stdresid)) + geom_boxplot() + ggtitle("Residuals plot of model predictors") + geom_hline(yintercept = 0)

# Point range graph for isometric size
lygus_meas_F$Species <- reorder(lygus_meas_F$Species, pca_F$isosize, FUN=mean)

size_plot_F <- ggplot(lygus_meas_F, aes(x = Species, y = pca_F$isosize, color = Species)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  ggtitle('Point range graph\nfor isometric size, females') +
  scale_color_manual(values=palette_S) +
  theme(legend.text = element_text(face = "italic"), axis.text.x = element_text(face = "italic"))

# Tukey's test
TukeyHSD(aov(mod_size_F))

both_size_plots <- plot_grid(size_plot_M, size_plot_F, rel_widths = c(0.40, 0.58))
#ggsave('both_size_plots.tiff', plot = both_size_plots, path = './figs', scale = 1.5, width = 20, height = 10, device='tiff', dpi=300,  units = c("cm"))

## 5b. PCA Ratio Spectrum ---------

pcaRS <- function (x, pc=1, bootrep=500, barcol="blue", barlwd=2.5, linecol="black", linelwd=0.5, labelsize=0.9, labelfont=1, nosize=0.7, nofont=1, maina="PCA Ratio Spectrum for PC", mainb=pc, mainc="", suba="bars = 68% confidence intervals based on ", subb=bootrep, subc=" bootstrap replicates", main=paste(maina, mainb, mainc, sep=""), sub=paste(suba, subb, subc, sep="")) {
  trait.names <- colnames(x)
  Y <- log(x)
  Y <- scale(Y, center=TRUE, scale=FALSE)
  Sigma <- cov(Y)
  p <- length(trait.names)
  n <- nrow(Y)
  I <- diag(1, p, p)
  a0 <- as.vector(rep(1,p)/p)
  P <- I-a0%*%solve(t(a0)%*%a0)%*%t(a0)
  Sigma.1 <- P%*%Sigma%*%P
  u <- eigen(Sigma.1)$vectors[, pc]
  rep.boot=bootrep
  U.boot <- matrix(0, nrow=rep.boot, ncol=p)
  for (i in c(1:rep.boot)){
    index.boot=sample(c(1:n), n, replace=TRUE)
    Y.boot <- Y[index.boot, ]
    Sigma.boot <- P%*%cov(Y.boot)%*%P
    u.boot <- eigen(Sigma.boot)$vectors[, pc]
    if (u%*%u.boot<0){
      u.boot <- -u.boot
    }
    U.boot[i, ] <- u.boot
  }
  mean.u <- apply(U.boot, 2, mean)
  sd.u  <- apply(U.boot, 2, sd)
  mean.minus.sd <- mean.u-sd.u
  mean.plus.sd <- mean.u+sd.u
  m <- min(mean.u); M <- max(mean.u)
  u.sorted <- sort(mean.u, decreasing=TRUE)
  index <- sort(mean.u, index.return=TRUE, decreasing=TRUE)$ix
  trait.names.sorted <- trait.names[index]
  mean.minus.sd <- mean.minus.sd[index]
  mean.plus.sd <- mean.plus.sd[index]
  plot(c(0, 0), c(m, M), xlab="", ylab="", type="n", cex=1, asp=1, main=main, sub=sub, cex.lab=1, axes=FALSE)
  lines(c(0, 0), c(m, M), col=linecol, lwd=linelwd)
  text(0, u.sorted[1], labels=signif(u.sorted[1], digits=2), pos=3, srt=0, col="black", font=nofont, cex=nosize)
  text(0, u.sorted[p], labels=signif(u.sorted[p], digits=2), pos=1, srt=0, col="black", font=nofont, cex=nosize)
  for (k in 1:p){
    lines(c(mean.minus.sd[k]-u.sorted[k], mean.plus.sd[k]-u.sorted[k]), c(u.sorted[k], u.sorted[k]), col=barcol, lwd=barlwd)
    if(k%%2==0){
      text(max(sd.u), u.sorted[k], labels=trait.names.sorted[k], pos=4, srt=0, col="black", font=labelfont, cex=labelsize)
    }
    if(k%%2==1){
      text(-max(sd.u), u.sorted[k], labels=trait.names.sorted[k], pos=2, srt=0, col="black", font=labelfont, cex=labelsize)
    }}
}

dev.off()
pcaRS(lygus_meas_M[, -c(1:4)], pc=1, mainc=", Lygus males")
pcaRS(lygus_meas_M[, -c(1:4)], pc=2, mainc=", Lygus males")

pcaRS(lygus_meas_F[, -c(1:4)], pc=1, mainc=", Lygus females")
pcaRS(lygus_meas_F[, -c(1:4)], pc=2, mainc=", Lygus females")


## 5c. LDA Ratio Extractor for paires ----

ldaRE <- function(x, g, rd=1, kmax=2) {
  
  
  X <- data.frame(x)
  groups <- g
  species.names <- levels(groups)
  nb.species <- length(species.names)
  
  trait.names <- colnames(X)
  p <- length(trait.names)
  X1 <- (X[groups==species.names[1],1:p]) 
  X2 <- (X[groups==species.names[2],1:p])
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  Y1 <- as.matrix(log(X1))
  Y2 <- as.matrix(log(X2))
  m1 <- colMeans(Y1)
  m2 <- colMeans(Y2)
  m12 <- as.vector(m1-m2)
  m12 <- m12/sqrt(sum(m12^2))
  Y.pooled <- rbind(scale(Y1,center=TRUE,scale=FALSE),scale(Y2,center=TRUE,scale=FALSE))
  Sigma <- cov(Y.pooled)
  w <- solve(Sigma)%*%m12
  standard.distance <- function(x,y){
    x <- as.vector(x)
    y <- as.vector(y)
    D <- abs(mean(x)-mean(y))/(sd(c(scale(x,center=TRUE,scale=FALSE),c(scale(y,center=TRUE,scale=FALSE)))))
    return(D)
  }
  D.tot <- standard.distance(Y1%*%w,Y2%*%w)
  I <- diag(1,p,p)
  a0 <- as.vector(rep(1,p)/p)
  M.k <- matrix(a0,ncol=1)
  k.max <- kmax
  best.ratio <- matrix(rep(0,2*k.max),ncol=2)
  for (k in c(1:k.max)){
    P.k <- I-M.k%*%solve(t(M.k)%*%M.k)%*%t(M.k)
    Sigma.k <- P.k%*%Sigma%*%P.k
    Sigma.k <- 0.5*(Sigma.k+t(Sigma.k))
    ev <- eigen(Sigma.k)$values[1:(p-k)]
    V <- eigen(Sigma.k)$vectors
    Lambda.plus <- diag(c(1/ev,rep(0,k))) 
    Sigma.plus <- V%*%Lambda.plus%*%t(V)
    w.k <- Sigma.plus%*%P.k%*%m12
    w.k <- as.matrix(w.k,ncol=1)
    if (k==rd){ 
      D.bij=standard.distance(Y1%*%w.k,Y2%*%w.k)
      D.shape <- standard.distance(Y1%*%w.k,Y2%*%w.k)/D.tot
      D.size <- standard.distance(Y1%*%a0,Y2%*%a0)/D.tot
      delta <- D.size/(D.size+D.shape)
    }
    c.max <- 0
    for (i in c(1:(p-1))){  
      for (j in c((1+i):p)){
        b.ij <- rep(0,p)
        b.ij[i] <- 1
        b.ij[j] <- -1
        b.ij <- as.matrix(b.ij,ncol=1)  
        c.ij <- (t(b.ij)%*%Sigma%*%w.k)^2/(t(b.ij)%*%Sigma%*%b.ij)
        if (c.ij>c.max){
          best.b.ij <- b.ij 
          c.max <- c.ij
          best.ratio[k,1] <- i
          best.ratio[k,2] <- j 
        }    
      }}
    print(paste("Ratio no. ", as.character(k), ": ", trait.names[best.ratio[k,1]], "/", trait.names[best.ratio[k,2]], sep=""))
    M.k <- cbind(M.k,Sigma%*%best.b.ij)
  }
  D.bij <- paste("D.bij of ratio no. ", as.character(rd), ": ", trait.names[best.ratio[rd,1]], "/", trait.names[best.ratio[rd,2]], " = ", as.character(D.bij), sep="")
  D.shape <- paste("D.shape of ratio no. ", as.character(rd), ": ", trait.names[best.ratio[rd,1]], "/", trait.names[best.ratio[rd,2]], " = ", as.character(D.shape), sep="")
  D.size <- paste("D.size of ratio no. ", as.character(rd), ": ", trait.names[best.ratio[rd,1]], "/", trait.names[best.ratio[rd,2]], " = ", as.character(D.size), sep="")
  delta <- paste("delta of ratio no. ", as.character(rd), ": ", trait.names[best.ratio[rd,1]], "/", trait.names[best.ratio[rd,2]], " = ", as.character(delta), sep="")
  list(D.bij=D.bij, D.shape=D.shape, D.size=D.size, delta=delta)
}

sort_pratensis_gemellatus <- lygus_meas_F[, 1]=='Lygus pratensis'|lygus_meas_F[, 1]=='Lygus gemellatus'

ldaRE(lygus_meas_F[sort_pratensis_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_pratensis_gemellatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_pratensis_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_pratensis_gemellatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_pratensis_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_pratensis_gemellatus, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_pratensis_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_pratensis_gemellatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_pratensis_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_pratensis_gemellatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_pratensis_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_pratensis_gemellatus, 1]), rd=3, kmax=3)


sort_punctatus_gemellatus <- lygus_meas_F[, 1]=='Lygus punctatus'|lygus_meas_F[, 1]=='Lygus gemellatus'

ldaRE(lygus_meas_F[sort_punctatus_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_punctatus_gemellatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_punctatus_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_punctatus_gemellatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_punctatus_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_punctatus_gemellatus, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_punctatus_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_punctatus_gemellatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_punctatus_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_punctatus_gemellatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_punctatus_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_punctatus_gemellatus, 1]), rd=3, kmax=3)


sort_rugulipennis_gemellatus <- lygus_meas_F[, 1]=='Lygus rugulipennis'|lygus_meas_F[, 1]=='Lygus gemellatus'

ldaRE(lygus_meas_F[sort_rugulipennis_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_rugulipennis_gemellatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_rugulipennis_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_rugulipennis_gemellatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_rugulipennis_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_rugulipennis_gemellatus, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_rugulipennis_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_rugulipennis_gemellatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_rugulipennis_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_rugulipennis_gemellatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_rugulipennis_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_rugulipennis_gemellatus, 1]), rd=3, kmax=3)


sort_wagneri_gemellatus <- lygus_meas_F[, 1]=='Lygus wagneri'|lygus_meas_F[, 1]=='Lygus gemellatus'

ldaRE(lygus_meas_F[sort_wagneri_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_gemellatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_wagneri_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_gemellatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_wagneri_gemellatus, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_gemellatus, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_wagneri_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_gemellatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_wagneri_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_gemellatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_wagneri_gemellatus, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_gemellatus, 1]), rd=3, kmax=3)


sort_punctatus_pratensis <- lygus_meas_F[, 1]=='Lygus punctatus'|lygus_meas_F[, 1]=='Lygus pratensis'
ldaRE(lygus_meas_F[sort_punctatus_pratensis, -c(1:4)], droplevels(lygus_meas_F[sort_punctatus_pratensis, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_punctatus_pratensis, -c(1:4)], droplevels(lygus_meas_F[sort_punctatus_pratensis, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_punctatus_pratensis, -c(1:4)], droplevels(lygus_meas_F[sort_punctatus_pratensis, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_punctatus_pratensis, -c(1:4)], droplevels(lygus_meas_M[sort_punctatus_pratensis, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_punctatus_pratensis, -c(1:4)], droplevels(lygus_meas_M[sort_punctatus_pratensis, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_punctatus_pratensis, -c(1:4)], droplevels(lygus_meas_M[sort_punctatus_pratensis, 1]), rd=3, kmax=3)


sort_rugulipennis_pratensis <- lygus_meas_F[, 1]=='Lygus rugulipennis'|lygus_meas_F[, 1]=='Lygus pratensis'

ldaRE(lygus_meas_F[sort_rugulipennis_pratensis, -c(1:4)], droplevels(lygus_meas_F[sort_rugulipennis_pratensis, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_rugulipennis_pratensis, -c(1:4)], droplevels(lygus_meas_F[sort_rugulipennis_pratensis, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_rugulipennis_pratensis, -c(1:4)], droplevels(lygus_meas_F[sort_rugulipennis_pratensis, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_rugulipennis_pratensis, -c(1:4)], droplevels(lygus_meas_M[sort_rugulipennis_pratensis, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_rugulipennis_pratensis, -c(1:4)], droplevels(lygus_meas_M[sort_rugulipennis_pratensis, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_rugulipennis_pratensis, -c(1:4)], droplevels(lygus_meas_M[sort_rugulipennis_pratensis, 1]), rd=3, kmax=3)


sort_wagneri_pratensis <- lygus_meas_F[, 1]=='Lygus wagneri'|lygus_meas_F[, 1]=='Lygus pratensis'

ldaRE(lygus_meas_F[sort_wagneri_pratensis, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_pratensis, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_wagneri_pratensis, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_pratensis, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_wagneri_pratensis, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_pratensis, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_wagneri_pratensis, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_pratensis, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_wagneri_pratensis, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_pratensis, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_wagneri_pratensis, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_pratensis, 1]), rd=3, kmax=3)


sort_rugulipennis_punctatus <- lygus_meas_F[, 1]=='Lygus rugulipennis'|lygus_meas_F[, 1]=='Lygus punctatus'

ldaRE(lygus_meas_F[sort_rugulipennis_punctatus, -c(1:4)], droplevels(lygus_meas_F[sort_rugulipennis_punctatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_rugulipennis_punctatus, -c(1:4)], droplevels(lygus_meas_F[sort_rugulipennis_punctatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_rugulipennis_punctatus, -c(1:4)], droplevels(lygus_meas_F[sort_rugulipennis_punctatus, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_rugulipennis_punctatus, -c(1:4)], droplevels(lygus_meas_M[sort_rugulipennis_punctatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_rugulipennis_punctatus, -c(1:4)], droplevels(lygus_meas_M[sort_rugulipennis_punctatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_rugulipennis_punctatus, -c(1:4)], droplevels(lygus_meas_M[sort_rugulipennis_punctatus, 1]), rd=3, kmax=3)


sort_wagneri_punctatus <- lygus_meas_F[, 1]=='Lygus wagneri'|lygus_meas_F[, 1]=='Lygus punctatus'

ldaRE(lygus_meas_F[sort_wagneri_punctatus, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_punctatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_wagneri_punctatus, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_punctatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_wagneri_punctatus, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_punctatus, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_wagneri_punctatus, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_punctatus, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_wagneri_punctatus, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_punctatus, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_wagneri_punctatus, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_punctatus, 1]), rd=3, kmax=3)


sort_wagneri_rugulipennis <- lygus_meas_F[, 1]=='Lygus wagneri'|lygus_meas_F[, 1]=='Lygus rugulipennis'

ldaRE(lygus_meas_F[sort_wagneri_rugulipennis, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_rugulipennis, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_F[sort_wagneri_rugulipennis, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_rugulipennis, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_F[sort_wagneri_rugulipennis, -c(1:4)], droplevels(lygus_meas_F[sort_wagneri_rugulipennis, 1]), rd=3, kmax=3)

ldaRE(lygus_meas_M[sort_wagneri_rugulipennis, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_rugulipennis, 1]), rd=1, kmax=3)
ldaRE(lygus_meas_M[sort_wagneri_rugulipennis, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_rugulipennis, 1]), rd=2, kmax=3)
ldaRE(lygus_meas_M[sort_wagneri_rugulipennis, -c(1:4)], droplevels(lygus_meas_M[sort_wagneri_rugulipennis, 1]), rd=3, kmax=3)