####### Quantifying Genetic Variation in Reproductive Traits in Scots Pine #######

####### R code for data analysis and figure production ########

library(dplyr)
library(lme4)
library(lmeTest)
library(ggplot2)

cones<-read.csv('cone_data.csv')

#### Preliminary Correlations ###

cone_corrs<-read.csv('Preliminary_Traits.csv')
cone_dat<-as.matrix(cone_corrs)
data <- na.omit(cone_corrs)

library(Hmisc)

pcc<-rcorr(cone_dat, type=c("pearson"))
p<-pcc$P
r<-pcc$r
write.table(p, "/Users/jessshortman/Desktop/Dissertation/Neat Correlations/p_cone.txt", sep="\t")
write.table(r, "/Users/jessshortman/Desktop/Dissertation/Neat Correlations/r_cone.txt", sep="\t")

#### Principal Component Analysis ####

library(devtools)

# Add region column
cones <- cones %>%
  mutate(Region = case_when(
    pop %in% c("AM", "BW", "GA", "GC", "GE", "GL", "MG", "RD", 'SO') ~ "Central Atlantic",
    pop %in% c("BE", "CC", "CG", "CR", "LC", "SD") ~ "Hyper-Oceanic",
    pop %in% c("AB", "AC", "BB","GD", "GT", "RM") ~ "Cairngorms"
  ))

# Subset for traits 
cones_pca<-subset(cones[,1:20])
any(is.na(cones_pca))
cones_pca<-na.omit(cones_pca)
head(cones_pca)
cones_pca<-subset(cones_pca[,-2])

# Scale traits 
traits_scaled <- scale(cones_pca)

# PCA analysis 
pca_result <- prcomp(traits_scaled, center = TRUE, scale. = TRUE)
summary(pca_result)
pca_result

pca_df<-as.data.frame(pca_result$x)

# Add family, population and region back to PCA results 
cones_rm<-subset(cones[,1:20])
cones_rm<-na.omit(cones_rm)
cones_rm <- cones_rm %>%
  mutate(Region = case_when(
    pop %in% c("AM", "BW", "GA", "GC", "GE", "GL", "MG", "RD", 'SO') ~ "Central Atlantic",
    pop %in% c("BE", "CC", "CG", "CR", "LC", "SD") ~ "Hyper-Oceanic",
    pop %in% c("AB", "AC", "BB","GD", "GT", "RM") ~ "Cairngorms"
  ))
pca_df$Family<-cones_rm$fam
pca_df$Population<-cones_rm$pop
pca_df$Region<-cones_rm$Region

# Group by family 
family_pca_region <- pca_df %>%
  group_by(Family, Region, Population) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = 'drop'
  )

# Assign colour to each region
colour<-c('Cairngorms'= '#619CFF', "Central Atlantic"= '#F8766D', 'Hyper-Oceanic' =  '#F7C300')
family_pca_region <- family_pca_region %>%
  mutate(colour = case_when(
    Region %in% c("Cairngorms") ~ "#619CFF",
    Region %in% c("Central Atlantic") ~ "#F8766D",
    Region %in% c("Hyper-Oceanic") ~ "#F7C300"
  ))  

library(ggpubr)
library(ggforce)
library(cluster)
library(ggfortify)

# PCA plot clustered by region
ggplot(family_pca_region, aes(x = PC1, y = PC2, colour = Region, fill=colour)) +
  geom_point(size = 2) +
  stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95)+
  scale_fill_identity()+
  labs(x = "Principal Component 1 (22%)",
       y = "Principal Component 2 (14%)",
       color = "Region") +
  theme_minimal() +
  scale_colour_manual(values=colour)


#### Genetic Variation ####

### Cone length ###
model_cone_length<-lmer(cone_length_mm~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones)
summary(model_cone_length)

### Stalk length ###
model_stalk_length<-lmer(stalk_length_mm~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones)
summary(model_stalk_length)

### Scale length ###
model_scale_length<-lmer(scale_length_mm~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones)
summary(model_scale_length)

### Cone angle ###
model_cone_angle<-lmer(cone_angle~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones)
summary(model_cone_angle)

### Hue ###
model_hue<-lmer(colour..hue. ~ (1|pop) + (1|fam) + (1|block)+ (1|Order), data=cones)
summary(model_hue)

### Value ###
model_value<-lmer(value ~ (1|pop) + (1|fam) + (1|block)+ (1|Order), data=cones)
summary(model_value)
 
### Chroma ###
model_chroma<-lmer(chroma ~ (1|pop) + (1|fam) + (1|block)+ (1|Order), data=cones)
summary(model_chroma)

### Profile ###
model_profile<-lmer(profile ~ (1|pop) + (1|fam) + (1|block)+ (1|Order), data=cones)
summary(model_profile)

### Openness ###
model_openness<-lmer(openness ~ (1|pop) + (1|fam) + (1|block)+ (1|Order), data=cones)
summary(model_openness)

### Cone Ratio ###
model_cone_ratio<-lmer(cone_ratio ~ (1|pop) + (1|fam) + (1|block)+ (1|Order), data=cones)
summary(model_cone_ratio)

### Seed Viability ###
viability<-read.csv('viability.csv')

model_viability<-lmer(seed_viability ~ (1|fam)+ (1|pop)+(1|block)+ (1|Order), data=viability)
summary(model_viability)

model_viability_b <- glmer(seed_viability~ (1|fam)+ (1|pop)+(1|block)+ (1|Order), data=viability, family = binomial)
summary(model_viability)

AIC(model_viability, model_viability_b)

### Scale Number ###
model_scale_number<-lmer(scale_number ~ (1|fam)+ (1|pop)+(1|block)+ (1|id), data=cones)
summary(model_scale_number)

### Seed number ##
model_seed_number<-lmer(seed_number ~ (1|fam)+ (1|pop)+(1|block)+ (1|Order), data=viability)
summary(model_seed_number)

### Cone Number ###
library(insight)

subset_data <- cones %>%
  select(block, pop, fam, cone_number, height_m, growth_rate, long) %>%  
  distinct(block, pop, fam, .keep_all = TRUE)  

subset_data$log_cones<- log(subset_data$cone_number)
subset_data <- subset_data[is.finite(subset_data$log_cones), ]

model_cone_number <- lmer(cone_number ~ (1|fam) + (1|pop) +(1|block) , data = subset_data)
summary(model_cone_number)

model_cone_number_p <- glmer(cone_number ~ (1|fam) + (1|pop) +(1|block) , family='poisson',data = subset_data)
summary(model_cone_number_p)

model_cone_number_l <- lmer(log(cone_number) ~ (1|fam) + (1|pop) +(1|block) , data = subset_data)
summary(model_cone_number_l)

model_cone_number_nb<-glmer.nb(cone_number ~ (1|fam) + (1|pop) +(1|block) , data = subset_data)
summary(model_cone_number_nb)

AIC(model_cone_number,model_cone_number_p, model_cone_number_l, model_cone_number_nb)

#### Heritability ####

### Bar plot ###

# Create data frame 
trait<-c('Cone Count','Cone Length', 'Stalk Length', 'Scale Length', 'Scale Count',
         'Ratio', 'Angle', 'Hue', 'Value', 'Chroma', 'Profile', 'Openness', 'Seed Count', 'Seed Viability') 

half_sib<-c(0.36,0.63,0.39,0.39,0.50,0.68,0.00,0.12,0.07,0.06,0.59,0.46,0.29,0.27)

heritability<-data.frame(trait,half_sib)
heritability$half_sib_se<-c(0.15,0.16,0.15,0.15,0.15,0.16,0.13,0.13,
                            0.13,0.13,0.16,0.26,0.25,0.25)

heritability$admixture<-c(0.27,
                          0.47,
                          0.30,
                          0.29,
                          0.38,
                          0.51,
                          0.00,
                          0.09,
                          0.05,
                          0.05,
                          0.44,
                          0.34,
                          0.22,
                          0.20)

heritability$admixture_se<-c(0.11,
                             0.12,
                             0.11,
                             0.11,
                             0.11,
                             0.12,
                             0.09,
                             0.10,
                             0.10,
                             0.10,
                             0.12,
                             0.20,
                             0.19,
                             0.19)

heritability$full_sib<-c(0.18,
                         0.31,
                         0.20,
                         0.20,
                         0.25,
                         0.34,
                         0.00,
                         0.06,
                         0.03,
                         0.03,
                         0.30,
                         0.23,
                         0.14,
                         0.14)

heritability$full_sib_se<-c(0.07,
                            0.08,
                            0.07,
                            0.07,
                            0.08,
                            0.08,
                            0.06,
                            0.07,
                            0.07,
                            0.07,
                            0.08,
                            0.13,
                            0.12,
                            0.12)

heritability<-rename(heritability, 'Trait' = trait)

library(tidyverse)
long_data <- heritability %>%
  pivot_longer(cols = c(full_sib, half_sib, admixture),
               names_to = "Heritability_Type",
               values_to = "Heritability_Value") %>%
  pivot_longer(cols = c(full_sib_se, half_sib_se, admixture_se),
               names_to = "SE_Type",
               values_to = "SE_Value") %>%
  filter(substr(Heritability_Type, 1, 3) == substr(SE_Type, 1, 3))

# Order traits 
trait_order<-heritability$Trait
long_data$Trait <- factor(long_data$Trait, levels = trait_order)

# Order relatedness scenarios 
heritability_order<-c('half_sib', 'admixture', 'full_sib')
long_data$Heritability_Type<-factor(long_data$Heritability_Type, levels=heritability_order)

# Make bar plot 
ggplot(data=long_data,
       aes(x=Trait, y=Heritability_Value, fill=Heritability_Type))+
  scale_fill_manual(values = c('#619CFF', '#F8766D', '#F7C300'))+
  geom_bar(stat='identity', position='dodge')+
  geom_errorbar(aes(ymin=pmax(Heritability_Value - SE_Value, 0), ymax=Heritability_Value+SE_Value), width=0.5,
                position=position_dodge(.9)) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"))+
  labs(x='Trait', y= "Narrow-Sense Heritability", fill='Relatedness')+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, size = 10))

### Parent-Offspring Regressions ###

# read in maternal traits 
mat_cones<-read.table('8a86cb01-46d5-4620-a2be-c15c7616fbe8/data/coneseedtraits_v2.txt', sep='\t', header=T)

# cone length #

# put parent and offspring traits into one data frame 
off_cone_length<-cones$cone_length_mm
Family<-cones$fam
reg_data<-data.frame(Family, off_cone_length)
head(reg_data)

mean2 <- mat_cones %>%
  group_by(Family) %>%
  summarise(mean2 = mean(Le, na.rm = TRUE))

mean <- mean %>%
  left_join(mean2, by = "Family")

# run regression model
model_reg<-lm(mean~mean2, data=mean)
summary(model_reg)

# make plot 
CL<-ggplot(data=mean,
           aes(x=mean2, y=mean))+
  geom_point()+
  geom_abline(intercept = intercept, slope = slope,               
              color = 'red', size = 0.5)+
  theme_minimal()+
  labs(y='Mean Offspring Cone Length (mm)', x='Mean Maternal Cone Length (mm)')

# cone ratio #
off_cone_rat<-cones$cone_ratio
Family<-cones$fam
reg_data2<-data.frame(Family, off_cone_rat)
head(reg_data2)

# calculate maternal cone ratio
mat_cones$cone_ratio<-mat_cones$Wi/mat_cones$Le

cone_rat <- cones %>%
  group_by(fam) %>%
  summarise(mean_off_rat = mean(cone_ratio, na.rm = TRUE))

mat_cones$fam<-mat_cones$Family

cone_par<-mat_cones %>%
  group_by(fam) %>%
  summarise(mean_par_rat = mean(cone_ratio, na.rm = TRUE))

cone_rat <- cone_rat %>%
  left_join(cone_par, by = "fam")

model_reg_2<-lm(mean_off_rat~mean_par_rat, data=cone_rat)
summary(model_reg_2)

CR<-ggplot(data=cone_rat,
       aes(x=mean_par_rat, y=mean_off_rat))+
  geom_point()+
  geom_abline(intercept = intercept3, slope = slope3,               
              color = 'red', size = 0.5)+
  theme_minimal()+
  labs(y='Mean Offspring Cone Ratio', x='Mean Maternal Cone Ratio')

# seed number #
reg_data$seed_number<-cones$seed_number


seed_num <- reg_data %>%
  group_by(Family) %>%
  summarise(off_num = mean(seed_number, na.rm = TRUE))

seed_par<-mat_cones %>%
  group_by(Family) %>%
  summarise(par_num = mean(SN, na.rm = TRUE))

seed_num <- seed_num %>%
  left_join(seed_par, by = "Family")

model_reg_3<-lm(off_num~par_num, data=seed_num)
summary(model_reg_3)

SN<-ggplot(data=seed_num,
           aes(x=par_num, y=off_num))+
  geom_point()+
  geom_abline(intercept = intercept3, slope = slope3,               
              color = 'red', size = 0.5)+
  theme_minimal()+
  labs(y='Mean Offspring Seed Number', x='Mean Maternal Seed Number')

# Seed Viability #

reg_data$seed_viability<-cones$seed_viability

seed_via <- reg_data %>%
  group_by(Family) %>%
  summarise(off_via = mean(seed_viability, na.rm = TRUE))

mat_cones$seed_viability<-mat_cones$SV/100

par_via<-mat_cones %>%
  group_by(Family) %>%
  summarise(par_via = mean(seed_viability, na.rm = TRUE))

seed_via <- seed_via %>%
  left_join(par_via, by = "Family")

seed_via <- na.omit(seed_via)

model_reg_4<-lm(off_via~par_via, data=seed_via)
summary(model_reg_4)

CV<-ggplot(data=seed_via,
           aes(x=par_via, y=off_via))+
  geom_point()+
  geom_abline(intercept = intercept4, slope = slope4,               
              color = 'red', size = 0.5)+
  theme_minimal()+
  labs(y='Mean Offspring Seed Viability', x='Mean Maternal Seed Viability')

ggarrange(CL, CR, SN, SV,
          ncol = 2, nrow = 2,
          labels=c('A','B','C','D'),
          label.x= 0.85)

#### BLUPs ####

library(anyLib)
anyLib(c("apercu", "corpcor", "breedR", "RColorBrewer", "psych", "qqman"))

# Run BLUP models 
blups.trait_cn<-remlf90(fixed=cone_number~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_cl<-remlf90(fixed=cone_length_mm~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_stl<-remlf90(fixed=stalk_length_mm~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_scl<-remlf90(fixed=scale_length_mm~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_scn<-remlf90(fixed=scale_number~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_a<-remlf90(fixed=cone_angle~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_cr<-remlf90(fixed=cone_ratio~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_p<-remlf90(fixed=profile~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_o<-remlf90(fixed=openness~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_h<-remlf90(fixed=colour..hue.~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_v<-remlf90(fixed=value~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_ch<-remlf90(fixed=chroma~ 1 + pop+ block, random = ~ fam , data=cones)
blups.trait_sen<-remlf90(fixed=seed_number~ 1 + pop+ block, random = ~ fam , data=viability)
blups.trait_sev<-remlf90(fixed=seed_viability~ 1 + pop+ block, random = ~ fam , data=viability)

# Extract random effects 
cone_number<-blups.trait_cn$ranef 
cone_number<-as.data.frame(cone_number)
cone_number<-cone_number$fam.value

cone_length<-blups.trait_cl$ranef
cone_length<-as.data.frame(cone_length)
cone_length<-cone_length$fam.value

blup_cone<-data.frame(cone_number, cone_length)
  
  # Repeat for remaining traits # 

# Run genetic correlations 
blup_corr<-as.matrix(blup_values)
pcc<-rcorr(blup_corr, type=c("pearson"))
p<-pcc$P
r<-pcc$r
write.table(p, "/Users/jessshortman/Desktop/Dissertation/Final Data/p_blup_cones.txt", sep="\t")
write.table(r, "/Users/jessshortman/Desktop/Dissertation/Final Data/r_blup_cones.txt", sep="\t")

# Make heat map 

# Read in correlation results 
half_corr<-read.csv('BLUPs/Plot Data/gen_corr')

names(half_corr) <- c('X',"Cone Count", "Cone Length", "Stalk Length", 'Scale Length', 'Scale Count', 'Angle', 'Hue', 'Value', 'Chroma', "Profile ", 'Ratio', 'Openness', 'Seed Number', 'Seed Viability')
trait_order <- half_corr$X
trait_order

correlation_long <- melt(half_corr)
colnames(correlation_long) <- c("Trait1", "Trait2", "Correlation")

correlation_long$Trait1 <- factor(correlation_long$Trait1, levels = trait_order)
correlation_long$Trait2 <- factor(correlation_long$Trait2, levels = trait_order)

gen<-ggplot(correlation_long, aes(x = Trait1, y = Trait2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  labs(x = "Traits", y = "Traits")

# phenotypic correlations

## All populations ##
seed_corr<-read.csv('Seed_corr.csv')
cone_corr<-read.csv('Cones_corr.csv')

cone_dat <- na.omit(cone_corr)

cone_dat<-as.matrix(cone_corr)
seed_dat<-as.matrix(seed_corr)

pcc<-rcorr(cone_dat, type=c("pearson"))
p<-pcc$P
r<-pcc$r
write.table(p, "/Users/jessshortman/Desktop/Dissertation/Neat Correlations/p_cone_full.txt", sep="\t")
write.table(r, "/Users/jessshortman/Desktop/Dissertation/Neat Correlations/r_cone_full.txt", sep="\t")

pcc<-rcorr(seed_dat, type=c("pearson"))
p<-pcc$P
r<-pcc$r
write.table(p, "/Users/jessshortman/Desktop/Dissertation/Neat Correlations/p_seed_full.txt", sep="\t")
write.table(r, "/Users/jessshortman/Desktop/Dissertation/Neat Correlations/r_seed_full.txt", sep="\t")

# heat maps
phen_corr<-read.csv('BLUPs/Plot Data/phen_corr.csv')

names(phen_corr) <- c('X',"Cone Count", "Cone Length", "Stalk Length", 'Scale Length', 'Scale Count', 'Angle', 'Hue', 'Value', 'Chroma', "Profile ", 'Ratio', 'Openness', 'Seed Number', 'Seed Viability')
trait_order <- phen_corr$X
trait_order

correlation_long <- melt(phen_corr)
colnames(correlation_long) <- c("Trait1", "Trait2", "Correlation")

correlation_long$Trait1 <- factor(correlation_long$Trait1, levels = trait_order)
correlation_long$Trait2 <- factor(correlation_long$Trait2, levels = trait_order)

phen<-ggplot(correlation_long, aes(x = Trait1, y = Trait2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  labs(x = "Traits", y = "Traits")

ggarrange(gen, phen, nrow=1, ncol=2, labels=c('A','B'), label.x=0.8)

#### Environmental Variables ####

# Test for environmental correlations 
longitude<-cones$long
latitude<-cones$lat
temp<-cones$mean_july_temp
prec<-cones$annual_precip
alt<-cones$altitude

heicor<-data.frame(height,growth,longitude,temp,precip,alt)
heicor<-na.omit(heicor)
heicor<-cor(heicor)

# Environmental models with sequential removal of environmental predictors compared to a null model

# Cone Length # 
model_env_cone_length<-lmer(cone_length_mm ~long+height_mm+(1|pop) + (1|fam) + (1|block) + (1|Order), data=cones2, REML = FALSE)
summary(model_env_cone_length)

model_env_cone_length<-lmer(cone_length_mm ~(1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE)
summary(model_env_cone_length)

# Stalk Length #
model_env_stalk_length<-lmer(stalk_length_mm ~  long + lat + height_m+growth_rate+(1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE)

model_env_stalk_length<-lmer(stalk_length_mm ~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE)
summary(model_env_stalk_length)

# Scale Length #
model_env_scale_length<-lmer(scale_length_mm ~ annual_precip+ height_m+(1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE)
summary(model_env_scale_length)

model_env_scale_length<-lmer(scale_length_mm ~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE)
summary(model_env_scale_length)

# Cone Angle # 
model_env_cone_angle<-lmer(cone_angle~ height_m+ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML=FALSE)
summary(model_env_cone_angle)  

model_env_cone_angle<-lmer(cone_angle~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML=FALSE)
summary(model_env_cone_angle)  

# Hue #
model_env_hue<-lmer(colour..hue.~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE)
summary(model_env_hue) 

# Value #
model_env_value<-lmer(value~ height_m + growth_rate+(1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_value) 

model_env_value<-lmer(value~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_value) 

# Chroma # 
model_env_chroma<-lmer(chroma~ annual_precip+ height_m + growth_rate+ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_chroma) 

model_env_chroma<-lmer(chroma~ growth_rate+ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_chroma) 

# Profile#
model_env_profile <-lmer(profile~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_profile) 

# Openness #
model_env_openness <-lmer(openness~ height_m+(1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_openness) 

model_env_openness <-lmer(openness~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_openness) 

# Cone Ratio #
model_env_ratio <-lmer(cone_ratio ~ long + (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_ratio)

# Seed Number 
model_env_seed_number <-lmer(seed_number ~ height_m +(1|pop) + (1|fam) + (1|block) + (1|Order), data=viability, REML = FALSE )
summary(model_env_seed_number)

model_env_seed_number <-lmer(seed_number ~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=viability, REML = FALSE )
summary(model_env_seed_number)

# Scale Number 
model_env_scale_number <-lmer(scale_number ~ long+ height_m +(1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_scale_number)

model_env_scale_number <-lmer(scale_number ~ (1|pop) + (1|fam) + (1|block) + (1|Order), data=cones, REML = FALSE )
summary(model_env_scale_number)

# Viability #
model_env_via<-lmer(seed_viability~(1|pop) + (1|fam) + (1|block) + (1|Order), data=viability, REML = FALSE )
summary(model_env_via)

# Cone Number #
model_env_cone_number<-lmer(cone_number ~ height_m  + (1|fam)+ (1|block)+ (1|pop), data=subset_data, REML=FALSE)
summary(model_env_cone_number)

model_env_cone_number<-lmer(cone_number ~ (1|fam)+ (1|block)+ (1|pop), data=subset_data, REML=FALSE)
summary(model_env_cone_number)

## Plots ##

## Cone Length Longitude ##

summary_data <- cones %>%
  group_by(long) %>%
  summarise(
    y_mean = mean(cone_length_mm, na.rm = TRUE),                 
    y_se = sd(cone_length_mm, na.rm = TRUE) / sqrt(sum(!is.na(cone_length_mm))),  
    y_sd = sd(cone_length_mm, na.rm = TRUE)                       
  )

df_summary <- cones %>%
  group_by(pop) %>%
  summarise(long = unique(long))

summary_data <- summary_data %>%
  left_join(df_summary, by = "long")

le_long<-ggplot(data=summary_data,
             aes(x=long, y=y_mean))+
  geom_errorbar(aes(ymin = y_mean - y_se, ymax = y_mean + y_se, width=0.15))+
  geom_smooth(method = "lm", col = "red", size=0.5,se = FALSE, fullrange = TRUE)+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"))+
  scale_x_continuous(breaks = c(-5.5, -5, -4.5, -4, -3.5,-3))+
  labs(x='Longitude (°)', y= "Cone Length (mm)")+geom_text(data = summary_data,
                                                           aes(label = pop, y = y_mean + y_se +0.25 ),
                                                           vjust =3, hjust = -0.5, size = 3)

## Stalk length longitude ##

summary_data_stalk <- cones %>%
  group_by(long) %>%
  summarise(
    y_mean = mean(stalk_length_mm, na.rm = TRUE),                  
    y_se = sd(stalk_length_mm, na.rm = TRUE) / sqrt(sum(!is.na(stalk_length_mm))),  
    y_sd = sd(stalk_length_mm, na.rm = TRUE)                      
  )

df_summary <- cones %>%
  group_by(pop) %>%
  summarise(long = unique(long))

summary_data_stalk <- summary_data_stalk %>%
  left_join(df_summary, by = "long")

stlong<-ggplot(data=summary_data_stalk,
               aes(x=long, y=y_mean))+
  geom_errorbar(aes(ymin = y_mean - y_se, ymax = y_mean + y_se, width=0.1))+
  geom_smooth(method = "lm", col = "red", size=0.5,se = FALSE, fullrange = TRUE)+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"))+
  scale_x_continuous(breaks = c(-5.5, -5, -4.5, -4, -3.5,-3))+
  labs(x='Longitude (°)', y= "Stalk Length (mm)")+
  geom_text(aes(label = pop), hjust = -0.2, vjust = 2, size = 3)+
  scale_y_continuous(limits = c(8, 10.5))

ggarrange(le_long, st_long, ncol=2, nrow=1, labels=c('A','B'),label.x= 0.9, label.y=0.98)


