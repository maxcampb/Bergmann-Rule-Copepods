
## Packages used

library(tidyverse)
library(viridis)
library(lme4)
library(lmerTest)
library(sp)
library(gstat)
library(patchwork)
library(optimx)


# Exploratory data analysis -----------------------------------------------

# Load the processed data (available)
load("Data/CopeData.rda")

# Look at spearman correlations and distributions
ggpairs_plot <- GGally::ggpairs(CopeData, columns = c("MLC", "SST", "Sqrt_chl", "Asin_omni",
                                                      "Oxygen", "Lat_abs", "SST_mon"), progress = FALSE, 
                                upper = list(continuous = GGally::wrap("cor", method = "spearman", stars = FALSE))
)

### Look at the transformations ###

with(CopeData, {
  par(mfrow = c(2,4))
  
  # Mean length
  hist(MLC, xlab = "Mean Length of Copepods", prob = TRUE,
       main = "", cex.lab = 0.9, col = "grey",
       border = "black") # Heavily right skewed
  pu <- par("usr")
  text(pu[2], pu[4], "(A)", adj = c(1.5,1), cex = 1.5)
  
  # Mean length transformed
  hist(log(MLC), xlab = "Log(Mean Length of Copepods)", prob = TRUE,
       main = "", cex.lab = 0.9, col = "grey",
       border = "black") # much better
  pu <- par("usr")
  text(pu[2], pu[4], "(B)", adj = c(1.5,1), cex = 1.5)
  
  # Chl-a
  hist(Chl, xlab = "Chl-a", main = "", prob = TRUE, 
       cex.lab = 0.9, col = "grey", 
       border = "black") # Heavily right skewed
  pu <- par("usr")
  text(pu[2], pu[4], "(C)", adj = c(1.5,1), cex = 1.5)
  
  # Chl-a transformed
  hist(sqrt(Chl), xlab = "sqrt(Chl-a)", main = "", prob = TRUE, 
       cex.lab = 0.9, col = "grey", 
       border = "black") # much better
  pu <- par("usr")
  text(pu[2], pu[4], "(D)", adj = c(1.5,1), cex = 1.5)
  
  # Omni_prop
  hist(Omni_prop, xlab = "Omnivore proportion", main = "", prob = TRUE,
       cex.lab = 0.9, col = "grey", 
       border = "black") # Heavily right skewed
  pu <- par("usr")
  text(pu[2], pu[4], "(G)", adj = c(1.5,1), cex = 1.5)
  
  # Omni_prop transformed
  hist(Asin_omni, xlab = "arcsine(sqrt(Omnivore proportion))", main = "", prob = TRUE, 
       cex.lab = 0.9, col = "grey", 
       border = "black") # Heavily right skewed
  pu <- par("usr")
  text(pu[2], pu[4], "(H)", adj = c(1.5,1), cex = 1.5)
  
})


### Comparing size by diet ###

size_diet_df$Diet[size_diet_df$Diet == "herb"] <- "omni"
size_diet_df <- size_diet_df %>% filter(Diet %in% c("omni", "carn"))

# Histograms for each group
hist(size_diet_df$Size_resolved[size_diet_df$Diet == "omni"])
hist(size_diet_df$Size_resolved[size_diet_df$Diet == "carn"])

# Density plot
plot(density(size_diet_df$Size_resolved[size_diet_df$Diet == "omni"], bw = 0.37), 
     col = "darkgreen", main = "", xlab = "Mean Length of Copepod (mm)", 
     lwd = 1.5, cex.axis=1.2, cex.lab = 1.3)
lines(density(size_diet_df$Size_resolved[size_diet_df$Diet == "carn"], bw = 0.37),
      lwd = 1.5, col = "darkred")
legend(7, 0.3, legend=c("Omnivores", "Carnivores"),
       col=c("darkgreen", "darkred"), lty = c(1,1), lwd = c(1.5, 1.5))


# Step 1: WO18 data models to compare Oxygen, SST and Latitude ------------

# Filter out missing data
CopeData_WO18 <- CopeData %>% filter(rowSums(is.na(CopeData)) == 0)

# Run competing models
SST_mod <- glmer(MLC ~ SST_mon + Sqrt_chl + Asin_omni +
                   (1|Survey) + (1 + TowDays|Survey:Tow_No) + ( 1 |Longhurst),
                 data = CopeData_WO18, family = Gamma(link = log),
                 control = glmerControl(optimizer = "bobyqa"))

Oxy_mod <- glmer(MLC ~ Oxygen + Sqrt_chl + Asin_omni +
                   (1|Survey) + (1 + TowDays|Survey:Tow_No) + ( 1 |Longhurst),
                 data = CopeData_WO18, family = Gamma(link = log),
                 control = glmerControl(optimizer = "bobyqa"))

Lat_mod <- glmer(MLC ~ Lat_abs + Sqrt_chl + Asin_omni + 
                   (1|Survey) + (1 + TowDays|Survey:Tow_No) + ( 1 |Longhurst),
                 data = CopeData_WO18, family = Gamma(link = log),
                 control = glmerControl(optimizer = "bobyqa"))

# Compare models based on BIC
anova(SST_mod, Oxy_mod, Lat_mod)

# Get pseudo-R2 value
MuMIn::r.squaredGLMM(SST_mod)
MuMIn::r.squaredGLMM(Lat_mod)
MuMIn::r.squaredGLMM(Oxy_mod)


# Step 2: Higher resolution data models -----------------------------------

# Full model
hr_mod <- glmer(MLC ~ SST + Sqrt_chl + Asin_omni + 
                  (1|Survey) + (1 + TowDays|Survey:Tow_No) + ( 1 |Longhurst),
                data = CopeData, family = Gamma(link = log),
                control = glmerControl(optimizer = "optimx", optCtrl=list(method = "nlminb")))

cat("All Fixed effects\n")
MuMIn::r.squaredGLMM(hr_mod)


### Backward removals ###

## Fixed effect removals
hr_mod1 <- update(hr_mod, . ~ . -SST)
cat("No SST\n")
MuMIn::r.squaredGLMM(hr_mod1)

hr_mod2 <- update(hr_mod, . ~ . -Sqrt_chl)
cat("No Chl-a\n")
MuMIn::r.squaredGLMM(hr_mod2)

hr_mod3 <- update(hr_mod, . ~ . -Asin_omni)
cat("No Omnivore proportion\n")
MuMIn::r.squaredGLMM(hr_mod3)

## Random effect removals
hr_mod4 <- update(hr_mod, . ~ . -(1|Survey))
cat("No Survey\n")
MuMIn::r.squaredGLMM(hr_mod5)

hr_mod5 <- update(hr_mod, . ~ . -(1 + TowDays|Survey:Tow_No) + (1|Survey:Tow_No))
cat("No Slope for TowDays\n")
MuMIn::r.squaredGLMM(hr_mod6)

hr_mod6 <- update(hr_mod, . ~ . -(1 + TowDays|Survey:Tow_No))
cat("No Tow Number\n")
MuMIn::r.squaredGLMM(hr_mod7)

hr_mod7 <- update(hr_mod, . ~ . -(1 | Longhurst))
cat("No Longhurst Province\n")
MuMIn::r.squaredGLMM(hr_mod8)


# Model selection based on BIC
anova(hr_mod, hr_mod1, hr_mod2, hr_mod3, hr_mod4, hr_mod5, hr_mod6, hr_mod7)


# Checking assumptions ----------------------------------------------------

### Normality ###

## QQplot
qqnorm(resid(hr_mod, type = "deviance"), main = "")
qqline(resid(hr_mod, type = "deviance"))


### Heteroskedasticity ###

## Residual plot
ggplot(CopeData) + aes(fitted(hr_mod), resid(hr_mod, type = "pearson")) + 
  geom_hex(bins = 80) + theme_bw() + 
  scale_fill_viridis(begin = 0, end = .9, option = "C", limit = range(c(1,1200))) +
  theme(plot.title = element_text(hjust = 0.5), ) + 
  labs(fill = "Frequency") + ylab("Residual") + xlab("Fitted Value")

## Residuals by predictors

# SST vs residuals
ggplot(CopeData) + aes(SST, resid(hr_mod, type = "pearson")) + 
  geom_hex(bins = 80) + theme_bw() + 
  scale_fill_viridis(begin = 0, end = 1, option = "C", limit = range(c(1,400))) + 
  labs(title = "Hexplot of residuals for SST") + 
  theme(plot.title = element_text(hjust = 0.5))

# Chl vs residuals
ggplot(CopeData) + aes(sqrt(Chl), resid(hr_mod, type = "pearson")) + 
  geom_hex(bins = 80) + theme_bw() + 
  scale_fill_viridis(begin = 0, end = 1, option = "C", limit = range(c(1,1000))) +
  labs(title = "Hexplot of residuals for chl") + 
  theme(plot.title = element_text(hjust = 0.5))

# Omni_prop vs residuals
ggplot(CopeData) + aes(asin_omni, resid(hr_mod, type = "pearson")) + 
  geom_hex(bins = 80) + theme_bw() + 
  scale_fill_viridis(begin = 0, end = 1, option = "C", limit = range(c(1,1000))) +
  labs(title = "Hexplot of residuals for omnivore proportion") +
  theme(plot.title = element_text(hjust = 0.5))


### Temporal autocorrelation ###

# Take subset of the data
tows <- with(CopeData, unique(Tow_No)[1:60])
CopeDataSubset <- subset(CopeData, Tow_No %in% unique(Tow_No)[1:60])

# Plot residuals by tow over time WITHOUT slope for tow days
ggplot(data = CopeDataSubset, aes(Date, resid(hr_mod6)[CopeData$Tow_No %in% tows])) + 
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  theme_bw() + labs(title = "Without slope for tow_days") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~ Tow_No, ncol = 10, scales = "free_x")

# Plot residuals by tow over time WITH slope for tow days
ggplot(data = CopeDataSubset, aes(Date, resid(hr_mod)[CopeData$Tow_No %in% tows])) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw() + 
  labs(title = "With slope for tow_days") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~ Tow_No, ncol = 10, scales = "free_x")


### Spatial autocorrelation ###

## Spatial variograms

# Data format for variograms
autocorData <- data.frame(CopeData$Lon, CopeData$Lat, 
                          resids=resid(hr_mod)) # Extract residuals
names(autocorData) <- c('Lon', 'Lat', 'resids')
coordinates(autocorData) <- c('Lon', 'Lat') # Convert lon and lat to coordinates

# Look at variogram in all directions together
varmodel <- variogram(resids~1,data=autocorData, cutoff = 10)
plot(varmodel, main = "Variogram of all directions combined")

# Look at variogram in all directions separately
varmodel<-variogram(resids~1, data=autocorData, alpha=c(0,45,90,135), cutoff = 10)
plot(varmodel, main = "Variogram of all directions separately") 


# Check for spatial auto-correllation
autocorData <- data.frame(Lon = CopeData$Lon, 
                          Lat = CopeData$Lat, 
                          resids = resid(hr_mod)) %>%  # Extract residuals of GLMM
  within({
    signres <- sign(resids)
  })


## Plot the residuals in space
world <- map_data("world")

# Change longitude so it matches up with the world map
autocorData$Lon[autocorData$Lon < (-170)] <- autocorData$Lon[autocorData$Lon < (-170)] + 360

# Bubble plot WITH random effects
ggplot(data = autocorData, aes(x = Lon, y = Lat)) + 
  geom_map(data=world, map=world, aes(x = long, y = lat, map_id=region),
           color="white", fill="gray94", size=0.08) + 
  geom_point(aes(size = abs(resids), color = sign(resids)), shape = 1,
             alpha = 0.4) + 
  scale_size_continuous(range=c(.1,4)) + 
  scale_colour_gradient(low = "springgreen3", high = "magenta3") + 
  ylab(NULL) + xlab(NULL) + 
  annotate("text", x = -190, y = 90, label = "(b)", size = 9) +
  guides(colour = "none", size = guide_legend(title = "Magnitude"))


# Model without random effects
hr_mod8 <- glm(MLC ~ SST + Sqrt_chl + Asin_omni,  
               data = CopeData, family = Gamma(link = log))

# Use residuals from the GLM
autocorData$resids <- resid(hr_mod8) # Extract residuals of GLM
autocorData$signres<- sign(autocorData$resids)

# Bubble plot WITHOUT random effects
ggplot(data = autocorData, aes(x = Lon, y = Lat)) + 
  geom_map(data=world, map=world, aes(x = long, y = lat, map_id=region),
           color="white", fill="gray94", size=0.08) + 
  geom_point(aes(size = abs(resids), color = sign(resids)), shape = 1,
             alpha = 0.4)  + 
  scale_size_continuous(range=c(.1,4)) + 
  scale_colour_gradient(low = "springgreen3", high = "magenta3") + 
  ylab(NULL) + xlab(NULL) + 
  annotate("text", x = -190, y = 90, label = "(a)", size = 9) +
  guides(colour = "none", size = guide_legend(title = "Magnitude"))


# Effect sizes and calculations -------------------------------------------

### Proportion of variance explained by random effects ###

# conditional - marginal R^2
RE_var <- 0.4979803 - 0.09215131

# From the summary objects
0.0081754/(0.0438111 + 0.0952570 + 0.0001882 + 0.0081754) * RE_var # survey/all
(0.0438111 + 0.0952570)/(0.0438111 + 0.0952570 + 0.0001882 + 0.0081754) * RE_var # tow effects /all
0.0001882/(0.0438111 + 0.0952570 + 0.0001882 + 0.0081754) * RE_var # longhurst/all


### Interpretations of fixed effects ###

# Extract Coefficients
coef_tab <- coefficients(summary(hr_mod))

# Convert to the natural scale
exp(coef_tab[, "Estimate"])


### Omnivore Proportion effect sizes ###

# Setup a prediction dataframe
df <- head(CopeData, 2)
df <- df %>% select(MLC, SST,  Chl,  Asin_omni) %>% 
  within({
    SST <- mean(CopeData$SST)
    Sqrt_chl <- mean(CopeData$Sqrt_chl)
    Omni_prop <- c(1,0)
    Asin_omni <- asin(sqrt(Omni_prop))
    Tow_No <- NULL
    Longhurst <- NULL
  })

# Make predictions based on our model
y_hat <- predict(hr_mod, df, re.form = NA, type = "response")

# Percentage change in length
1 - y_hat[1] / y_hat[2]

# Absolute change
y_hat[2] - y_hat[1]

convert_to_weight <- function(x){
  0.03493 * x ^ 2.9878
}

# Percentage change in mass
ww_y_hat <- convert_to_weight(y_hat)
1 - ww_y_hat[1] / ww_y_hat[2]


### SST effect sizes ###

df <- df %>% select(MLC, SST,  Chl,  Asin_omni) %>% 
  within({
    SST <- c(min(CopeData$SST), max(CopeData$SST))
    Sqrt_chl <- mean(CopeData$Sqrt_chl)
    Asin_omni <- mean(CopeData$Asin_omni)
    Tow_No <- NULL
    Longhurst <- NULL
  })

# Make predictions based on our model
y_hat <- predict(hr_mod, df, re.form = NA, type = "response")

summary(hr_mod)

# Calculate the proportional change
1-y_hat[2]/y_hat[1]

# Absolute change
y_hat[2] - y_hat[1]

# length per degree (linear approximation)
(y_hat[1] - y_hat[2])/ diff(df$SST)

# Percentage change in Mass
ww_y_hat <- convert_to_weight(y_hat)
1 - ww_y_hat[2] / ww_y_hat[1]

# Percentage change in Mass / degree
(1 - ww_y_hat[2] / ww_y_hat[1])/(29.966+1.781)


### CHL-a effects sizes ###

df <- df %>% select(MLC, SST,  Chl,  Asin_omni) %>% 
  within({
    SST <- mean(CopeData$SST)
    Chl <- c(0.018596824,	9.410631) # removed 2 influential obs one on each end
    Sqrt_chl <- sqrt(Chl)
    Asin_omni <- mean(CopeData$Asin_omni)
    Tow_No <- NULL
    Longhurst <- NULL
  })

# Make predictions based on our model
y_hat <- predict(hr_mod, df, re.form = NA, type = "response")

# Calculate the proportional change
1- y_hat[2]/y_hat[1]

# Absolute change
y_hat[2] - y_hat[1]

# length per degree
(y_hat[1] - y_hat[2])/ diff(df$Sqrt_chl)

# Percentage change in Mass
ww_y_hat <- convert_to_weight(y_hat)
1 - ww_y_hat[2] / ww_y_hat[1]


# Chl predictions for climate change --------------------------------------

# Chlorophyll changes
old_NPP <- 52.1 # As per Bopp et al. 2013

# As per Bopp et al. 2013 (percent changes)
RCP8.5_per <- -8.6 
RCP2.6_per <- -2.0 

# Calculate future estimates of NPP
new_NPP_RCP8.5 <- old_NPP + old_NPP * RCP8.5_per/100
new_NPP_RCP2.6 <- old_NPP + old_NPP * RCP2.6_per/100

# Conversion to Chl from Maranon et al. 2014
convert_to_chl <- function(x) { 10^((log10(x) - 1.58)/1.29) }

# Absolute changes
abs_chl_change_RCP8.5<- convert_to_chl(new_NPP_RCP8.5) - convert_to_chl(old_NPP)
abs_chl_change_RCP2.6<- convert_to_chl(new_NPP_RCP2.6) - convert_to_chl(old_NPP)


### RCP8.5 predictions ###

newdat <- with(CopeData, data.frame( Sqrt_chl = c(sqrt(convert_to_chl(old_NPP)),
                                                  sqrt(convert_to_chl(new_NPP_RCP8.5))),
                                     SST  = mean(SST),
                                     Asin_omni = mean(Asin_omni)))

# Make predictions based on our model
y_hat <- predict(hr_mod, newdata = newdat, re.form = NA, type = "response")

# Percentage change in Mass
ww <- convert_to_weight(y_hat)
ww[2] / ww[1]


### RCP2.6 predictions ###

newdat <- with(CopeData, data.frame( Sqrt_chl = c(sqrt(convert_to_chl(old_NPP)),
                                                  sqrt(convert_to_chl(new_NPP_RCP2.6))),
                                     SST  = mean(SST),
                                     Asin_omni = mean(Asin_omni)))

# Make predictions based on our model
y_hat <- predict(hr_mod, newdata = newdat, re.form = NA, type = "response")

# Percentage change in Mass
ww <- convert_to_weight(y_hat)
ww[2] / ww[1]


# Combined changes
1-1*1.012*0.927 #RCP8.5
1-1*1.003*0.984 #RCP2.6


# Plots -------------------------------------------------------------------

# Ben Bolker's function for computing confidence intervals based on fixed effects only
easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
  # baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  # fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
  beta <- fixef(model) # fixed-effects coefficients
  V <- vcov(model)     # variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) # std errors of predictions
  # inverse-link function
  linkinv <- family(model)$linkinv
  # construct 95% Normal CIs on the link scale and
  # transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low=pred0-crit*pred.se,
                conf.high=pred0+crit*pred.se))
}


### SST plot ###

# Make data that spans the SST range
newdat <- with(CopeData, data.frame( SST = seq(min(SST), max(SST), length.out = 100),
                                     Sqrt_chl  = mean(Sqrt_chl),
                                     Asin_omni = mean(Asin_omni)))

# Make predictions
y_pred <- predict(hr_mod,re.form=NA,newdata=newdat,type="response")
# Estimate the confidence intervals using Ben Bolker's CI 
# function (based on uncertainty in fixed effects only)
conf_int <- easyPredCI(hr_mod,newdata = newdat)

# Bind data and predictions
newdat <- cbind(newdat, y_pred, conf_int)

# Specify parameters of plotting
annotate_sz <- 14
point_sz <- 2.4
line_sz <- 2
error_sz <- 1.3

# Make effect plot
ggplot(data = CopeData, mapping = aes(x = SST, y = fitted(hr_mod))) + 
  geom_ribbon(data = newdat, mapping = aes( y = y_pred,  ymin = conf.low, ymax =conf.high,
                                            x = SST),  alpha = .9, fill = "grey80") +
  geom_point(alpha = 0.2, size = point_sz) +
  geom_line(data = newdat, aes(x = SST, y = y_pred), col = "dodgerblue", size = line_sz) +
  geom_line(data = newdat, aes(x = SST, y = conf.low), col = "dodgerblue", 
            size = error_sz, lty = "dashed") +
  geom_line(data = newdat, aes(x = SST, y = conf.high), col = "dodgerblue", 
            size = error_sz, lty = "dashed") +
  theme_bw() + annotate("text", x = 0, y = 6, label = "(d)", size = annotate_sz) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), 
        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank()) + xlab("SST (°C)") + 
  ylab("Mean length of copepod (mm)") +
  scale_x_continuous(breaks = seq(-0, 30, by = 5)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) 


### Chl-a ###

# Make data that spans the Chl-a range
newdat <- with(CopeData, data.frame( Sqrt_chl = seq(min(Sqrt_chl), 
                                                    max(Sqrt_chl), length.out = 100),
                                     SST  = mean(SST),
                                     Asin_omni = mean(Asin_omni)))

# Make predictions
y_pred <- predict(hr_mod,re.form=NA,newdata=newdat,type="response")
# Estimate the confidence intervals using Ben Bolker's CI 
# function (based on uncertainty in fixed effects only)
conf_int <- easyPredCI(hr_mod,newdata = newdat)

# Bind data and predictions
newdat <- cbind(newdat, y_pred, conf_int)

# Make effect plot
ggplot(data = CopeData, mapping = aes(x = Chl, y = fitted(hr_mod))) + 
  geom_ribbon(data = newdat, mapping = aes( y = y_pred,  ymin = conf.low,
                                            ymax=conf.high, x = Sqrt_chl^2),  
              alpha = .9, fill = "grey80") +
  geom_point(alpha = 0.2, size = point_sz) +
  geom_line(data = newdat, aes(x = Sqrt_chl^2, y = y_pred), col = "dodgerblue", 
            size = line_sz) +
  geom_line(data = newdat, aes(x = Sqrt_chl^2, y = conf.low), col = "dodgerblue", 
            size = error_sz, lty = "dashed") +
  geom_line(data = newdat, aes(x = Sqrt_chl^2, y = conf.high), col = "dodgerblue", 
            size = error_sz, lty = "dashed") +
  theme_bw() + annotate("text", x = 0.54, y = 6, label = "(e)", size = annotate_sz) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), 
        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank()) + xlab("Chl-a") + 
  ylab("Mean length of copepod (mm)") +
  scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) 


### Omnivore proportion ###

# Make data that spans the range of proportion of omnivores
newdat <- with(CopeData, data.frame( Asin_omni = seq(min(Asin_omni), 
                                                     max(Asin_omni), length.out = 100),
                                     SST  = mean(SST),
                                     Sqrt_chl = mean(Sqrt_chl)))

# Make predictions
y_pred <- predict(hr_mod,re.form=NA,newdata=newdat,type="response")
# Estimate the confidence intervals using Ben Bolker's CI 
# function (based on uncertainty in fixed effects only)
conf_int <- easyPredCI(hr_mod,newdata = newdat)

# Bind data and predictions
newdat <- cbind(newdat, y_pred, conf_int)

# Inverse of the arcsine sqrt transformation
inv_arcsine <- function(x) {sin(x)^2}

# Make effect plot
ggplot(data = CopeData, mapping = aes(x = inv_arcsine(Asin_omni), y = fitted(hr_mod))) + 
  geom_ribbon(data = newdat, mapping = aes( y = y_pred,  ymin = conf.low,
                                            ymax=conf.high, x = inv_arcsine(Asin_omni)),  
              alpha = .9, fill = "grey80") +
  geom_point(alpha = 0.2, size = point_sz) +
  geom_line(data = newdat, aes(x = inv_arcsine(Asin_omni), y = y_pred), 
            col = "dodgerblue", size = line_sz) +
  geom_line(data = newdat, aes(x = inv_arcsine(Asin_omni), y = conf.low), 
            col = "dodgerblue", size = error_sz, lty = "dashed") +
  geom_line(data = newdat, aes(x = inv_arcsine(Asin_omni), y = conf.high), 
            col = "dodgerblue", size = error_sz, lty = "dashed") +
  theme_bw() + annotate("text", x = 0.05, y = 6, label = "(g)", size = annotate_sz) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), 
        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank()) + xlab("Proportion of Omnivores") + 
  ylab("Mean length of copepod (mm)") +
  scale_x_continuous(breaks = seq(0, 1, by = .2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) 


### Survey ###

# Extract random effects
REs <- ranef(hr_mod, condVar = TRUE) 

# Extract variance
qq <- attr(ranef(hr_mod, condVar = TRUE)$Survey, "postVar")

# Extract intercepts
rand.interc <- REs$Survey

# Create a dataframe for plotting
df <- data.frame(Intercepts = REs$Survey[,1],
                 sd.interc = 2*sqrt(qq[,,1:length(qq)]),
                 lev.names = rownames(rand.interc)) %>% 
  arrange(Intercepts) %>% 
  within({  # Reorder levels
    lev.names <- factor(as.character(lev.names),
                        as.character(lev.names))
  })

# Make random effect plot
ggplot(df, aes(lev.names, Intercepts)) +
  geom_hline(yintercept=0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin=Intercepts-sd.interc,
                    ymax=Intercepts+sd.interc), 
                width = 0,color="black") +
  geom_point(aes(color = lev.names), size = 9) +
  guides(size = 16, shape = "none", 
         color = guide_legend(reverse = TRUE , override.aes = list(size=8))) + 
  theme_bw() + annotate("text", y = -.3, x = 4, label = "(b)", size = annotate_sz) +
  coord_flip() + ylab("Intercept") + xlab("") + 
  scale_color_manual(values=c("green4", "#E69F00", "purple", "royalblue1")) +
  scale_y_continuous(breaks = seq(-.2, .4, by = .2)) +
  theme(axis.text.x=element_text(size=10), axis.title.x=element_text(size=13),
        axis.text.y=element_blank(), panel.grid.minor=element_blank(), 
        axis.ticks.y = element_blank(), legend.position=c(.7,.2),
        legend.title=element_blank(), legend.text = element_text(size=26))


### Tow slope and intercept ###

# Create a dataframe for plotting
towDf <- as.data.frame(REs$`Survey:Tow_No`) %>% setNames(c("Intercept", "Slope"))

# Specify frequency breaks
my_breaks <- c(1,5,30,175,1000)

# Make random effect hexplot
ggplot(towDf, aes(Intercept, Slope)) +
  geom_hex(bins = 30) + theme_bw() + 
  scale_fill_viridis(begin = 0, end = 1, option = "D", limit = range(c(1,1100)), 
                     trans = "log", breaks = my_breaks, labels = my_breaks)  +
  ylab("Slope") + xlab("Intercept") + labs(fill = "Frequency") +
  guides(size = "none", shape = "none", 
         color = guide_legend(reverse = TRUE , override.aes = list(size=3))) + 
  annotate("text", y = .2, x = -.33, label = "(c)", size = annotate_sz) +
  scale_y_continuous(breaks = seq(-0.2, .4, by = .1), limits = c(-.25, .25 )) +
  scale_x_continuous(breaks = seq(-0.3, .6, by = .3))


### Longhurst Province ###

qq <- attr(ranef(hr_mod, condVar = TRUE)$Longhurst, "postVar")

# Extract intercepts
rand.interc <- REs$Longhurst

# Make a dataframe for plotting
df <- data.frame(Intercepts = REs$Longhurst[,1],
                 sd.interc = 2*sqrt(qq[,,1:length(qq)]),
                 lev.names = factor(rownames(rand.interc))) %>% 
  arrange(Intercepts) %>% 
  within({  # Reorder levels
    lev.names <- factor(as.character(lev.names),
                        as.character(lev.names))
  })

ggplot(df, aes(lev.names, Intercepts))+ 
  geom_hline(yintercept=0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin=Intercepts-sd.interc, 
                    ymax=Intercepts+sd.interc),
                width = 0,color="black") +
  geom_point(color = "black", size = 4.5) +
  guides(size=FALSE,shape=FALSE) + theme_bw() +
  theme(axis.text.x=element_text(size=10), 
        axis.title.x=element_text(size=13),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  coord_flip() + 
  ylab("Intercept") + xlab("") +
  annotate("text", y = -.15, x = 28, label = "(a)", size = annotate_sz)

### CPR map ###

# Get map off world
world <- map_data("world")

# Change longitude to fit map
CopeData$Lon[CopeData$Lon < (-170)] <- CopeData$Lon[CopeData$Lon < (-170)] + 360

# Plot the surveys on a world map
ggplot(data = CopeData, aes(x = Lon, y = Lat))  + coord_quickmap() +
  theme_bw() + scale_color_manual(values=c("green4", "#E69F00", "purple", "royalblue1")) +
  geom_map(data=world, map=world, aes(x = long, y = lat, map_id=region), 
           color="white", fill="gray80", size=0.08) + 
  geom_point(data = CopeData, aes(color = Survey), size = .3) +
  theme( legend.title = element_text(size=16, face="bold"), 
         legend.text = element_text(size = 12), 
         legend.background = element_rect(color="black", size=.3),
         legend.position=c(.13, .4)) + 
  guides(color = guide_legend(override.aes = list(size=3))) + 
  ylab(NULL) + xlab(NULL) + xlim(-170, 190) + ylim(-85,85) +
  scale_x_continuous(breaks = seq(-180, 180, by = 30), 
                     labels = paste0(seq(-180, 180, by = 30), "º")) +
  scale_y_continuous(breaks = seq(-90, 90, by = 30), 
                     labels = paste0(seq(-90, 90, by = 30), "º")) 
