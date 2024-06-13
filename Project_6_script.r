library(tidyverse)
library(lme4)
library(pwr)

df <- read.csv("C:\\Users\\Damian\\Documents\\SZKOŁA\\SEMESTR 10\\nuno\\Fluge_etal_2019_data.csv")

### PROPORTION IN RESPONSE AFTER 24 MONTHS #################################################

# We group the data by treatment_group so 0 is for placebo and 1 is rituximab
treatment_groups <- df %>%
  group_by(treatment_group)

# Then we calculate the proportion of responses based on already existing column response
prop_responded_by_group <- treatment_groups %>%
  summarise(
    prop_responded = mean(response == 1)
  )

prop_responded_by_group

prop_responded_df <- prop_responded_by_group %>%
  ungroup()

treatment_0_prop <- prop_responded_df$prop_responded[1]
treatment_1_prop <- prop_responded_df$prop_responded[2]

n_patients_per_group <- df %>%
  count(treatment_group)

treatment_0_count <- n_patients_per_group$n[1]

treatment_1_count <- n_patients_per_group$n[2]

# We use Wilson Score Interval to caluclate the uncertainty for each proportion
treatment_0_CI <- BinomCI(treatment_0_prop*treatment_0_count, treatment_0_count, conf.level = 0.95, sides = "two.sided",
                          method = "wilson")
treatment_0_CI

treatment_1_CI <- BinomCI(treatment_1_prop*treatment_1_count, treatment_1_count, conf.level = 0.95, sides = "two.sided",
                          method = "wilson")
treatment_1_CI

# Now we create a new column which represents response by the difference between FSS_24m and FSS_0m
# If the difference is negative then the response is 1
# If the difference is positive then the response is 0
df <- df %>%
  mutate(
    response_FSS = ifelse(FSS_24m - FSS_0m < 0, 1, 0)  
  )

# After creating such variable we continue as we did previously and finally use Wilson Score Interval to calculate the uncertainty
treatment_groups <- df %>%
  group_by(treatment_group)

prop_responded_by_group <- treatment_groups %>%
  summarise(
    prop_responded = mean(response_FSS == 1)
  )

prop_responded_by_group

prop_responded_df <- prop_responded_by_group %>%
  ungroup()

treatment_0_prop <- prop_responded_df$prop_responded[1]
treatment_1_prop <- prop_responded_df$prop_responded[2]

n_patients_per_group <- df %>%
  count(treatment_group)

print(n_patients_per_group)

treatment_0_count <- n_patients_per_group$n[1]

treatment_1_count <- n_patients_per_group$n[2]

treatment_0_CI <- BinomCI(treatment_0_prop*treatment_0_count, treatment_0_count, conf.level = 0.95, sides = "two.sided",
                          method = "wilson")
treatment_0_CI

treatment_1_CI <- BinomCI(treatment_1_prop*treatment_1_count, treatment_1_count, conf.level = 0.95, sides = "two.sided",
                          method = "wilson")
treatment_1_CI

### MODEL FOR FSS DIFFERENCE AFTER 24 MONTHS ###############################################

  #data cleaning and preprocessing
  #we create variables representing nominal and percentage differences between
  #FSS at 2 timestamps
  # NAs in "Infeksjon" variable are not imputed, we treat this as third category
  df <- df %>% mutate(
    FSS_diff_24m = FSS_24m - FSS_0m,
    FSS_diff_24m_rel = FSS_diff_24m / FSS_0m,
    Center = ď.żCenter,
    Infection = ifelse(is.na(Infeksjon), 2, Infeksjon) %>%  as.factor(),
  )
  
  #there are 2 missing variables in 2 columns in total - replacing with a mean
  df[is.na(df$Stepsmean_0m), ]$Stepsmean_0m <- mean(df$Stepsmean_0m, na.rm=T)
  df[is.na(df$Height), ]$Height <- mean(df$Height, na.rm=T)
  
  # a) BASE MODEL ------------------------------------------------------------------------
  
    #fitting the model. After experiments we noticed that FSS_0m is potentially important,
    # but its effect is not linear. Patients with highest severity (>60) seem to perform best
    model <- lmer(FSS_diff_24m ~ Age_inclusion + Sex + treatment_group + Infection + 
                  factor(SykVarig) +
                  Weight_0m + Height + log(Stepsmean_0m) + I(FSS_0m > 60) + (1 | Center),
                data = df)
    summary(model)
    
    car::Anova(model, type = 3) 
    #only I(FSS_0m > 60) significant
    #Age_inclusion has the second lowest p-value, patients perform worse as they age
    step(model)
    #step AIC confirms that model with only 1 variable FSS_0m is best
  
    residuals(model) %>% hist
    #the residuals are skewed, we try to apply box-cox transformation to improve that
    
    
  # b) BOX-COX TRANSFORMATION -------------------------------------------------------------

    # fitting a model without random effect to estimate optimal lambda
    # target variable shifted up to ensure positivity
    model_bc <- lm(FSS_diff_24m + 48 ~ Age_inclusion + Sex + treatment_group + Infection + 
                  factor(SykVarig) +
                  Weight_0m + Height + log(Stepsmean_0m) + I(FSS_0m > 60), 
                data = df_model)
    bc <- MASS::boxcox(model_bc, lambda = seq(0,5,0.1))
    
    #lambda = 4 looks optimal
    lambda = 4

    #creating new variable -> normalizing
    df_model <- df_model %>% 
      mutate(
        FSS_diff_24m_bc = ((FSS_diff_24m + 48)^lambda-1) / (lambda*100000)
      )
    df_model$FSS_diff_24m_bc %>% range

    model_bc_final <- lmer(FSS_diff_24m_bc ~ Age_inclusion + Sex + treatment_group + Infection + 
                     factor(SykVarig) +
                     Weight_0m + Height + log(Stepsmean_0m) + I(FSS_0m > 60) + (1 | Center),
                   data = df_model)
    summary(model_bc_final)

    residuals(model_bc_final) %>% hist
    #residuals fairly normal

    car::Anova(model_bc_final, type = 3)
    step(model_bc_final)
    #however findings are exactly the same - only 1 significant variable


    
### POWER ANALYSIS ###################################################################### 
    
    # Sample sizes for treatment and control groups
    n_0 <- 74
    n_1 <- 77
    
    #Expected proportions for both groups
    p_0 <- 0.25
    p_1 <- 0.50
    
    #We're examining effect size for proportion, so we need to calculate Cohen's h
    h <- ES.h(p_1, p_0)
    h #by rule of thumb, 0.52 is a moderate effect
    
    #The test in question is proportion test for 2 samples of different sizes,
    # so we use pwr.2p.test function from pwr package, two-sided test
    pwr <- pwr.2p2n.test(h = h, n1 = n_0, n2 = n_1, sig.level = 0.05, power = NULL, alternative = "two.sided")
    pwr$power
    plot(pwr)
    #we get a power of approx. 0.90
    
    #let's calculate now the power function which is a relationship between power 
    #and effect size
    
    # Define a range of expected responses to a treatment
    props <- seq(0.25, 1, by = 0.01)
    hs <- ES.h(props, 0.25)
    
    index_h = which(hs > h)[1] - 1
    
    # Calculate statistical power for each effect size
    power_values <- sapply(hs, function(h) {
      pwr.2p2n.test(h = h, n1 = n_0, n2 = n_1, sig.level = 0.05, power = NULL, alternative = "two.sided")$power
    })
    
    
    # Power plot for parameter
    plot(x = props, y = power_values,
         type = "l", col = "blue", lwd = 1,
         xlab = "Effect Size", ylab = "Statistical Power",
         main = "Power Function Plot", ylim = c(0,1))
    abline(h = power_values[index_h], col = "red", lty = 2)  # Add a line for 80% power
    abline(v = 0.5, col = "blue", lty = 2)
    
    # Power plot for effect size
    plot(x = hs, y = power_values,
         type = "l", col = "blue", lwd = 1,
         xlab = "Effect Size", ylab = "Statistical Power",
         main = "Power Function Plot", ylim = c(0,1))
    abline(h = power_values[index_h], col = "red", lty = 2)  # Add a line for 80% power
    abline(v = h, col = "blue", lty = 2)
    
    
    # SIMULATION ------------------------------------------------------------
    pvalues <- numeric(10000)
    for(i in 1:length(pvalues)){
      sample1 <- rbinom(n_1, 1, 0.5)
      sample0 <- rbinom(n_0, 1, 0.25)
      test <- prop.test(x = c(sum(sample1), sum(sample0)), n = c(n_1, n_0))
      pvalues[i] <- test$p.value
    }
    mean(pvalues < 0.05)
    #slightly lower than before, = 0.86
    #but the order of magnitude is correct, perhaps there are differences between
    #pwr.2p2n.test test and prop.test
    
    



# SVA_plot(typ = "conv",       
#          dane = df_model %>% mutate(E = 1),
#          zmienna = "FSS_0m",         
#          #predykcje = c("pred"),  
#          Y = 'FSS_diff_24m',                  
#          ETR = "E",                
#          n_szkod = "E",
#          exp_order = F,
#          paczkuj = T,
#          n_bin = 12,
#          top_n = 100,
#          inter = T,       #czy pokazywac przekroj w roznych orkesach czasu
#          zmienna_int = "treatment_group", log_transform = F)

