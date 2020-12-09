##################packages##################
library(readxl)
library(MASS)
library(ggplot2)
library(lme4)
library(readxl)
library(ggpubr)
##################load data##################
ne_data <- read_excel("Trachy Novel Environment.xlsx", sheet = "Calculations", col_types = c("text", "text", "text", "numeric", "numeric", "numeric", "numeric", "numeric"))
View(ne_data)

##################box and whisker plot##################
ggplot(ne_data, aes(x = species, y = tempo, color = colony)) +
  geom_boxplot(color="black", width = 0.3) +
  geom_jitter(position = position_jitter(0.1)) +
  scale_color_manual(values=c("#FFD800", "#A07C00", "#CB4949", "#EC0808", "#00480A", "#00FF67", "#447140", "#85EA2B")) +
  scale_x_discrete(breaks=c("BR", "LB"), labels=c("T. zeteki", "T. cf. zeteki")) +
  ggtitle("Tempo") +
  theme(plot.title = element_text(hjust = 0.5))

##################figure out which GLMM to use##################
#models
tempo.lmm <- lmer(tempo ~ species + (1|colony), data = ne_data, REML=T)
tempo.glmm <- glmer(tempo ~ species + (1|colony), data = ne_data, family = gaussian(link = log))
tempo.glmm2 <- glmer(tempo ~ species + (1|colony), data = ne_data, family = Gamma(link = log))
tempo.glmm3 <- glmer(tempo ~ species + (1|colony), data = ne_data, family = Gamma(link = identity))
tempo.glmm4 <- glmer(tempo ~ species + (1|colony), data = ne_data, family = Gamma(link = inverse))
#compare AIC values to determine which model to use
AIC(tempo.lmm, tempo.glmm, tempo.glmm2, tempo.glmm3, tempo.glmm4)
#best model
best_summary<-summary(tempo.glmm4)

##################check fit of GLMM with Gamma distribution and inverse link function##################
#test normality of residuals Shap-Wilk, qqplot, densplot
tempo.resid <- residuals(tempo.glmm4)
shapiro.test(tempo.resid)
ggqqplot(tempo.resid, main="QQ-plot of Best GLMM")
ggdensity(tempo.resid, main="Density Plot of Tempo", xlab="Tempo")

#test for homoscedasticity
plot(fitted(tempo.glmm4), residuals(tempo.glmm4), main="GLMM Residuals", xlab="Fitted", ylab="Residuals")
abline(h = 0, b = 2)
ggplot(fitted(tempo.glmm4), residuals(tempo.glmm4))
lines(smooth.spline(fitted(tempo.glmm4), residuals(tempo.glmm4)))


#check GLMM against a reduced model without random effect
tempo.glm <- glm(tempo ~ species, data = ne_data, family = Gamma(link = inverse))
tempo.colony.p <- pchisq(as.numeric(-2 * logLik(tempo.glm) + 2 * logLik(tempo.glmm4)), df=1, lower.tail=F)/2


################output tables########################
library(stargazer)
all_models<-stargazer(tempo.glm,tempo.glmm,tempo.glmm2,tempo.glmm3,tempo.glmm4,summary=TRUE,rownames=TRUE,align=TRUE,type="html",out="deg.html")

stargazer(glmm_final,rownames=TRUE,align=TRUE,type="html",out="glmm4.html")