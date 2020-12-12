####load data and packages####

#data_full <- readRDS("data/working_full.RDS")

data_full <- readRDS("data/working.RDS")

library(tidyr)
library(dplyr)
library(ggplot2)
library(GGally)
library(plm)
library(stargazer)
library(lmtest)
library(sandwich)
library(margins)
library(corrplot)
library(ggeffects)
library(splines)
library(ggpubr)

####visualize and explore####

attach(data_full)

ggplot(data = data_full, aes(x=log(poverty19+1),e_polity2))+
  geom_point()+
  geom_smooth()+
  theme_light()

hist(log(poverty19))
hist(poverty19)
hist(log(poverty_nat))
hist(log(GDP_pc_ppp))
hist(log(population_WB))
hist(log(gini_WB))
hist(v2x_regime)

gplots::plotmeans(e_polity2~year, data=data_full, barcol = "gray")
gplots::plotmeans(poverty19~year, data=data_full, ylim=c(-5, 50), barcol = "gray")

ggplot(data_full, aes(y=e_polity2, x = poverty19, color=c_name), show.legend = FALSE)+
  geom_smooth(method=lm, se=FALSE)+
  theme(legend.position = "none")+
  theme_light()


####creat pca regime index####

library("haven")
library("factoextra")
library("psych")

regimes <- data_full[!is.na(data_full$e_polity2) & !is.na(data_full$v2x_regime) & 
                       !is.na(data_full$e_fh_status),]

pca_res <- prcomp(regimes[c("e_polity2", "v2x_regime", "e_fh_status")], center = TRUE, scale = TRUE)
summary(pca_res)

eig.value <- get_eigenvalue(pca_res)
eig.value

pca_res$rotation[,1:2]

fviz_contrib(pca_res, choice = "var", axes = 1, top = 10)
fviz_contrib(pca_res, choice = "var", axes = 2, top = 10)

fviz_pca_var(pca_res, col.var = "contrib", repel = TRUE)          

eig.value
fviz_eig(pca_res, addlabels = TRUE)

regimes$regime <- 0.5814937*regimes$e_polity2 + 
  0.5783873*regimes$v2x_regime - 0.5721305*regimes$e_fh_status

plot(regimes$regime~regimes$e_polity2)

regimes <- regimes[c("c_name", "year", "c_code", "regime")]

data_full <- plyr::join_all(list(data_full, regimes), by = c("c_code", "year"), type = "full")

rm(regimes, eig.value, pca_res)

stargazer(data_full)


paneldata <- pdata.frame(data_full, index = c("c_name", "year"))

paneldata$polity2_lag <- lag(paneldata$e_polity2)
paneldata$regime_lag <- lag(paneldata$regime)

ggpairs(data_full[4:19], theme="light")

M <- cor(data_full[4:19],use="complete.obs")
corrplot(M, method = "circle")

####hypothesis 1, lagged dependent variable####

h1.n.p.l <- plm(polity2_lag~poverty_nat+log(GDP_pc_ppp)+log(population_WB)+
                  v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="pooling", data=paneldata)

summary(h1.n.p.l)


h1.19.p.l <- plm(polity2_lag~log(poverty19+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="pooling", data=paneldata)

summary(h1.19.p.l)

h1.32.p.l <- plm(polity2_lag~log(poverty32+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="pooling", data=paneldata) 

summary(h1.32.p.l)

h1.55.p.l <- plm(polity2_lag~log(poverty55+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="pooling", data=paneldata) 

summary(h1.55.p.l)

h1.19g.p.l <- plm(polity2_lag~log(poverty_gap_19+1)+log(GDP_pc_ppp)+
                    log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="pooling", data=paneldata)

summary(h1.19g.p.l)

h1.32g.p.l <- plm(polity2_lag~log(poverty_gap_32+1)+log(GDP_pc_ppp)+
                    log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="pooling", data=paneldata) 

summary(h1.32g.p.l)

h1.55g.p.l <- plm(polity2_lag~log(poverty_gap_55+1)+log(GDP_pc_ppp)+
                    log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="pooling", data=paneldata) 

summary(h1.55g.p.l)

se1pl <- c(list(sqrt(diag(vcovHC(h1.n.p.l, type = "HC1"))),
               sqrt(diag(vcovHC(h1.19.p.l, type = "HC1"))),
               sqrt(diag(vcovHC(h1.32.p.l, type = "HC1"))),
               sqrt(diag(vcovHC(h1.55.p.l, type = "HC1"))),
               sqrt(diag(vcovHC(h1.19g.p.l, type = "HC1"))),
               sqrt(diag(vcovHC(h1.32g.p.l, type = "HC1"))),
               sqrt(diag(vcovHC(h1.55g.p.l, type = "HC1")))))

stargazer(h1.n.p.l, h1.19.p.l, h1.32.p.l, h1.55.p.l, h1.19g.p.l, h1.32g.p.l,
          h1.55g.p.l, se = se1pl)


####one-way####

h1.n.o.l <- plm(polity2_lag~poverty_nat+log(GDP_pc_ppp)+log(population_WB)+
                  v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="within", data=paneldata)

summary(h1.n.o.l)
coeftest(h1.n.o.l, vcov = vcovHC, type = "HC3")

h1.19.o.l <- plm(polity2_lag~log(poverty19+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata)

summary(h1.19.o.l)
coeftest(h1.19.o.l, vcov = vcovHC, type = "HC3")

h1.32.o.l <- plm(polity2_lag~log(poverty32+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata) 

summary(h1.32.o.l)
coeftest(h1.32.o.l, vcov = vcovHC, type = "HC3")

h1.55.o.l <- plm(polity2_lag~log(poverty55+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata) 

summary(h1.55.o.l)
coeftest(h1.55.o.l, vcov = vcovHC, type = "HC3")

h1.19g.o.l <- plm(polity2_lag~log(poverty_gap_19+1)+log(GDP_pc_ppp)+log(population_WB)+
                    v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata)

summary(h1.19g.o.l)
coeftest(h1.19g.o.l, vcov = vcovHC, type = "HC3")

h1.32g.o.l <- plm(polity2_lag~log(poverty_gap_32+1)+log(GDP_pc_ppp)+log(population_WB)+
                    v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata) 

summary(h1.32g.o.l)
coeftest(h1.32g.o.l, vcov = vcovHC, type = "HC3")

h1.55g.o.l <- plm(polity2_lag~log(poverty_gap_55+1)+log(GDP_pc_ppp)+log(population_WB)+
                    v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata) 

summary(h1.55g.o.l)
coeftest(h1.19.o.l, vcov = vcovHC, type = "HC3")

pFtest(h1.n.o.l, h1.n.p.l)
pFtest(h1.19.o.l, h1.19.p.l)
pFtest(h1.32.o.l, h1.32.p.l)
pFtest(h1.55.o.l, h1.55.p.l)
pFtest(h1.19g.o.l, h1.19g.p.l)
pFtest(h1.32g.o.l, h1.32g.p.l)
pFtest(h1.55g.o.l, h1.55g.p.l)

se1ol <- c(list(sqrt(diag(vcovHC(h1.n.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.19.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.32.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.55.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.19g.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.32g.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.55g.o.l, type = "HC1")))))

stargazer(h1.n.o.l, h1.19.o.l, h1.32.o.l, h1.55.o.l, h1.19g.o.l, h1.32g.o.l,
          h1.55g.o.l, se = se1ol)

####random####

h1.n.r.l <- plm(polity2_lag~poverty_nat+log(GDP_pc_ppp)+log(population_WB)+
                  v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="random", data=paneldata)

summary(h1.n.r.l)
coeftest(h1.n.r.l, vcov = vcovHC, type = "HC3")

h1.19.r.l <- plm(polity2_lag~log(poverty19+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="random", data=paneldata)

summary(h1.19.r.l)
coeftest(h1.19.r.l, vcov = vcovHC, type = "HC3")

h1.32.r.l <- plm(polity2_lag~log(poverty32+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="random", data=paneldata) 

summary(h1.32.r.l)
coeftest(h1.32.r.l, vcov = vcovHC, type = "HC3")

h1.55.r.l <- plm(polity2_lag~log(poverty55+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="random", data=paneldata) 

summary(h1.55.r.l)
coeftest(h1.55.r.l, vcov = vcovHC, type = "HC3")

h1.19g.r.l <- plm(polity2_lag~log(poverty_gap_19+1)+log(GDP_pc_ppp)+log(population_WB)+
                    v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="random", data=paneldata)

summary(h1.19g.r.l)
coeftest(h1.19g.r.l, vcov = vcovHC, type = "HC3")

h1.32g.r.l <- plm(polity2_lag~log(poverty_gap_32+1)+log(GDP_pc_ppp)+log(population_WB)+
                    v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="random", data=paneldata) 

summary(h1.32g.r.l)
coeftest(h1.32g.r.l, vcov = vcovHC, type = "HC3")

h1.55g.r.l <- plm(polity2_lag~log(poverty_gap_55+1)+log(GDP_pc_ppp)+log(population_WB)+
                    v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="random", data=paneldata) 

summary(h1.55g.r.l)
coeftest(h1.19.r.l, vcov = vcovHC, type = "HC3")

phtest(h1.n.o.l, h1.n.r.l)
phtest(h1.19.o.l, h1.19.r.l)
phtest(h1.32.o.l, h1.32.r.l)
phtest(h1.55.o.l, h1.55.r.l)
phtest(h1.19g.o.l, h1.19g.r.l)
phtest(h1.32g.o.l, h1.32g.r.l)
phtest(h1.55g.o.l, h1.55g.r.l)

####two-ways####

h1.n.t.l <- plm(polity2_lag~poverty_nat+log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+
                  v2xeg_eqdr+log(gini_WB), 
                model="within", effect="twoways", data=paneldata)

summary(h1.n.t.l)
coeftest(h1.n.t.l, vcov = vcovHC, type = "HC3")

h1.19.t.l <- plm(polity2_lag~log(poverty19+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata)

summary(h1.19.t.l)
coeftest(h1.19.t.l, vcov = vcovHC, type = "HC3")

h1.32.t.l <- plm(polity2_lag~log(poverty32+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata) 

summary(h1.32.t.l)
coeftest(h1.32.t.l, vcov = vcovHC, type = "HC3")

h1.55.t.l <- plm(polity2_lag~log(poverty55+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata) 

summary(h1.55.t.l)
coeftest(h1.55.t.l, vcov = vcovHC, type = "HC3")

h1.19g.t.l <- plm(polity2_lag~log(poverty_gap_19+1)+log(GDP_pc_ppp)+log(population_WB)+
                    v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata)

summary(h1.19g.t.l)
coeftest(h1.19g.t.l, vcov = vcovHC, type = "HC3")

h1.32g.t.l <- plm(polity2_lag~log(poverty_gap_32+1)+log(GDP_pc_ppp)+log(population_WB)+
                    v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata) 

summary(h1.32g.t.l)
coeftest(h1.32g.t.l, vcov = vcovHC, type = "HC3")

h1.55g.t.l <- plm(polity2_lag~log(poverty_gap_55+1)+log(GDP_pc_ppp)+log(population_WB)+
                    v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata) 

summary(h1.55g.t.l)
coeftest(h1.55g.t.l, vcov = vcovHC, type = "HC3")

pFtest(h1.n.t.l, h1.n.o.l)
pFtest(h1.19.t.l, h1.19.o.l)
pFtest(h1.32.t.l, h1.32.o.l)
pFtest(h1.55.t.l, h1.55.o.l)
pFtest(h1.19g.t.l, h1.19g.o.l)
pFtest(h1.32g.t.l, h1.32g.o.l)
pFtest(h1.55g.t.l, h1.55g.o.l)

se1tl <- c(list(sqrt(diag(vcovHC(h1.n.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.19.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.32.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.55.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.19g.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.32g.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h1.55g.t.l, type = "HC1")))))

stargazer(h1.n.t.l, h1.19.t.l, h1.32.t.l, h1.55.t.l, h1.19g.t.l, h1.32g.t.l,
          h1.55g.t.l, se = se1tl)

##truncated sample####

dffp <- data.frame(paneldata)
dffp <- dffp[c("c_name", "c_code", "year", "poverty19", "poverty_gap_19", "polity2_lag",
               "GDP_pc_ppp", "population_WB", "v2regsupgroupssize", "v2xeg_eqdr",
               "gini_WB")]
dffp <- na.omit(dffp)

h2t1 <- lm(polity2_lag~log(poverty19+1)+log(GDP_pc_ppp)+log(population_WB)+
             v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+factor(c_name)+factor(year),
           data = dffp)
h2t2 <- lm(polity2_lag~log(poverty_gap_19+1)+log(GDP_pc_ppp)+log(population_WB)+
             v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+factor(c_name)+factor(year),
           data = dffp)
summary(h2t1)

dffp$c_dummy <- as.numeric(dffp$c_name)
#dffp$c_dummy <- fastDummies::dummy_cols(dffp$c_name)

y_pred1 <- h2t1$fitted.values
y_pred2 <- h2t2$fitted.values
panel1 <- data.frame(dffp, y_pred1)
merged <- panel1 %>%
  group_by(c_dummy)%>%
  summarize(., cor(polity2_lag, y_pred1))%>%
  merge(panel1, ., by="c_dummy")

merged$new <- ifelse(abs(merged$`cor(polity2_lag, y_pred1)`)<0.3,1,0)
fe_twoways_1 <- plm(polity2_lag~log(poverty19+1)+log(GDP_pc_ppp)+log(population_WB)+
                      v2regsupgroupssize+v2xeg_eqdr+log(gini_WB) , merged[merged$new == 0,],
                    index=c("c_name", "year"), effect = "twoways")
coeftest(fe_twoways_1, vcov = vcovHC, type = "HC3")

panel2 <- data.frame(dffp, y_pred2)
merged <- panel2 %>%
  group_by(c_dummy)%>%
  summarize(., cor(polity2_lag, y_pred2))%>%
  merge(panel1, ., by="c_dummy")

merged$new <- ifelse(abs(merged$`cor(polity2_lag, y_pred2)`)<0.3,1,0)
fe_twoways_2 <- plm(polity2_lag~log(poverty_gap_19+1)+log(GDP_pc_ppp)+log(population_WB)+
                      v2regsupgroupssize+v2xeg_eqdr+log(gini_WB) , merged[merged$new == 0,],
                    index=c("c_name", "year"), effect = "twoways")
coeftest(fe_twoways_2, vcov = vcovHC, type = "HC3")

se1trl <- c(list(sqrt(diag(vcovHC(fe_twoways_1, type = "HC1"))),
                sqrt(diag(vcovHC(fe_twoways_2, type = "HC1")))))

stargazer(fe_twoways_1, fe_twoways_2, se = se1trl)


####hypothesis 2: quadratic term, lagged dependent variable####

h2.n.p.l <- plm(polity2_lag~poverty_nat+I(poverty_nat^2)+log(GDP_pc_ppp)+
                  log(population_WB)+log(gini_WB), 
                model="pooling", data=paneldata)

summary(h2.n.p.l)


h2.19.p.l <- plm(polity2_lag~log(poverty19+1)+I(log(poverty19+1)^2)+log(GDP_pc_ppp)+
                   log(population_WB)+log(gini_WB), 
                 model="pooling", data=paneldata)

summary(h2.19.p.l)

h2.32.p.l <- plm(polity2_lag~log(poverty32+1)+I(log(poverty32+1))+log(GDP_pc_ppp)+
                   log(population_WB)+log(gini_WB), 
                 model="pooling", data=paneldata) 

summary(h2.32.p.l)

h2.55.p.l <- plm(polity2_lag~log(poverty55+1)+I(log(poverty55+1)^2)+log(GDP_pc_ppp)+
                   log(population_WB)+log(gini_WB), 
                 model="pooling", data=paneldata) 

summary(h2.55.p.l)

h2.19g.p.l <- plm(polity2_lag~log(poverty_gap_19+1)+I(log(poverty_gap_19+1)^2)+log(GDP_pc_ppp)
                  +log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="pooling", data=paneldata)

summary(h2.19g.p.l)

h2.32g.p.l <- plm(polity2_lag~log(poverty_gap_32+1)+I(log(poverty_gap_32+1)^2)+log(GDP_pc_ppp)+
                    log(population_WB)+log(gini_WB), 
                  model="pooling", data=paneldata) 

summary(h2.32g.p.l)

h2.55g.p.l <- plm(polity2_lag~log(poverty_gap_55+1)+I(log(poverty_gap_55+1)^2)+log(GDP_pc_ppp)+
                    log(population_WB)+log(gini_WB), 
                  model="pooling", data=paneldata) 

summary(h2.55g.p.l)

####one-way####

h2.n.o.l <- plm(polity2_lag~poverty_nat+I(poverty_nat^2)+
                  log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="within", data=paneldata)

summary(h2.n.o.l)
coeftest(h2.n.o.l, vcov = vcovHC, type = "HC3")


data_f <- data.frame(paneldata)

data_f$year_f <- factor(data_f$year)
data_f$c_name_f <- factor(data_f$c_name)

h2.n.o.l.slm <- lm(polity2_lag~poverty_nat+I(poverty_nat^2)+log(GDP_pc_ppp)+
                     log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                   data=data_f)


ggeffect(h2.n.o.l.slm, terms = "poverty_nat")
mydf <- ggeffect(h2.n.o.l.slm, terms = "poverty_nat")
#mydf1 <- ggpredict(h2.n.o.l.slm, terms = "poverty19")
h2n<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_nat")+
  ggtitle("Marginal effect of poverty_nat")+
  theme_light()


h2.19.o.l <- plm(polity2_lag~log(poverty19+1)+I((poverty19+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata)

summary(h2.19.o.l)
coeftest(h2.19.o.l, vcov = vcovHC, type = "HC3")

h2.19.o.l.slm <- lm(polity2_lag~log(poverty19+1)+I(log(poverty19+1)^2)+log(GDP_pc_ppp)+
                     log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                   data=data_f)


ggeffect(h2.19.o.l.slm, terms = "poverty19")
mydf <- ggeffect(h2.19.o.l.slm, terms = "poverty19")
#mydf1 <- ggpredict(h2.n.o.l.slm, terms = "poverty19")
h219<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty19")+
  ggtitle("Marginal effect of poverty19")+
  theme_light()

h2.32.o.l <- plm(polity2_lag~log(poverty32+1)+I(log(poverty32+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata) 

summary(h2.32.o.l)
coeftest(h2.32.o.l, vcov = vcovHC, type = "HC3")

h2.32.o.l.slm <- lm(polity2_lag~log(poverty32+1)+I(log(poverty32+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(h2.32.o.l.slm, terms = "poverty32")
mydf <- ggeffect(h2.32.o.l.slm, terms = "poverty32")
#mydf1 <- ggpredict(h2.n.o.l.slm, terms = "poverty19")
h232<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty32")+
  ggtitle("Marginal effect of poverty32")+
  theme_light()

h2.55.o.l <- plm(polity2_lag~log(poverty55+1)+I(log(poverty55+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata) 

summary(h2.55.o.l)
coeftest(h2.55.o.l, vcov = vcovHC, type = "HC3")

h2.55.o.l.slm <- lm(polity2_lag~log(poverty55+1)+I(log(poverty55+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(h2.55.o.l.slm, terms = "poverty55")
mydf <- ggeffect(h2.55.o.l.slm, terms = "poverty55")
#mydf1 <- ggpredict(h2.n.o.l.slm, terms = "poverty19")
h255<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty55")+
  ggtitle("Marginal effect of poverty55")+
  theme_light()

h2.19g.o.l <- plm(polity2_lag~log(poverty_gap_19+1)+I(log(poverty_gap_19+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata)

summary(h2.19g.o.l)
coeftest(h2.19g.o.l, vcov = vcovHC, type = "HC3")

h2.19g.o.l.slm <- lm(polity2_lag~log(poverty_gap_19+1)+I(log(poverty_gap_19+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(h2.19g.o.l.slm, terms = "poverty_gap_19")
mydf <- ggeffect(h2.19g.o.l.slm, terms = "poverty_gap_19")
#mydf1 <- ggpredict(h2.n.o.l.slm, terms = "poverty19")
h219g<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_19")+
  ggtitle("Marginal effect of poverty_gap_19")+
  theme_light()

h2.32g.o.l <- plm(polity2_lag~log(poverty_gap_32+1)+I(log(poverty_gap_32+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata) 

summary(h2.32g.o.l)
coeftest(h2.32g.o.l, vcov = vcovHC, type = "HC3")

h2.32g.o.l.slm <- lm(polity2_lag~log(poverty_gap_32+1)+I(log(poverty_gap_32+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                     data=data_f)


ggeffect(h2.32g.o.l.slm, terms = "poverty_gap_32")
mydf <- ggeffect(h2.32g.o.l.slm, terms = "poverty_gap_32")
#mydf1 <- ggpredict(h2.n.o.l.slm, terms = "poverty32")
h232g<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_32")+
  ggtitle("Marginal effect of poverty_gap_32")+
  theme_light()

h2.55g.o.l <- plm(polity2_lag~log(poverty_gap_55+1)+I(log(poverty_gap_55+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata) 

summary(h2.55g.o.l)
coeftest(h2.19.o.l, vcov = vcovHC, type = "HC3")

h2.55g.o.l.slm <- lm(polity2_lag~log(poverty_gap_55+1)+I(log(poverty_gap_55+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                     data=data_f)


ggeffect(h2.55g.o.l.slm, terms = "poverty_gap_55")
mydf <- ggeffect(h2.55g.o.l.slm, terms = "poverty_gap_55")
#mydf1 <- ggpredict(h2.n.o.l.slm, terms = "poverty55")
h255g<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_55")+
  ggtitle("Marginal effect of poverty_gap_55")+
  theme_light()

figure2 <- ggarrange(h2n, h219, h232, h255, h219g, h232g, h255g,
                    ncol = 2, nrow = 4)
figure2


pFtest(h2.n.o.l, h2.n.p.l)
pFtest(h2.19.o.l, h2.19.p.l)
pFtest(h2.32.o.l, h2.32.p.l)
pFtest(h2.55.o.l, h2.55.p.l)
pFtest(h2.19g.o.l, h2.19g.p.l)
pFtest(h2.32g.o.l, h2.32g.p.l)
pFtest(h2.55g.o.l, h2.55g.p.l)

se2ol <- c(list(sqrt(diag(vcovHC(h2.n.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.19.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.32.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.55.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.19g.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.32g.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.55g.o.l, type = "HC1")))))

stargazer(h2.n.o.l, h2.19.o.l, h2.32.o.l, h2.55.o.l, h2.19g.o.l, h2.32g.o.l,
          h2.55g.o.l, se = se2ol)





####two-ways####

h2.n.t.l <- plm(polity2_lag~poverty_nat+I(poverty_nat^2)+
                  log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="within", effect="twoways", data=paneldata)

summary(h2.n.t.l)
coeftest(h2.n.t.l, vcov = vcovHC, type = "HC3")

h2.19.t.l <- plm(polity2_lag~log(poverty19+1)+I(log(poverty19+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata)

summary(h2.19.t.l)
coeftest(h2.19.t.l, vcov = vcovHC, type = "HC3")

h2.32.t.l <- plm(polity2_lag~log(poverty32+1)+I(log(poverty32+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata) 

summary(h2.32.t.l)
coeftest(h2.32.t.l, vcov = vcovHC, type = "HC3")

h2.55.t.l <- plm(polity2_lag~log(poverty55+1)+I(log(poverty55+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata) 

summary(h2.55.t.l)
coeftest(h2.55.t.l, vcov = vcovHC, type = "HC3")

h2.19g.t.l <- plm(polity2_lag~log(poverty_gap_19+1)+I(log(poverty_gap_19+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata)

summary(h2.19g.t.l)
coeftest(h2.19g.t.l, vcov = vcovHC, type = "HC3")

h2.32g.t.l <- plm(polity2_lag~log(poverty_gap_32+1)+I(log(poverty_gap_32+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata) 

summary(h2.32g.t.l)
coeftest(h2.32g.t.l, vcov = vcovHC, type = "HC3")

h2.55g.t.l <- plm(polity2_lag~log(poverty_gap_55+1)+I(log(poverty_gap_55+1)^2)
                  +log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata) 

summary(h2.55g.t.l)
coeftest(h2.55g.t.l, vcov = vcovHC, type = "HC3")


se2tl <- c(list(sqrt(diag(vcovHC(h2.n.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.19.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.32.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.55.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.19g.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.32g.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h2.55g.t.l, type = "HC1")))))

stargazer(h2.n.t.l, h2.19.t.l, h2.32.t.l, h2.55.t.l, h2.19g.t.l, h2.32g.t.l,
          h2.55g.t.l, se = se2tl)

h2.n.t.l.slm <- lm(polity2_lag~log(poverty_nat+1)+I(log(poverty_nat+1)^2)+log(GDP_pc_ppp)+
                     log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                   data=data_f)


ggeffect(h2.n.t.l.slm, terms = "poverty_nat")
mydf <- ggeffect(h2.n.t.l.slm, terms = "poverty_nat")
#mydf1 <- ggpredict(h2.n.t.l.slm, terms = "poverty55")
h2nt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_nat")+
  ggtitle("Marginal effect of poverty_nat")+
  theme_light()

h2.19.t.l.slm <- lm(polity2_lag~log(poverty19+1)+I(log(poverty19+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(h2.19.t.l.slm, terms = "poverty19")
mydf <- ggeffect(h2.19.t.l.slm, terms = "poverty19")
#mydf1 <- ggpredict(h2.n.t.l.slm, terms = "poverty55")
h219t <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty19")+
  ggtitle("Marginal effect of poverty19")+
  theme_light()

h2.32.t.l.slm <- lm(polity2_lag~log(poverty32+1)+I(log(poverty32+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(h2.32.t.l.slm, terms = "poverty32")
mydf <- ggeffect(h2.32.t.l.slm, terms = "poverty32")
#mydf1 <- ggpredict(h2.n.t.l.slm, terms = "poverty55")
h232t <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty32")+
  ggtitle("Marginal effect of poverty32")+
  theme_light()

h2.55.t.l.slm <- lm(polity2_lag~log(poverty55+1)+I(log(poverty55+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(h2.55.t.l.slm, terms = "poverty55")
mydf <- ggeffect(h2.55.t.l.slm, terms = "poverty55")
#mydf1 <- ggpredict(h2.n.t.l.slm, terms = "poverty55")
h255t <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty55")+
  ggtitle("Marginal effect of poverty55")+
  theme_light()

h2.19g.t.l.slm <- lm(polity2_lag~log(poverty_gap_19+1)+I(log(poverty_gap_19+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)
                     +c_name_f+year_f,
                     data=data_f)


ggeffect(h2.19g.t.l.slm, terms = "poverty_gap_19")
mydf <- ggeffect(h2.19g.t.l.slm, terms = "poverty_gap_19")
#mydf1 <- ggpredict(h2.n.t.l.slm, terms = "poverty19g")
h219gt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_19")+
  ggtitle("Marginal effect of poverty_gap_19")+
  theme_light()

h2.32g.t.l.slm <- lm(polity2_lag~log(poverty_gap_32+1)+I(log(poverty_gap_32+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+
                       c_name_f+year_f,
                     data=data_f)


ggeffect(h2.32g.t.l.slm, terms = "poverty_gap_32")
mydf <- ggeffect(h2.32g.t.l.slm, terms = "poverty_gap_32")
#mydf1 <- ggpredict(h2.n.t.l.slm, terms = "poverty32g")
h232gt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_32")+
  ggtitle("Marginal effect of poverty_gap_32")+
  theme_light()

h2.55g.t.l.slm <- lm(polity2_lag~log(poverty_gap_55+1)+I(log(poverty_gap_55+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)
                     +c_name_f+year_f,
                     data=data_f)


ggeffect(h2.55g.t.l.slm, terms = "poverty_gap_55")
mydf <- ggeffect(h2.55g.t.l.slm, terms = "poverty_gap_55")
#mydf1 <- ggpredict(h2.n.t.l.slm, terms = "poverty55g")
h255gt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_55")+
  ggtitle("Marginal effect of poverty_gap_55")+
  theme_light()

figure2t <- ggarrange(h2nt, h219t, h232t, h255t, h219gt, h232gt, h255gt,
                     ncol = 2, nrow = 4)
figure2t



#####hypothesis 3: interaction effect, lagged dependent variable####

h3.n.p.l <- plm(polity2_lag~poverty_nat+poverty_nat:log(gini_WB)+log(GDP_pc_ppp)+
                  log(population_WB)+log(gini_WB), 
                model="pooling", data=paneldata)

summary(h3.n.p.l)


h3.19.p.l <- plm(polity2_lag~log(poverty19+1)+log(poverty19+1):log(gini_WB)+log(GDP_pc_ppp)+
                   log(population_WB)+log(gini_WB), 
                 model="pooling", data=paneldata)

summary(h3.19.p.l)

h3.32.p.l <- plm(polity2_lag~log(poverty32+1)+I(log(poverty32+1))+log(GDP_pc_ppp)+
                   log(population_WB)+log(gini_WB), 
                 model="pooling", data=paneldata) 

summary(h3.32.p.l)

h3.55.p.l <- plm(polity2_lag~log(poverty55+1)+log(poverty55+1):log(gini_WB)+log(GDP_pc_ppp)+
                   log(population_WB)+log(gini_WB), 
                 model="pooling", data=paneldata) 

summary(h3.55.p.l)

h3.19g.p.l <- plm(polity2_lag~log(poverty_gap_19+1)+log(poverty_gap_19+1):log(gini_WB)
                  +log(GDP_pc_ppp)
                  +log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="pooling", data=paneldata)

summary(h3.19g.p.l)

h3.32g.p.l <- plm(polity2_lag~log(poverty_gap_32+1)+log(poverty_gap_32+1):log(gini_WB)
                  +log(GDP_pc_ppp)+
                    log(population_WB)+log(gini_WB), 
                  model="pooling", data=paneldata) 

summary(h3.32g.p.l)

h3.55g.p.l <- plm(polity2_lag~log(poverty_gap_55+1)+log(poverty_gap_55+1):log(gini_WB)
                  +log(GDP_pc_ppp)+
                    log(population_WB)+log(gini_WB), 
                  model="pooling", data=paneldata) 

summary(h3.55g.p.l)

####one-way####

h3.n.o.l <- plm(polity2_lag~poverty_nat+poverty_nat:log(gini_WB)+
                  log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="within", data=paneldata)

summary(h3.n.o.l)
coeftest(h3.n.o.l, vcov = vcovHC, type = "HC3")

h3.19.o.l <- plm(polity2_lag~log(poverty19+1)+log(poverty19+1):log(gini_WB)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata)

summary(h3.19.o.l)
coeftest(h3.19.o.l, vcov = vcovHC, type = "HC3")

h3.32.o.l <- plm(polity2_lag~log(poverty32+1)+log(poverty32+1):log(gini_WB)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata) 

summary(h3.32.o.l)
coeftest(h3.32.o.l, vcov = vcovHC, type = "HC3")

h3.55.o.l <- plm(polity2_lag~log(poverty55+1)+log(poverty55+1):log(gini_WB)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata) 

summary(h3.55.o.l)
coeftest(h3.55.o.l, vcov = vcovHC, type = "HC3")

h3.19g.o.l <- plm(polity2_lag~log(poverty_gap_19+1)+log(poverty_gap_19+1):log(gini_WB)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata)

summary(h3.19g.o.l)
coeftest(h3.19g.o.l, vcov = vcovHC, type = "HC3")

h3.32g.o.l <- plm(polity2_lag~log(poverty_gap_32+1)+log(poverty_gap_32+1):log(gini_WB)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata) 

summary(h3.32g.o.l)
coeftest(h3.32g.o.l, vcov = vcovHC, type = "HC3")

h3.55g.o.l <- plm(polity2_lag~log(poverty_gap_55+1)+log(poverty_gap_55+1):log(gini_WB)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata) 

summary(h3.55g.o.l)
coeftest(h3.19.o.l, vcov = vcovHC, type = "HC3")

se3ol <- c(list(sqrt(diag(vcovHC(h3.n.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.19.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.32.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.55.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.19g.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.32g.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.55g.o.l, type = "HC1")))))

stargazer(h3.n.o.l, h3.19.o.l, h3.32.o.l, h3.55.o.l, h3.19g.o.l, h3.32g.o.l,
          h3.55g.o.l, se = se3ol)

h3.19.o.l.slm <- lm(polity2_lag~poverty_nat+poverty_nat*log(gini_WB)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(h3.n.o.l.slm, terms = "poverty_nat")
mydf <- ggeffect(h3.n.o.l.slm, terms = "poverty_nat")
#mydf1 <- ggpredict(h3.n.o.l.slm, terms = "poverty55")
h3no <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_nat")+
  ggtitle("Marginal effect of poverty_nat")+
  theme_light()

h3.19.o.l.slm <- lm(polity2_lag~log(poverty19+1)+log(poverty19+1)*log(gini_WB)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(h3.19.o.l.slm, terms = "poverty19")
mydf <- ggeffect(h3.19.o.l.slm, terms = "poverty19")
#mydf1 <- ggpredict(h3.n.o.l.slm, terms = "poverty55")
h319o <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty19")+
  ggtitle("Marginal effect of poverty19")+
  theme_light()

h3.32.o.l.slm <- lm(polity2_lag~log(poverty32+1)+log(poverty32+1)*log(gini_WB)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(h3.32.o.l.slm, terms = "poverty32")
mydf <- ggeffect(h3.32.o.l.slm, terms = "poverty32")
#mydf1 <- ggpredict(h3.n.o.l.slm, terms = "poverty55")
h332o <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty32")+
  ggtitle("Marginal effect of poverty32")+
  theme_light()

h3.55.o.l.slm <- lm(polity2_lag~log(poverty55+1)+log(poverty55+1)*log(gini_WB)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(h3.55.o.l.slm, terms = "poverty55")
mydf <- ggeffect(h3.55.o.l.slm, terms = "poverty55")
#mydf1 <- ggpredict(h3.n.o.l.slm, terms = "poverty55")
h355o <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty55")+
  ggtitle("Marginal effect of poverty55")+
  theme_light()

h3.19g.o.l.slm <- lm(polity2_lag~log(poverty_gap_19+1)+log(poverty_gap_19+1)*log(gini_WB)
                     +log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                     data=data_f)


ggeffect(h3.19g.o.l.slm, terms = "poverty_gap_19")
mydf <- ggeffect(h3.19g.o.l.slm, terms = "poverty_gap_19")
#mydf1 <- ggpredict(h3.n.o.l.slm, terms = "poverty19g")
h319go <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_19")+
  ggtitle("Marginal effect of poverty_gap_19")+
  theme_light()

h3.32g.o.l.slm <- lm(polity2_lag~log(poverty_gap_32+1)+log(poverty_gap_32+1)*log(gini_WB)
                     +log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                     data=data_f)


ggeffect(h3.32g.o.l.slm, terms = "poverty_gap_32")
mydf <- ggeffect(h3.32g.o.l.slm, terms = "poverty_gap_32")
#mydf1 <- ggpredict(h3.n.o.l.slm, terms = "poverty32g")
h332go <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_32")+
  ggtitle("Marginal effect of poverty_gap_32")+
  theme_light()

h3.55g.o.l.slm <- lm(polity2_lag~log(poverty_gap_55+1)+log(poverty_gap_55+1)*log(gini_WB)
                     +log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                     data=data_f)


ggeffect(h3.55g.o.l.slm, terms = "poverty_gap_55")
mydf <- ggeffect(h3.55g.o.l.slm, terms = "poverty_gap_55")
#mydf1 <- ggpredict(h3.n.o.l.slm, terms = "poverty55g")
h355go <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_55")+
  ggtitle("Marginal effect of poverty_gap_55")+
  theme_light()

figure3o <- ggarrange(h3no, h319o, h332o, h355o, h319go, h332go, h355go,
                      ncol = 2, nrow = 4)
figure3o





####two-ways####

h3.n.t.l <- plm(polity2_lag~poverty_nat+poverty_nat:log(gini_WB)+
                  log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="within", effect="twoways", data=paneldata)

summary(h3.n.t.l)
coeftest(h3.n.t.l, vcov = vcovHC, type = "HC3")
summary(margins(h3.n.t.l))

h3.19.t.l <- plm(polity2_lag~log(poverty19+1)+log(poverty19+1):log(gini_WB)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata)

summary(h3.19.t.l)
coeftest(h3.19.t.l, vcov = vcovHC, type = "HC3")

h3.32.t.l <- plm(polity2_lag~log(poverty32+1)+log(poverty32+1):log(gini_WB)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata) 

summary(h3.32.t.l)
coeftest(h3.32.t.l, vcov = vcovHC, type = "HC3")

h3.55.t.l <- plm(polity2_lag~log(poverty55+1)+log(poverty55+1):log(gini_WB)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata) 

summary(h3.55.t.l)
coeftest(h3.55.t.l, vcov = vcovHC, type = "HC3")

h3.19g.t.l <- plm(polity2_lag~log(poverty_gap_19+1)+log(poverty_gap_19+1):log(gini_WB)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata)

summary(h3.19g.t.l)
coeftest(h3.19g.t.l, vcov = vcovHC, type = "HC3")

h3.32g.t.l <- plm(polity2_lag~log(poverty_gap_32+1)+log(poverty_gap_32+1):log(gini_WB)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata) 

summary(h3.32g.t.l)
coeftest(h3.32g.t.l, vcov = vcovHC, type = "HC3")

h3.55g.t.l <- plm(polity2_lag~log(poverty_gap_55+1)+log(poverty_gap_55+1):log(gini_WB)
                  +log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata) 

summary(h3.55g.t.l)
coeftest(h3.55g.t.l, vcov = vcovHC, type = "HC3")


se3tl <- c(list(sqrt(diag(vcovHC(h3.n.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.19.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.32.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.55.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.19g.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.32g.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(h3.55g.t.l, type = "HC1")))))

stargazer(h3.n.t.l, h3.19.t.l, h3.32.t.l, h3.55.t.l, h3.19g.t.l, h3.32g.t.l,
          h3.55g.t.l, se = se3tl)



h3.n.t.l.slm <- lm(polity2_lag~poverty_nat+poverty_nat*log(gini_WB)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(h3.n.t.l.slm, terms = "poverty_nat")
mydf <- ggeffect(h3.n.t.l.slm, terms = "poverty_nat")
#mydf1 <- ggpredict(h3.n.t.l.slm, terms = "poverty55")
h3nt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_nat")+
  ggtitle("Marginal effect of poverty_nat")+
  theme_light()

h3.19.t.l.slm <- lm(polity2_lag~log(poverty19+1)+log(poverty19+1)*log(gini_WB)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(h3.19.t.l.slm, terms = "poverty19")
mydf <- ggeffect(h3.19.t.l.slm, terms = "poverty19")
#mydf1 <- ggpredict(h3.n.t.l.slm, terms = "poverty55")
h319t <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty19")+
  ggtitle("Marginal effect of poverty19")+
  theme_light()

h3.32.t.l.slm <- lm(polity2_lag~log(poverty32+1)+log(poverty32+1)*log(gini_WB)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(h3.32.t.l.slm, terms = "poverty32")
mydf <- ggeffect(h3.32.t.l.slm, terms = "poverty32")
#mydf1 <- ggpredict(h3.n.t.l.slm, terms = "poverty55")
h332t <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty32")+
  ggtitle("Marginal effect of poverty32")+
  theme_light()

h3.55.t.l.slm <- lm(polity2_lag~log(poverty55+1)+log(poverty55+1)*log(gini_WB)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(h3.55.t.l.slm, terms = "poverty55")
mydf <- ggeffect(h3.55.t.l.slm, terms = "poverty55")
#mydf1 <- ggpredict(h3.n.t.l.slm, terms = "poverty55")
h355t <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty55")+
  ggtitle("Marginal effect of poverty55")+
  theme_light()

h3.19g.t.l.slm <- lm(polity2_lag~log(poverty_gap_19+1)+log(poverty_gap_19+1)*log(gini_WB)
                     +log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)
                     +c_name_f+year_f,
                     data=data_f)


ggeffect(h3.19g.t.l.slm, terms = "poverty_gap_19")
mydf <- ggeffect(h3.19g.t.l.slm, terms = "poverty_gap_19")
#mydf1 <- ggpredict(h3.n.t.l.slm, terms = "poverty19g")
h319gt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_19")+
  ggtitle("Marginal effect of poverty_gap_19")+
  theme_light()

h3.32g.t.l.slm <- lm(polity2_lag~log(poverty_gap_32+1)+log(poverty_gap_32+1)*log(gini_WB)
                     +log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)
                     +c_name_f+year_f,
                     data=data_f)


ggeffect(h3.32g.t.l.slm, terms = "poverty_gap_32")
mydf <- ggeffect(h3.32g.t.l.slm, terms = "poverty_gap_32")
#mydf1 <- ggpredict(h3.n.t.l.slm, terms = "poverty32g")
h332gt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_32")+
  ggtitle("Marginal effect of poverty_gap_32")+
  theme_light()

h3.55g.t.l.slm <- lm(polity2_lag~log(poverty_gap_55+1)+log(poverty_gap_55+1)*log(gini_WB)
                     +log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)
                     +c_name_f+year_f,
                     data=data_f)


ggeffect(h3.55g.t.l.slm, terms = "poverty_gap_55")
mydf <- ggeffect(h3.55g.t.l.slm, terms = "poverty_gap_55")
#mydf1 <- ggpredict(h3.n.t.l.slm, terms = "poverty55g")
h355gt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_55")+
  ggtitle("Marginal effect of poverty_gap_55")+
  theme_light()

figure3t <- ggarrange(h3nt, h319t, h332t, h355t, h319gt, h332gt, h355gt,
                      ncol = 2, nrow = 4)
figure3t



####robustness####

r.19.t.l <- plm(regime_lag~log(poverty19+1)+log(GDP_pc_ppp)+log(population_WB)+
                  v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="within", effect="twoways", data=paneldata)

summary(r.19.t.l)
coeftest(r.19.t.l, vcov = vcovHC, type = "HC3")

r.19g.t.l <- plm(regime_lag~log(poverty_gap_19+1)+log(GDP_pc_ppp)+log(population_WB)+
                   v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata)

summary(r.19g.t.l)
coeftest(r.19g.t.l, vcov = vcovHC, type = "HC3")

dffp <- data.frame(paneldata)
dffp <- dffp[c("c_name", "c_code", "year", "poverty19", "poverty_gap_19", "regime_lag",
               "GDP_pc_ppp", "population_WB", "v2regsupgroupssize", "v2xeg_eqdr",
               "gini_WB")]
dffp <- na.omit(dffp)

rt1 <- lm(regime_lag~log(poverty19+1)+log(GDP_pc_ppp)+log(population_WB)+
            v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+factor(c_name)+factor(year),
          data = dffp)
rt2 <- lm(regime_lag~log(poverty_gap_19+1)+log(GDP_pc_ppp)+log(population_WB)+
            v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+factor(c_name)+factor(year),
          data = dffp)
summary(rt1)

dffp$c_dummy <- as.numeric(dffp$c_name)
#dffp$c_dummy <- fastDummies::dummy_cols(dffp$c_name)

y_pred1 <- rt1$fitted.values
y_pred2 <- rt2$fitted.values
panel1 <- data.frame(dffp, y_pred1)
merged <- panel1 %>%
  group_by(c_dummy)%>%
  summarize(., cor(regime_lag, y_pred1))%>%
  merge(panel1, ., by="c_dummy")

merged$new <- ifelse(abs(merged$`cor(regime_lag, y_pred1)`)<0.3,1,0)
fe_twoways_1 <- plm(regime_lag~log(poverty19+1)+log(GDP_pc_ppp)+log(population_WB)+
                      v2regsupgroupssize+v2xeg_eqdr+log(gini_WB) , merged[merged$new == 0,],
                    index=c("c_name", "year"), effect = "twoways")
coeftest(fe_twoways_1, vcov = vcovHC, type = "HC3")

panel2 <- data.frame(dffp, y_pred2)
merged <- panel2 %>%
  group_by(c_dummy)%>%
  summarize(., cor(regime_lag, y_pred2))%>%
  merge(panel1, ., by="c_dummy")

merged$new <- ifelse(abs(merged$`cor(regime_lag, y_pred2)`)<0.3,1,0)
fe_twoways_2 <- plm(regime_lag~log(poverty_gap_19+1)+log(GDP_pc_ppp)+log(population_WB)+
                      v2regsupgroupssize+v2xeg_eqdr+log(gini_WB) , merged[merged$new == 0,],
                    index=c("c_name", "year"), effect = "twoways")
coeftest(fe_twoways_2, vcov = vcovHC, type = "HC3")

se1trl <- c(list(sqrt(diag(vcovHC(fe_twoways_1, type = "HC1"))),
                 sqrt(diag(vcovHC(fe_twoways_2, type = "HC1")))))

stargazer(fe_twoways_1, fe_twoways_2, se = se1trl)

r2.n.o.l <- plm(regime_lag~poverty_nat+I(poverty_nat^2)+
                  log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="within", data=paneldata)

summary(r2.n.o.l)
coeftest(r2.n.o.l, vcov = vcovHC, type = "HC3")


data_f <- data.frame(paneldata)

data_f$year_f <- factor(data_f$year)
data_f$c_name_f <- factor(data_f$c_name)

r2.n.o.l.slm <- lm(regime_lag~poverty_nat+I(poverty_nat^2)+log(GDP_pc_ppp)+
                     log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                   data=data_f)


ggeffect(r2.n.o.l.slm, terms = "poverty_nat")
mydf <- ggeffect(r2.n.o.l.slm, terms = "poverty_nat")
#mydf1 <- ggpredict(r2.n.o.l.slm, terms = "poverty19")
r2n<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_nat")+
  ggtitle("Marginal effect of poverty_nat")+
  theme_light()


r2.19.o.l <- plm(regime_lag~log(poverty19+1)+I((poverty19+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata)

summary(r2.19.o.l)
coeftest(r2.19.o.l, vcov = vcovHC, type = "HC3")

r2.19.o.l.slm <- lm(regime_lag~log(poverty19+1)+I(log(poverty19+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(r2.19.o.l.slm, terms = "poverty19")
mydf <- ggeffect(r2.19.o.l.slm, terms = "poverty19")
#mydf1 <- ggpredict(r2.n.o.l.slm, terms = "poverty19")
r219<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty19")+
  ggtitle("Marginal effect of poverty19")+
  theme_light()

r2.32.o.l <- plm(regime_lag~log(poverty32+1)+I(log(poverty32+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata) 

summary(r2.32.o.l)
coeftest(r2.32.o.l, vcov = vcovHC, type = "HC3")

r2.32.o.l.slm <- lm(regime_lag~log(poverty32+1)+I(log(poverty32+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(r2.32.o.l.slm, terms = "poverty32")
mydf <- ggeffect(r2.32.o.l.slm, terms = "poverty32")
#mydf1 <- ggpredict(r2.n.o.l.slm, terms = "poverty19")
r232<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty32")+
  ggtitle("Marginal effect of poverty32")+
  theme_light()

r2.55.o.l <- plm(regime_lag~log(poverty55+1)+I(log(poverty55+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", data=paneldata) 

summary(r2.55.o.l)
coeftest(r2.55.o.l, vcov = vcovHC, type = "HC3")

r2.55.o.l.slm <- lm(regime_lag~log(poverty55+1)+I(log(poverty55+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                    data=data_f)


ggeffect(r2.55.o.l.slm, terms = "poverty55")
mydf <- ggeffect(r2.55.o.l.slm, terms = "poverty55")
#mydf1 <- ggpredict(r2.n.o.l.slm, terms = "poverty19")
r255<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty55")+
  ggtitle("Marginal effect of poverty55")+
  theme_light()

r2.19g.o.l <- plm(regime_lag~log(poverty_gap_19+1)+I(log(poverty_gap_19+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata)

summary(r2.19g.o.l)
coeftest(r2.19g.o.l, vcov = vcovHC, type = "HC3")

r2.19g.o.l.slm <- lm(regime_lag~log(poverty_gap_19+1)+I(log(poverty_gap_19+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                     data=data_f)


ggeffect(r2.19g.o.l.slm, terms = "poverty_gap_19")
mydf <- ggeffect(r2.19g.o.l.slm, terms = "poverty_gap_19")
#mydf1 <- ggpredict(r2.n.o.l.slm, terms = "poverty19")
r219g<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_19")+
  ggtitle("Marginal effect of poverty_gap_19")+
  theme_light()

r2.32g.o.l <- plm(regime_lag~log(poverty_gap_32+1)+I(log(poverty_gap_32+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata) 

summary(r2.32g.o.l)
coeftest(r2.32g.o.l, vcov = vcovHC, type = "HC3")

r2.32g.o.l.slm <- lm(regime_lag~log(poverty_gap_32+1)+I(log(poverty_gap_32+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                     data=data_f)


ggeffect(r2.32g.o.l.slm, terms = "poverty_gap_32")
mydf <- ggeffect(r2.32g.o.l.slm, terms = "poverty_gap_32")
#mydf1 <- ggpredict(r2.n.o.l.slm, terms = "poverty32")
r232g<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_32")+
  ggtitle("Marginal effect of poverty_gap_32")+
  theme_light()

r2.55g.o.l <- plm(regime_lag~log(poverty_gap_55+1)+I(log(poverty_gap_55+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", data=paneldata) 

summary(r2.55g.o.l)
coeftest(r2.19.o.l, vcov = vcovHC, type = "HC3")

r2.55g.o.l.slm <- lm(regime_lag~log(poverty_gap_55+1)+I(log(poverty_gap_55+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f,
                     data=data_f)


ggeffect(r2.55g.o.l.slm, terms = "poverty_gap_55")
mydf <- ggeffect(r2.55g.o.l.slm, terms = "poverty_gap_55")
#mydf1 <- ggpredict(r2.n.o.l.slm, terms = "poverty55")
r255g<-ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_55")+
  ggtitle("Marginal effect of poverty_gap_55")+
  theme_light()

figure2 <- ggarrange(r2n, r219, r232, r255, r219g, r232g, r255g,
                     ncol = 2, nrow = 4)
figure2


pFtest(r2.n.o.l, r2.n.p.l)
pFtest(r2.19.o.l, r2.19.p.l)
pFtest(r2.32.o.l, r2.32.p.l)
pFtest(r2.55.o.l, r2.55.p.l)
pFtest(r2.19g.o.l, r2.19g.p.l)
pFtest(r2.32g.o.l, r2.32g.p.l)
pFtest(r2.55g.o.l, r2.55g.p.l)

se2ol <- c(list(sqrt(diag(vcovHC(r2.n.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.19.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.32.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.55.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.19g.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.32g.o.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.55g.o.l, type = "HC1")))))

stargazer(r2.n.o.l, r2.19.o.l, r2.32.o.l, r2.55.o.l, r2.19g.o.l, r2.32g.o.l,
          r2.55g.o.l, se = se2ol)





####two-ways####

r2.n.t.l <- plm(regime_lag~poverty_nat+I(poverty_nat^2)+
                  log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                model="within", effect="twoways", data=paneldata)

summary(r2.n.t.l)
coeftest(r2.n.t.l, vcov = vcovHC, type = "HC3")

r2.19.t.l <- plm(regime_lag~log(poverty19+1)+I(log(poverty19+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata)

summary(r2.19.t.l)
coeftest(r2.19.t.l, vcov = vcovHC, type = "HC3")

r2.32.t.l <- plm(regime_lag~log(poverty32+1)+I(log(poverty32+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata) 

summary(r2.32.t.l)
coeftest(r2.32.t.l, vcov = vcovHC, type = "HC3")

r2.55.t.l <- plm(regime_lag~log(poverty55+1)+I(log(poverty55+1)^2)+
                   log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                 model="within", effect="twoways", data=paneldata) 

summary(r2.55.t.l)
coeftest(r2.55.t.l, vcov = vcovHC, type = "HC3")

r2.19g.t.l <- plm(regime_lag~log(poverty_gap_19+1)+I(log(poverty_gap_19+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata)

summary(r2.19g.t.l)
coeftest(r2.19g.t.l, vcov = vcovHC, type = "HC3")

r2.32g.t.l <- plm(regime_lag~log(poverty_gap_32+1)+I(log(poverty_gap_32+1)^2)+
                    log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata) 

summary(r2.32g.t.l)
coeftest(r2.32g.t.l, vcov = vcovHC, type = "HC3")

r2.55g.t.l <- plm(regime_lag~log(poverty_gap_55+1)+I(log(poverty_gap_55+1)^2)
                  +log(GDP_pc_ppp)+log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB), 
                  model="within", effect="twoways", data=paneldata) 

summary(r2.55g.t.l)
coeftest(r2.55g.t.l, vcov = vcovHC, type = "HC3")


se2tl <- c(list(sqrt(diag(vcovHC(r2.n.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.19.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.32.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.55.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.19g.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.32g.t.l, type = "HC1"))),
                sqrt(diag(vcovHC(r2.55g.t.l, type = "HC1")))))

stargazer(r2.n.t.l, r2.19.t.l, r2.32.t.l, r2.55.t.l, r2.19g.t.l, r2.32g.t.l,
          r2.55g.t.l, se = se2tl)

r2.n.t.l.slm <- lm(regime_lag~log(poverty_nat+1)+I(log(poverty_nat+1)^2)+log(GDP_pc_ppp)+
                     log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                   data=data_f)


ggeffect(r2.n.t.l.slm, terms = "poverty_nat")
mydf <- ggeffect(r2.n.t.l.slm, terms = "poverty_nat")
#mydf1 <- ggpredict(r2.n.t.l.slm, terms = "poverty55")
r2nt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_nat")+
  ggtitle("Marginal effect of poverty_nat")+
  theme_light()

r2.19.t.l.slm <- lm(regime_lag~log(poverty19+1)+I(log(poverty19+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(r2.19.t.l.slm, terms = "poverty19")
mydf <- ggeffect(r2.19.t.l.slm, terms = "poverty19")
#mydf1 <- ggpredict(r2.n.t.l.slm, terms = "poverty55")
r219t <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty19")+
  ggtitle("Marginal effect of poverty19")+
  theme_light()

r2.32.t.l.slm <- lm(regime_lag~log(poverty32+1)+I(log(poverty32+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(r2.32.t.l.slm, terms = "poverty32")
mydf <- ggeffect(r2.32.t.l.slm, terms = "poverty32")
#mydf1 <- ggpredict(r2.n.t.l.slm, terms = "poverty55")
r232t <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty32")+
  ggtitle("Marginal effect of poverty32")+
  theme_light()

r2.55.t.l.slm <- lm(regime_lag~log(poverty55+1)+I(log(poverty55+1)^2)+log(GDP_pc_ppp)+
                      log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+c_name_f+year_f,
                    data=data_f)


ggeffect(r2.55.t.l.slm, terms = "poverty55")
mydf <- ggeffect(r2.55.t.l.slm, terms = "poverty55")
#mydf1 <- ggpredict(r2.n.t.l.slm, terms = "poverty55")
r255t <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty55")+
  ggtitle("Marginal effect of poverty55")+
  theme_light()

r2.19g.t.l.slm <- lm(regime_lag~log(poverty_gap_19+1)+I(log(poverty_gap_19+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)
                     +c_name_f+year_f,
                     data=data_f)


ggeffect(r2.19g.t.l.slm, terms = "poverty_gap_19")
mydf <- ggeffect(r2.19g.t.l.slm, terms = "poverty_gap_19")
#mydf1 <- ggpredict(r2.n.t.l.slm, terms = "poverty19g")
r219gt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_19")+
  ggtitle("Marginal effect of poverty_gap_19")+
  theme_light()

r2.32g.t.l.slm <- lm(regime_lag~log(poverty_gap_32+1)+I(log(poverty_gap_32+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)+
                       c_name_f+year_f,
                     data=data_f)


ggeffect(r2.32g.t.l.slm, terms = "poverty_gap_32")
mydf <- ggeffect(r2.32g.t.l.slm, terms = "poverty_gap_32")
#mydf1 <- ggpredict(r2.n.t.l.slm, terms = "poverty32g")
r232gt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_32")+
  ggtitle("Marginal effect of poverty_gap_32")+
  theme_light()

r2.55g.t.l.slm <- lm(regime_lag~log(poverty_gap_55+1)+I(log(poverty_gap_55+1)^2)+log(GDP_pc_ppp)+
                       log(population_WB)+v2regsupgroupssize+v2xeg_eqdr+log(gini_WB)
                     +c_name_f+year_f,
                     data=data_f)


ggeffect(r2.55g.t.l.slm, terms = "poverty_gap_55")
mydf <- ggeffect(r2.55g.t.l.slm, terms = "poverty_gap_55")
#mydf1 <- ggpredict(r2.n.t.l.slm, terms = "poverty55g")
r255gt <- ggplot(mydf, aes(x, predicted)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  ylab("")+
  xlab("poverty_gap_55")+
  ggtitle("Marginal effect of poverty_gap_55")+
  theme_light()

figure2t <- ggarrange(r2nt, r219t, r232t, r255t, r219gt, r232gt, r255gt,
                      ncol = 2, nrow = 4)
figure2t


pFtest(r2.19.t.l, r.19.t.l)
pFtest(r2.19g.t.l, r.19g.t.l)











