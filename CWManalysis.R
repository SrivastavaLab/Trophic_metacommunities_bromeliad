library(tidyverse)
library(car)
library(visreg)
library(lavaan)
library(semPlot)
library(MASS)

depth<-read.csv("Final_Water_Depth.csv", header=TRUE)
catchment<-read.csv("catchment_areas.csv", header=TRUE)

pit<-read.csv("pitilla_all_years.csv", header=TRUE)
pit1<-filter(pit, max_water!="NA")%>%mutate(culex=Culex_daumastocampa+Culex_daumastocampa_or_jenningsi_or_rejector+Culex_jenningsi,
                                            wyeomyia=Wyeomyia_abebela+Wyeomyia_abebela_or_ccircumcincta_or_melanopus+Wyeomyia_melanopus,
                                            mecist=Mecistogaster_modesta,
                                            culex_pres=1*(culex>0),
                                            wyeomyia_pres=1*(wyeomyia>0),
                                            mecist_pres=1*(mecist>0),
                                            prop.culex=culex/total_number_individuals,
                                            prop.wyeomyia=wyeomyia/total_number_individuals,
                                            prop.mecist=mecist/total_number_individuals,
                                            culex.null=total_number_individuals*0.01607057,
                                            wyeomyia.null=total_number_individuals*0.02532862,
                                            mecist.null=total_number_individuals*0.01545919,
                                            culex.resid=culex-total_number_individuals*0.01607057,
                                            wyeomyia.resid=wyeomyia-total_number_individuals*0.02532862,
                                            mecist.resid=mecist-total_number_individuals*0.01545919,
                                            culex.sign=sign(culex.resid),
                                            wyeomyia.sign=sign(wyeomyia.resid),
                                            mecist.sign=sign(mecist),
                                            culex.lresid=log(abs(culex.resid))*culex.sign,
                                            wyeomyia.lresid=log(abs(wyeomyia.resid))*wyeomyia.sign,
                                            mecist.lresid=log(abs(mecist.resid))*mecist.sign)


pitlong.raw<-select(pit1, max_water, culex, wyeomyia, mecist)%>%gather(species, abund, -max_water)
pitlong.null<-select(pit1, max_water, culex.null, wyeomyia.null, mecist.null)%>%gather(species, null.abund, -max_water)
pitlong.raw$null.abund<-pitlong.null$null.abund
pitlong.prop<-select(pit1, max_water, prop.culex, prop.wyeomyia, prop.mecist)%>%gather(species, abund, -max_water)
pitlong.occ<-select(pit1, max_water, culex_pres, wyeomyia_pres, mecist_pres)%>%gather(species, abund, -max_water)
pitlong.resid<-select(pit1, max_water, culex.resid, wyeomyia.resid, mecist.resid)%>%gather(species, abund, -max_water)
pitlong.lresid<-select(pit1, max_water, culex.lresid, wyeomyia.lresid, mecist.lresid)%>%gather(species, abund, -max_water)

#panel a

p <- ggplot(pitlong.raw, aes(log(max_water), (abund)^0.5, colour=species))
p + geom_point()+geom_smooth(se = FALSE, method = "glm", method.args = list(family = "poisson"))

p <- ggplot(pit1, aes(log(max_water), (total_number_individuals*0.015)^0.5))
p + geom_point()+geom_smooth(se = FALSE, method = "glm", method.args = list(family = "poisson"))

a1<-glm((total_number_individuals)^0.5~log(max_water), family=poisson, data =pit1)
summary(a1)
par(mfrow=c(2,2))
plot(a1)
a2<-glm((culex)^0.5~log_max_water, data= pit1, family=poisson); anova(a2, test="Chi")
a2<-glm((mecist)^0.5~log_max_water, data= pit1, family=poisson); anova(a2, test="Chi")
a2<-glm((wyeomyia)^0.5~log_max_water, data= pit1, family=poisson); anova(a2, test="Chi")

cor.test((pit1$culex)^0.5,(pit1$culex.null)^0.5)
cor.test((pit1$lculex.resid)^0.5,log(pit1$max_water))
cor.test((pit1$mecist)^0.5,(pit1$mecist.null)^0.5)
cor.test((pit1$lmecist.resid)^0.5,log(pit1$max_water))
cor.test((pit1$wyeomyia)^0.5,(pit1$wyeomyia.null)^0.5)
cor.test((pit1$lwyeomyia.resid)^0.5,log(pit1$max_water))

a2<-glm((culex)^0.5~I(culex.null^0.5), family=poisson, data =pit1)
p <- ggplot(pit1, aes((culex.null)^0.5, (culex)^0.5))
p + geom_point()+geom_smooth(se = FALSE, method = "glm", method.args = list(family = "poisson"))



sum(pit1$culex)/sum(pit1$total_number_individuals) #0.01607057
sum(pit1$wyeomyia)/sum(pit1$total_number_individuals) #0.02532862
sum(pit1$mecist)/sum(pit1$total_number_individuals) #sum(pit1$wyeomyia)/sum(pit1$total_number_individuals)

p <- ggplot(pitlong.resid, aes(log(max_water), log(abund+35), colour=species))
p + geom_point()+geom_smooth(se = FALSE, method = "glm", method.args = list(family = "poisson"))

p <- ggplot(pitlong.lresid, aes(log(max_water), abund, colour=species))
p + geom_point()+geom_smooth(se=FALSE)

pit1$log_max_water<-log(pit1$max_water)

c3<-glm(culex~culex.null, data= pit1, family=gaussian)
c2<-glm(culex.lresid~log_max_water, data= pit1, family=gaussian)
c1<-glm(culex.lresid~log_max_water+I(log_max_water^2)+I((wyeomyia)^0.5)+I((mecist)^0.5), data= pit1, family=gaussian)
c0<-glm(culex.lresid~I((wyeomyia)^0.5)+I((mecist)^0.5), data= pit1, family=gaussian)
par(mfrow=c(2,2)); plot(c1)
Anova(c1, type=2)
anova(c0,c1, test="Chi")

m3<-glm(mecist~mecist.null, data= pit1, family=gaussian)
m2<-glm(mecist.lresid~log_max_water, data= pit1, family=gaussian)
m1<-glm(mecist.lresid~log_max_water+I(log_max_water^2), data= pit1, family=gaussian)
par(mfrow=c(2,2)); plot(m1)
Anova(m1, type=2)
m0<-glm(mecist.lresid~1, data= pit1, family=gaussian)
anova(m0,m1, test="Chi")
par(mfrow=c(1,1));visreg(m1,  ylab="Mecistogaster residuals, log scale")

w3<-glm(wyeomyia~wyeomyia.null, data= pit1, family=gaussian)
w2<-glm(wyeomyia.lresid~log_max_water, data= pit1, family=gaussian)
w1<-glm(wyeomyia.lresid~log_max_water+I(log_max_water^2)+I((culex)^0.5)+I((mecist)^0.5), data= pit1, family=gaussian)
par(mfrow=c(2,2)); plot(w1)
w0<-glm(wyeomyia.lresid~I((culex)^0.5)+I((mecist)^0.5), data= pit1, family=gaussian)
Anova(w1, type=2)
anova(w0,w1, test="Chi")
w2<-glm(wyeomyia.lresid~log_max_water+I((culex)^0.5)+I((mecist)^0.5), data= pit1, family=gaussian)
Anova(w2, type=2)
summary(w1)
visreg(w1, ylab="Wyemoyia residuals, log scale")


#Here I considered putting the null assembly explanation in the model as a term
pit1$culex.null.disp<-(pit1$culex.null)/8.3
c2<-glm((culex/8.3)~culex.null.disp+log_max_water+I(log_max_water^2)+I((wyeomyia)^0.5)+I((mecist)^0.5), data= pit1, family=poisson)
par(mfrow=c(2,2)); plot(c2)
Anova(c2, type=2)
summary(c2)
visreg(c2)

m2<-glm(mecist~mecist.null+log_max_water+I(log_max_water^2), data= pit1, family=gaussian)
par(mfrow=c(2,2)); plot(m2)
Anova(m2, type=2)
summary(m2)
visreg(m2)

pit1$wyeomyia.null.disp<-(pit1$wyeomyia.null)/14.12125
w2<-glm(wyeomyia/14.12125~wyeomyia.null.disp+log_max_water+I(log_max_water^2)+I((culex)^0.5)+I((mecist)^0.5), data= pit1, family=poisson)
par(mfrow=c(2,2)); plot(w2)
Anova(w2, type=2)
anova(w2, test="Chi")
summary(w2)

par(mfrow=c(1,1))
plot(pit1$mecist~pit1$mecist.null)
plot(pit1$culex~pit1$culex.null)
plot(pit1$wyeomyia~pit1$wyeomyia.null)

###SEM
pit1$wyeomyia.scale<-pit1$wyeomyia/10
pit1$culex.scale<-pit1$culex/10
pit1$mecist.scale<-pit1$mecist/10

mod_random <- ' # direct effect
culex.scale~log_max_water
mecist.scale~log_max_water
wyeomyia.scale~log_max_water
wyeomyia.lresid~log_max_water
culex.lresid~log_max_water
mecist.lresid~log_max_water
'

fit <- sem(mod_random, data = pit1)
varTable(fit) #variance is fine

mod_wyeomyia1<- ' # direct effect
wyeomyia.lresid~c*log_max_water

# mediator
mecist.lresid ~ a * log_max_water
wyeomyia.lresid ~ b* mecist.lresid
culex.lresid~d*log_max_water+f*mecist.lresid
wyeomyia.lresid ~ e* culex.lresid

# indirect effect (a*b)
ab := a*b
de := d*e
afe :=a*f*e

# total effect
total := c + (a*b)+(d*e)+(a*f*e)
'

fit <- sem(mod_wyeomyia, data = pit1)
varTable(fit)
summary(fit)
semPaths(fit)

mod_wyeomyia2<- ' # direct effect
wyeomyia.lresid~c*log_max_water

# mediator
mecist.scale ~ a * log_max_water
wyeomyia.lresid ~ b* mecist.scale
culex.scale~d*log_max_water+f*mecist.scale
wyeomyia.lresid ~ e* culex.scale

# indirect effect (a*b)
ab := a*b
de := d*e
afe :=a*f*e

# total effect
total := c + (a*b)+(d*e)+(a*f*e)
'
fit2 <- sem(mod_wyeomyia2, data = pit1)
varTable(fit2)
summary(fit2)
standardizedsolution(fit2)
semPaths(fit2)

mod_culex1<- ' # direct effect
culex.lresid~c*log_max_water

# mediator
mecist.scale ~ a * log_max_water
culex.lresid ~ b* mecist.scale
wyeomyia.scale~d*log_max_water+f*mecist.scale
culex.lresid ~ e* wyeomyia.scale

# indirect effect (a*b)
ab := a*b
de := d*e
afe :=a*f*e

# total effect
total := c + (a*b)+(d*e)+(a*f*e)
'
fit3 <- sem(mod_culex1, data = pit1)
varTable(fit3)
summary(fit3)
standardizedsolution(fit3)
semPaths(fit3)

mod_mecist1<- ' # direct effect
mecist.lresid~c*log_max_water
'

fit4 <- sem(mod_mecist1, data = pit1)
varTable(fit4)
summary(fit4)
standardizedsolution(fit4)
semPaths(fit4)

####Drought risk
#41190- 41264 (oct 8 2012 to Dec 21 2012) or 41532 to 41556 ( sept 15 to oct 9) this is sept 15 to dec 22, so about 2 wks before surveys to 3 wks (22 days) after
#22 days for Wyemyia smithii larval development
out3wk<-depth%>%filter((Julian.date<=41264)|((Julian.date>=41532)&(Julian.date>=41556)))%>%group_by(ID)%>%summarise(threewk= sum(Average>5))%>%mutate(prop3wk=threewk/39)
#max is 39


#M. modesta adults in Pacific zone seen wk 7 to wk 22, Hedstrom Sahlen, which is day 49 to 154
#41275 is Jan 1 so adult 41324 to 41429
out9mo<-depth%>%filter((Julian.date<=41324)|(Julian.date>=41429))%>%group_by(ID)%>%summarise(ninemo= sum(Average>5))%>%mutate(prop9mo=ninemo/132)
out9mo2<-depth%>%filter((Julian.date<=41324)|(Julian.date>=41429))%>%group_by(ID)%>%summarise(ninemo= sum(Average>-1))

outall<-left_join(out3wk, out9mo)%>%left_join(catchment)
outall$d9mo<-(1-(outall$prop9mo)^259)
outall$d3wk<-(1-(outall$prop3wk)^22)
outall$log_max_water<-log(outall$Max.volume)

m9mo<-glm(d9mo~log_max_water, family=binomial, data = outall); anova(m9mo)
visreg.m9mo<-visreg(m9mo, scale="response", rug=FALSE)
m3wk<-glm(d3wk~log_max_water, family=binomial, data = outall)
visreg.m3wk<-visreg(m3wk, scale="response", rug=FALSE)

plot(outall$d9mo~log(outall$Max.volume))
plot(outall$d3wk~log(outall$Max.volume))

dose.p(m9mo, cf = c(1,2), p = 1:3/4) #LD50 log vol: 4.682794 se: 0.2367013
exp(4.682794+1.96*0.2367013) #171.8683
exp(4.682794-1.96*0.2367013) # 67.95593
dose.p(m3wk, cf = c(1,2), p = 1:3/4) #LD50 log vol: 3.524957 se: 0.4605310
exp(3.524957+1.96*0.4605310) #84
exp(3.524957-1.96*0.4605310) #14
#
CWManalysis_list <- list(
  pitlong.raw = pitlong.raw,
  outall = outall,
  m9mo = m9mo,
  m3wk = m3wk,
  visreg.m9mo = visreg.m9mo,
  visreg.m3wk = visreg.m3wk
)

saveRDS(CWManalysis_list, "CWManalysis.rds")
