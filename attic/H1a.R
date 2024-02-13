_HiStOrY_V2_
source("logisticModel.R")
yy
model=glm(CD3.Pos ~ CD3.TIntensity + CD20.TIntensity, data=train)
model
summary(model)
modelA=glm(CD3.Pos ~ CD3.TIntensity + CD20.TIntensity, data=train)
modelB=glm(CD20.Pos ~ CD3.TIntensity + CD20.TIntensity, data=train)
modelA
modelB
modelA
yy %>% filter(DPos!="1_1")
yy %>% filter(DPos=="1_1")
test=yy %>% filter(DPos=="1_1")
predict(model,dat=test)
predict(modelA,dat=test)
cbind(predict(modelA,dat=test),predict(modelB,dat=test))
cbind(A=predict(modelA,dat=test),B=predict(modelB,dat=test)) %>% tibble
modelA=glm(CD3.Pos ~ CD3.TIntensity, data=train)
modelB=glm(CD20.Pos ~ CD20.TIntensity, data=train)
cbind(A=predict(modelA,dat=test),B=predict(modelB,dat=test)) %>% tibble
cbind(A=predict(modelA,dat=test),B=predict(modelB,dat=test)) %>% as_tibble
pp=cbind(A=predict(modelA,dat=test),B=predict(modelB,dat=test)) %>% as_tibble
ggplot(pp,aes(A,B)) + geom_point()
model
modelA=glm(CD3.Pos ~ CD3.TIntensity + CD20.TIntensity, data=train, family = binomial(link='logit')))
modelA=glm(CD3.Pos ~ CD3.TIntensity + CD20.TIntensity, data=train, family = binomial(link='logit'))
train
modelA=glm(CD3.Pos ~ CD3.TIntensity, data=train, family = binomial(link='logit'))
yy
yy %>% ggplot(aes(CD3.TIntensity,CD3.Pos))
yy %>% ggplot(aes(CD3.TIntensity,CD3.Pos)) + geom_point()
yy %>% ggplot(aes(CD3.TIntensity,CD3.Pos)) + geom_point() + geom_smooth(method = "glm",method.args = list(family = "binomial"))
modelA=glm(CD3.Pos ~ CD3.TIntensity , data=train, family = binomial(link='logit'))
modelA
predict(modelA,dat=test)
modelA=glm(CD3.Pos ~ CD3.TIntensity , data=train, family = binomial(link='logit'))
yy %>% ggplot(aes(CD3.TIntensity,CD3.Pos)) + geom_point() + geom_smooth(method = "glm",method.args = list(family = "binomial"))
test
test %>% select(CD3.TIntensity)
tt=test %>% select(CD3.TIntensity)
modelA=glm(CD3.Pos ~ CD3.TIntensity , data=train, family = binomial(link='logit'))
predict(modelA,newdat=test,type='response')
tibble(A=predict(modelA,newdat=test,type='response'))
modelB=glm(CD20.Pos ~ CD20.TIntensity , data=train, family = binomial(link='logit'))
tibble(A=predict(modelA,newdat=test,type='response'),B=predict(modelB,newdat=test,type='response'))
pp=tibble(A=predict(modelA,newdat=test,type='response'),B=predict(modelB,newdat=test,type='response'))
ggplot(pp,aes(A,B)) + geom_point()
ggplot(pp,aes(A/(1-A),B)) + geom_point()
ggplot(pp,aes(A/(1-A),B/(1-B))) + geom_point()
ggplot(pp,aes(A/(1-A),B/(1-B))) + geom_point() + scale_x_continuous(trans="log")
ggplot(pp,aes(A/(1-A),B/(1-B))) + geom_point() + scale_x_continuous(trans="log") + scale_y_continuous(trans="log")
savehistory("H1.R")
