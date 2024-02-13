require(tidyverse)
obj=readRDS("../data/HodgkinsV10dev/rda/v5_Exclusion/H16_4923___HaloObj_v10.7__210714_Exclusions_.rda")
atlas=readRDS("../data/HodgkinsV10dev/cellAtlas_Hodgkins_v10.7__210714_b_CTDv____ba5cf549_20210718_195242.rda")

mA="CD3"
mB="CD20"

uuids=obj$geom.data %>% filter(!Exclude) %>% pull(UUID)

mm=obj$marker.data %>%
  filter(UUID %in% uuids) %>%
  filter(Marker %in% c(mA,mB)) %>%
  select(UUID,Marker,Pos=Positive_Classification,TIntensity) %>%
  mutate(Pos=factor(Pos)) %>%
  left_join(select(obj$geom.data,UUID,Sample,SPOT))

ggplot(mm,aes(Pos,TIntensity,color=Pos)) + geom_violin() + facet_wrap(~Marker)

yy=mm %>%
  mutate(Pos=as.numeric(Pos)-1) %>%
  gather(Metric,Value,Pos,TIntensity) %>%
  unite(MarkerMetric,Marker,Metric,sep=".") %>%
  spread(MarkerMetric,Value) %>%
  mutate(DPos=cc(CD3.Pos,CD20.Pos)) %>%
  left_join(atlas %>% filter(UUID %in% uuids) %>% select(UUID,CellType))

ggplot(yy,aes(CD3.TIntensity,CD20.TIntensity)) + geom_hex(bins=25) + facet_wrap(~DPos) + theme_light()

train=yy %>% filter(CellType!="UNKNOWN")
test=yy %>% filter(DPos=="1_1")

modelA=glm(CD3.Pos ~ CD3.TIntensity , data=train, family = binomial(link='logit'))
modelB=glm(CD20.Pos ~ CD20.TIntensity , data=train, family = binomial(link='logit'))

#yy %>% ggplot(aes(X,Y)) + geom_point() + geom_smooth(method = "glm",method.args = list(family = "binomial"))

modelA2=glm(CD3.Pos ~ CD3.TIntensity + CD20.TIntensity , data=train, family = binomial(link='logit'))
modelB2=glm(CD20.Pos ~ CD20.TIntensity + CD3.TIntensity, data=train, family = binomial(link='logit'))

pp=tibble(A=predict(modelA,newdat=test,type='response'),B=predict(modelB,newdat=test,type='response'))
pp=pp %>% mutate(lOA=log(A/(1-A)),lOB=log(B/(1-B)),lOR=lOA-lOB)

pp2=tibble(A=predict(modelA2,newdat=test,type='response'),B=predict(modelB2,newdat=test,type='response'))
pp2=pp2 %>% mutate(lOA=log(A/(1-A)),lOB=log(B/(1-B)),lOR=lOA-lOB)
