#https://adv-r.hadley.nz/function-factories.html

#https://www.r-bloggers.com/2015/04/how-and-why-to-return-functions-in-r/


#stripGLMModel <- function(model) { ... ; model }
#https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/

wrapGLMModel <- function(model) {
  force(model)
  function(newd) {
    predict(model,newdata=newd,type='response')
  }
}

logisticFitter <- function(vars,yTarget,data) {
  formula <- paste(yTarget,
    paste(vars,collapse=' + '),sep=' ~ ')
  model <- glm(as.formula(formula),data,
               family=binomial(link='logit'))
#  model <- stripGLMModel(model)
  wrapGLMModel(model)
}

# Theory:
#   https://www.countbayesie.com/blog/2019/6/12/logistic-regression-from-bayes-theorem

#http://www.sthda.com/english/articles/36-classification-methods-essentials/151-logistic-regression-essentials-in-r/

#other refs:
# https://www.r-bloggers.com/2015/09/how-to-perform-a-logistic-regression-in-r/
# https://www.tidymodels.org/start/case-study/
# https://parsnip.tidymodels.org/index.html
#   https://parsnip.tidymodels.org/reference/logistic_reg.html
# https://stateofther.github.io/finistR2019/s-tidymodels.html

library(tidyverse)
library(caret)
# Load the data and remove NAs
data("PimaIndiansDiabetes2", package = "mlbench")
PimaIndiansDiabetes2 <- na.omit(PimaIndiansDiabetes2)
# Inspect the data
sample_n(PimaIndiansDiabetes2, 3)
# Split the data into training and test set
set.seed(123)
training.samples <- PimaIndiansDiabetes2$diabetes %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- PimaIndiansDiabetes2[training.samples, ]
test.data <- PimaIndiansDiabetes2[-training.samples, ]

#
# for binomial link='logit' is default but be explicit
#

model <- glm( diabetes ~., data = train.data, family = binomial)
model <- glm( diabetes ~ glucose, data = train.data, family = binomial(link='logit'))

mm=logisticFitter(c("glucose"),"diabetes",train.data)
mm(test.data)
