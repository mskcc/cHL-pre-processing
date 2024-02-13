#
# https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/
#
stripGlmLR = function(cm) {
  cm$y = c()
  cm$model = c()

  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()


  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()

  cm
}


#https://www.r-bloggers.com/2015/04/how-and-why-to-return-functions-in-r/

wrapGLMModel <- function(model) {
  force(model)
  ret=list(
    model=model,
    predict=function(newd) {
      predict(model,newdata=newd,type='response')
    }
  )
}

logisticFitter <- function(vars,yTarget,data) {
  formula <- paste(yTarget, paste(vars, collapse=' + '), sep=' ~ ')
  model <- glm(as.formula(formula), data, family=binomial(link='logit'))
#  model <- stripGLMModel(model)
  wrapGLMModel(model)
}
