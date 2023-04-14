get.best.model_C3 = function(bncmark){
  modelo = getBMRModels(bncmark)
  modelodf = as.data.frame(bncmark)
  
  # Check the model with best results(we asume here that we have only 1 algorithm)
  # We can change this according to algorithms used
  models_index = modelodf[
    order( modelodf[,4], modelodf[,5]),
  ]
  index = as.numeric(tail(rownames(models_index), 1))
  # Get the model
  best = getLearnerModel(modelo$dataset[[1]][[index]])
  return(best)
}

bm3 = get.best.model_C3(bncmark = bmr_3C)
test_task3 = makeClassifTask(data = Sriniv.predict_3C, target = "cluster")
test_task3 = normalizeFeatures(
  test_task3,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")
predict3 = predict(bm3, task = test_task3, type = "prob")


