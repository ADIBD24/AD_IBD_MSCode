# Define a function named prediction_3 that takes training data and a name as arguments
prediction_3 <- function(train_a, name) {
  # Set a random seed for reproducibility
  set.seed(123) 
  
  # Create a data partition for training and testing with a 75% training split
  train_index <- createDataPartition(train_a$AD, p = 0.75, list = FALSE)
  train_data <- train_a[train_index, ]
  test_data <- train_a[-train_index, ]
  
  # Fit a generalized linear model (GLM) with a binomial link
  a <- glm(AD ~ ., data = train_data, family = binomial(link = 'logit'))
  # Perform stepwise regression and print the summary of the model
  b <- stats::step(a)
  summary(b)
  
  # Perform an analysis of variance (ANOVA) test on the model and save the results
  argPlot <- anova(b, test = "Chisq")
  write.csv(argPlot, paste0(paste(title, name, sep = "-"), "-glmimportance.csv"))
  
  # Predict probabilities on the test set and create a ROCR prediction object
  prob <- predict(b, newdata = test_data, type = "response")
  pred <- ROCR::prediction(prob, test_data$AD)
  # Print the AUC (Area Under the Curve) value for the GLM model
  cat(name, " GLM performance AUC:", performance(pred, 'auc')@y.values[[1]], "\n")
  
  # Calculate the ROC curve and AUC using the pROC package
  roc3 <- pROC::roc(test_data$AD, prob)
  print(roc3$auc)
  
  # Plot the ROC curve using ggplot2 and save the plot as a PDF
  pdf(paste0(paste(title, name, sep = "-"), "-ROC.pdf"))
  # Combine ROC curves from different models for comparison and customize the plot
  ggroc(list("GLM AUC: 0.969" = roc3, "RF AUC: 0.842" = roc2, "CART AUC: 0.807" = roc1),
        aes = c("linetype", "color"), size = 1.0, alpha = .7, legacy.axes = TRUE) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), colour = "darkgrey", linetype = "dashed") +
    scale_colour_brewer(palette = "Dark2") + theme_bw() +
    ggtitle(paste0(paste(title, name, sep = "-"), " ROC curve")) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.75, 0.15))
  dev.off()
  
  # Save the ROC objects as an R data file
  save(roc1, roc2, roc3, file = paste0(paste(title, name, sep = "-"), "-ROC.Rdata"))
  
  # Output the confusion matrix for the GLM model
  perf <- performance(pred, 'acc', cutoffs = seq(0, 1, 0.01))
  cutoff <- perf@x.values[[1]][which.max(perf@y.values[[1]])]
  
  test_data$pred1 <- ifelse(prob > cutoff, "nonAD", "AD")
  c <- confusionMatrix(as.factor(test_data$pred1), as.factor(test_data$AD))
  print(c$table)
  cat(title, name, " glm acc: ", c$overall[1], "\n")
  
  # Perform 10-fold cross-validation for the GLM model
  set.seed(200)
  folds <- createFolds(y=train_a[,"AD"],k=10)
  
  max=0
  max_auc = 0 
  num=0
  num_auc=0
  
  for (i in 1:10) {
    fold_test <- train_a[folds[[i]],]
    fold_train <- train_a[-folds[[i]],]
    
    b1 <- glm(b$formula,data = fold_train,family = "binomial")
    fold_predict <- predict(b1,type = "response",newdata = fold_test)
    
    pred <- ROCR::prediction(fold_predict,fold_test$AD) 
    auc = performance(pred,'auc')@y.values[[1]] #AUC值
    
    fold_predict <- ifelse(fold_predict>cutoff,"nonAD","AD")
    fold_test$predict <- fold_predict
    
    c1 <- confusionMatrix(factor(fold_test$predict),fold_test$AD)
    fold_accuracy <- c1$overall[1]
    
    #fold_error <- as.numeric(fold_test[,15])-as.numeric(fold_test[,16])
    #fold_accuracy <- (nrow(fold_test)-sum(abs(fold_error)))/nrow(fold_test)
    
    if(fold_accuracy>max)
    {
      max=fold_accuracy  
      num=i
    }
    if(auc>max_auc)
    {
      max_auc=auc  
      num_auc=i
    }
  }
  
  cat(title, "GLM max ACC:",max," ",num,"\n")
  fold_test <- train_a[folds[[num]],]
  fold_train <- train_a[-folds[[num]],]
  
  b1 <- glm(b$formula,data = fold_train,family = "binomial")
  fold_predict <- predict(b1,type = "response",newdata = fold_test)
  
  fold_predict <- ifelse(fold_predict>cutoff,"nonAD","AD")
  fold_test$predict <- fold_predict
  
  c1 <- confusionMatrix(factor(fold_test$predict),fold_test$AD)
  print(c1$table)
  
  cat(title, "GLM max AUC:", max_auc," ",num_auc,"\n")
  
  fold_test <- train_a[folds[[num_auc]],]
  fold_train <- train_a[-folds[[num_auc]],]
  b1 <- glm(b$formula,data = fold_train,family = "binomial")
  fold_predict <- predict(b1,type = "response",newdata = fold_test)
  pred <- ROCR::prediction(fold_predict,fold_test$AD)
  performance(pred,'auc')@y.values 
  
  roc3<-pROC::roc(fold_test$AD, as.numeric(fold_predict)) 
  roc3$auc
  
  # Perform Random Forest analysis and plot the ROC curve
  # Print the title for the analysis
  cat(title, "Perform RF: ", "\n")
  
  # Set a seed for reproducibility
  set.seed(123)
  
  # Split the data into training and testing sets (75% training, 25% testing)
  train_index <- createDataPartition(train_a$AD, p=0.75, list=F)
  train_data <- train_a[train_index,]  # Training data
  train_data_group <- train_a$AD[train_index]  # Group/Class labels for training data
  
  test_data <- train_a[-train_index,]  # Testing data
  test_data_group <- train_a$AD[-train_index]  # Group/Class labels for testing data
  
  # Print the dimensions of the training data
  dim(train_data)
  
  # Convert the AD column to a factor for both training and testing data
  train_data$AD = as.factor(train_data$AD)
  test_data$AD = as.factor(test_data$AD)
  
  # Set up the control parameters for the random forest model
  train.control <- trainControl(method = "repeatedcv", number = 10, repeats = 10, classProbs=T, summaryFunction=twoClassSummary)
  
  # Train the random forest model
  ADNI_randomforest <- randomForest(AD ~ .,  # Formula specifying the response variable and predictors
                                    data = train_data,  # Training data
                                    ntree = 500,  # Number of trees to grow
                                    mtry = 3,  # Number of variables randomly sampled as candidates at each split
                                    importance = TRUE,  # Calculate variable importance
                                    proximity = TRUE)  # Calculate proximity measure
  
  # Obtain variable importance from the random forest model
  argPlot <- importance(ADNI_randomforest)
  
  
  # Commented out: Plot variable importance
  # varImpPlot(ADNI_randomforest, main = "variable importance")
  
  # Save the variable importance results to a CSV file
  write.csv(argPlot, paste0(paste(title, name, sep = "-"), "-rfimportance.csv"))
  
  # Predict on the testing data
  pre_ran <- predict(ADNI_randomforest, newdata = test_data)
  
  # Create a data frame combining the observed and predicted probabilities
  obs_p_ran = data.frame(prob = pre_ran, obs = test_data$AD)
  
  # Create a contingency table comparing the true values and predicted values
  table(test_data$AD, pre_ran, dnn = c("True", "Predict"))
  
  # Perform CART (Classification and Regression Trees) analysis and plot the ROC curve
  #set.seed(300)  #GWAS
  train_index <- createDataPartition(train_a$AD, p=0.75, list=F)
  train_data <- train_a[train_index,]
  test_data <- train_a[-train_index,]
  
  tree<-rpart(AD ~.,method = "class",data = train_data)
  argPlot <- as.data.frame(tree$variable.importance)
  write.csv(argPlot,paste0(paste(title,name,sep = "-"),"-rpartimportance.csv"))
  # pdf("plotcp.pdf")
  # plotcp(tree)
  # dev.off()
  #tree <- prune(tree, cp=0.031)
  pdf(paste0(paste(title,name,sep = "-"),"-rpartplot.pdf"))
  rpart.plot(tree)
  dev.off()
  
  test_data$pred<-predict(tree,test_data,type = "class")
  c <- confusionMatrix(test_data$pred,test_data$AD)
  print(c$table)
  pred <- ROCR::prediction(as.numeric(test_data$pred),test_data$AD)
  auc <- performance(pred,'auc')@y.values[[1]] #AUC值
  
  #https://blog.csdn.net/dege857/article/details/126599513
  roc1<-pROC::roc(as.ordered(test_data$AD) ,as.ordered(test_data$pred)) #GWAS / DEG
  roc1$auc
  # pdf(paste(title,"-CART-ROC.pdf"))
  # ggroc(roc1)
  # dev.off()
  
  cat(name," CART ACC: ", c$overall[1],"\n")
  cat(name," CART AUC: ", auc,"\n")
  
  cat("Perform CART 10cv: ","\n")
  
 ##10CV
  max=0
  max_auc = 0 
  num=0
  num_auc=0
  set.seed(200)  #adni set.seed(1)
  folds <- createFolds(y=train_a[,"AD"],k=10)
  
  for (i in 1:10) {
    fold_test <- train_a[folds[[i]],]
    fold_train <- train_a[-folds[[i]],]
    
    b1 <-rpart(AD~.,method = "class",data =  fold_train)
    fold_predict <- predict(b1,fold_test, type = "class")
    
    fold_test$predict <- fold_predict
    
    c1 <- confusionMatrix(factor(fold_test$predict),factor(fold_test$AD))
    fold_accuracy <- c1$overall[1]
    
    pred <- ROCR::prediction(as.numeric(fold_test$predict),fold_test$AD)
    auc <- performance(pred,'auc')@y.values[[1]] #AUC值
    #fold_error <- as.numeric(fold_test[,15])-as.numeric(fold_test[,16])
    #fold_accuracy <- (nrow(fold_test)-sum(abs(fold_error)))/nrow(fold_test)
    
    if(fold_accuracy>max)
    {
      max=fold_accuracy  
      num=i
    }
    if(auc>max_auc)
    {
      max_auc=auc  
      num_auc=i
    }
  }
  cat(name,title, "CART max ACC:",max," ",num,"\n")
  
  fold_test <- train_a[folds[[num]],]
  fold_train <- train_a[-folds[[num]],]
  
  b1 <-rpart(AD~.,method = "class",data =  fold_train)
  fold_predict <- predict(b1,fold_test, type = "class")
  
  fold_test$predict <- fold_predict
  
  c1 <- confusionMatrix(factor(fold_test$predict),factor(fold_test$AD))
  print(c1$table)
  
  cat(title, "CART max AUC:", max_auc," ",num_auc,"\n")
  fold_test <- train_a[folds[[num_auc]],]
  fold_train <- train_a[-folds[[num_auc]],]
  
  b1 <-rpart(AD~.,method = "class",data =  fold_train)
  fold_predict <- predict(b1,fold_test, type = "class")
  
  fold_test$predict <- fold_predict
  
  pred <- ROCR::prediction(as.numeric(fold_test$predict),fold_test$AD)
  auc <- performance(pred,'auc')@y.values[[1]] #AUC值
  roc1<-pROC::roc(as.ordered(fold_test$AD) ,as.ordered(fold_test$predict)) #ADIBD-KMEALL-DEG-MT 0.582/ ADIBD-KMEALL-DEG 0.5697/ADIBD-KMEALLDEG-MT+ALL 0.6667
  roc1$auc

  # Output the script execution process to a text file
  sink(paste0(paste(title, name, sep = "-"), "-output.txt"))
  
  # Execute the prediction_3 function
  prediction_3(train_a, name)
  
  # Close the sink to stop writing to the text file
  sink()
  
  # Print the completion message
  cat("PROCESS:", name, title, "done!!", "\n")
}

# Print the dimensions of the training data
dim(train_a)

# Print the start message
cat("PROCESS:", name, title, "begin!!", "\n")

# Call the prediction_3 function with the training data and the name
prediction_3(train_a, name)

# Print the end message
cat("PROCESS:", name, title, "done!!", "\n")
