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
  
  # Perform 10-fold cross-validation for the GLM model and print the maximum accuracy and AUC values
  # ...
  # (The code for 10-fold cross-validation is omitted for brevity.)
  
  # Perform Random Forest analysis and plot the ROC curve
  # ...
  # (The code for Random Forest analysis is omitted for brevity.)
  
  # Perform CART (Classification and Regression Trees) analysis and plot the ROC curve
  # ...
  # (The code for CART analysis is omitted for brevity.)
  
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