# Select Build, Build and reload to build and lode into the R-session.

mylm <- function(formula, data = list(), contrasts = NULL, ...){
  # Extract model matrix & responses
  mf <- model.frame(formula = formula, data = data)
  X  <- model.matrix(attr(mf, "terms"), data = mf, contrasts.arg = contrasts)
  y  <- model.response(mf)
  terms <- attr(mf, "terms")


  # Add code here to calculate coefficients, residuals, values, etc...

  #a)
  #Finding coefficient
  Bhat = solve(t(X)%*%X)%*%t(X)%*%y

  #b)

  #First we need to find an estimate for sigma squared
  fitted_values = X %*% Bhat
  residuals=y-X%*%Bhat
  sigma_Squared=as.numeric(1/(length(y)-ncol(X)+1)*t(residuals)%*%residuals)

  #Using this to find covariance matrix of coefficients
  cov_Bhat=sigma_Squared*solve(t(X)%*%X)
  SE_Bhat=as.numeric(sqrt(diag(cov_Bhat)))

  #Finding z-values and p-values
  z_Bhat=rep(0,length(Bhat))
  p_Value=rep(0,length(Bhat))
  for (i in 1:length(Bhat)){
    z_Bhat[i]=Bhat[i]/SE_Bhat[i]
    p_Value[i]=2*pnorm(abs(z_Bhat[i]), lower.tail=FALSE)
  }

  #d)
  #Calculating the SSE and the degrees of freedom
  SSE <- sum((t(residuals)%*%residuals))
  n <- length(y)
  p <- length(Bhat)
  df <- (n - (p+1))

  #Calculating SST
  ressidual_H0 <- y-mean(y)
  SST <- sum(t(ressidual_H0)%*%ressidual_H0)


  #Testing the significance of the regression with a  χ2-test
  F_statistic <- ((df)*(SST-SSE))/(SSE*(p-1))

  F_pvalue <- 1 - pchisq(F_statistic*(p-1), df=(p-1))


  #Find the critical values for both tests #still assuming the significance level of 0.05
  critical_z <- qnorm(1- (0.05/2))
  critical_chi <- qchisq(1- 0.05,df)

  #Calculate R^2
  R_suared <- 1 - SSE/SST


  # and store the results in the list est
  est <- list(terms = terms,
              model = mf,
              coefficients=Bhat,
              fitted_values = fitted_values,
              residuals = residuals,
              SE=SE_Bhat,
              z=z_Bhat,
              p=p_Value,
              SSE=SSE,
              SST=SST,
              critical_z = critical_z,
              critical_chi = critical_chi,
              F_statistic=F_statistic,
              F_pvalue = F_pvalue,
              R_suared=R_suared,
              df = df)

  # Store call and formula used
  est$call <- match.call()
  est$formula <- formula

  # Set class name. This is very important!
  class(est) <- 'mylm'

  # Return the object with all results
  return(est)
}

print.mylm <- function(object, ...){
  # Code here is used when print(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Info about object\n\n')
  cat('Coefficients\n')
  print(round(object$coefficients,3))
}

summary.mylm <- function(object, ...){
  # Code here is used when summary(object) is used on objects of class "mylm"
  # Useful functions include cat, print.default and format
  cat('Summary of object\n\n')

  #Creating summary table
  summary_Table=data.frame(matrix(c(round(object$coefficients, 5),
                                    round(object$SE, 5),
                                    round(object$z, 5),
                                    object$p),
                                  nrow=length(object$coefficients)))
  colnames(summary_Table) <- c('Estimate', 'Std. Error', 'z-value', 'p-value')
  rownames(summary_Table) <- rownames(object$coefficients)

  summary_F=data.frame(matrix(c(object$F_statistic, object$F_pvalue),nrow = 1,byrow = TRUE))
  colnames(summary_F) <- c('F-statistic','p-value (chi square)')

  print(summary_Table)
  print(summary_F)

#  if (F_pvalue < 0.05) {
#    print("The regression is significant")
#  } else {
#    print("The regression is not significant")
#  }
  cat("R-squared:")
  print(object$R_suared)


}

plot.mylm <- function(object, ...){
  # Code here is used when plot(object) is used on objects of class "mylm"

  #c)
  #Creating a plot function
  plot(object$fitted_values, object$residuals, xlab="Fitted Values", ylab="Raw Residuals")
  abline(a=0, b=0, col="blue", lty=3)
}



# This part is optional! You do not have to implement anova
anova.mylm <- function(object, ...){
  # Code here is used when anova(object) is used on objects of class "mylm"

  # Components to test
  comp <- attr(object$terms, "term.labels")

  # Name of response
  response <- deparse(object$terms[[2]])

  # Fit the sequence of models
  txtFormula <- paste(response, "~", sep = "")
  model <- list()
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula <- paste(txtFormula, comp[numComp])
    }
    else{
      txtFormula <- paste(txtFormula, comp[numComp], sep = "+")
    }
    formula <- formula(txtFormula)
    model[[numComp]] <- lm(formula = formula, data = object$model)
  }

  # Print Analysis of Variance Table
  cat('Analysis of Variance Table\n')
  cat(c('Response: ', response, '\n'), sep = '')
  cat('          Df  Sum sq X2 value Pr(>X2)\n')
  for(numComp in 1:length(comp)){
    # Add code to print the line for each model tested
  }

  return(model)

}
