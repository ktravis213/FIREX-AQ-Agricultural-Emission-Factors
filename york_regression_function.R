#2020-04-30 Hannah Halliday EPA 

#Formualtion of the York Regression from York 2004 and Cantrell 2008
 
# Inputs ------------------------------------------------------------------

#data = dataframe containing all variables
#x = xvariable name, as string
#y = yvariable name, as string
#variance_x = x variable uncertainty (name, string). This should be the precision; Uncertainty data should be a column in the data frame, with an uncertainty for each observation.  
#variance_y = y variable uncertainty (name, string). Same restrincts and the x uncertainty. 
#ncycles = How many iterations of the algorithm should be done before recording a value 
#error_correlation = T or F. If false, we assume no error correlation between x and y, and ri is ignored. If T, a pointwise named column of the error correlation values should be included. 
#ri = pointwise error correlation values, name, string. Leave out when error_correlation = F 
#track_convergence = T or F. If set to T, will break for loop when the change in the calculated slope is below  1e-15; the number of cycles undertaken before convergence is recorded. If set to F, or if the change is slope is always >  1e-15, the for loop will terminate after ncycles. 

# Outputs -----------------------------------------------------------------

#Slope = slope of regression (y_units/x_units)
#Intercept = Intercept of regression (y_units)
#r_sq = Coefficenct of Determination. Calcualted from cor(x,y)^2. 
#uncertainty_slope = slope uncertinaty value
#uncertainty_intercept = intercept uncertainty value
#S = weighted sum of deviations from the best-fit line (y direction only, form from Cantrell's paper)  
#gofit = Goodness of Fit metric. See York and Cantrell papers for interpretation of results
#converged_cycles = how many cycles of the algorithm ran before the change in the calculated slope was < 1e-15? 

# Notes -------------------------------------------------------------------

#Helpful note: the function can handle NA data with no problem, I only calculate the regression with complete cases in the dataframe. 
#However, if there are a lot of zeros in a one of the variables, the track_convergence should be set to F; a lot of zero values can break the convergence check at the end of the For loop. 
#X, y, variance_x, and variance_y should be named columns in a dataframe. 
#If you are assuming that the correlatuions in the errors is zero, set error_correlation = F and ignore ri.
#Otherwise, include a named column that includes the pointwise error correlation values. 

#Example set up of the function:
#york_regression(data = data_cut, x = "co2_ppm", y = "co_ppm", variance_x = "co2_precision_ppm", variance_y = "co_precision_ppm")
#Ri is assumed to be zero; will run until solution is converged or 100 cycles have elapsed; convergence of solution is tracked. 

# Bug Tracker -------------------------------------------------------------

#Bug Tracker:
#2020-05-11 - Changed Ri to a pointwise value in the regression.


# Function ----------------------------------------------------------------

york_regression <- function(data = datacheck,  x = "variable_1", y = "variable_2", variance_x = "variance1", variance_y = "variance2",
                            error_correlation = F, ri = "error_correlation_variable", ncycles = 100, track_convergence = T) 
{
  
  #Set up the dataset that we'll use for these correlations:
  if(error_correlation == F) #We are assuming that Error Correlation is Zero: 
  {
    data_function <- data[,c(x, y, variance_x, variance_y)]      #Call by name, set in the function
    names(data_function) <- c("x", "y", "x_variance", "y_variance" ) 
    #Create a column of zeros: 
    data_function$ri <- 0 
  }
  else #if we are giving the function pointwise correlation values: 
  {
    data_function <- data[,c(x, y, variance_x, variance_y, ri)]      #Call by name, set in the function
    names(data_function) <- c("x", "y", "x_variance", "y_variance", "ri" )  #then rename for simplicity
  }

  #only keep the complete cases, to prevent issues with the calculation x-bar and y-bar. Worked out on 11-28-2017 with the chinese dataset.
  data_function <- data_function[which(complete.cases(data_function)),]
  n_obs <- nrow(data_function) #How many pairwise obs do we have?

  #Calculate the intial guess for the slope:
  #print(data_function)
  m_initial <- as.numeric(lm(data = data_function, y~x)$coefficients[2])

  for (i in 1:ncycles)
  {
    #Set the m value to m_initial for the first iteration, then change it.
    ifelse(i == 1, mloop <- m_initial, mloop <- m_new)
    #print(mloop)

    #Calculate the weights from the variances:
    weightx <- 1/(data_function$x_variance^2)
    weighty <- 1/(data_function$y_variance^2)

    #Calculate the alpha val:
    alpha_val <- sqrt(weightx * weighty)

    #Calculate combined weighting:
    Wi <- (weightx * weighty)/(weightx + (mloop^2)*weighty - 2*mloop*alpha_val*data_function$ri)

    #Calculate the weighted bar values:
    x_bar <- sum(Wi*data_function$x, na.rm = T)/sum(Wi, na.rm = T)
    y_bar <- sum(Wi*data_function$y, na.rm = T)/sum(Wi, na.rm = T)
    
    #Calculate the perturbation of each value:
    Ui <- data_function$x - x_bar
    Vi <- data_function$y - y_bar

    #calculate the B value:
    beta_val <- Wi*((Ui/weighty) + (mloop*Vi/weightx) - (mloop*Ui + Vi)*(data_function$ri/alpha_val))
    #beta_val <- Wi*((Ui/weighty) + (m_initial*Vi/weightx))

    #Then re-calculate the m value:
    m_new <- sum(Wi*beta_val*Vi, na.rm = T)/sum(Wi*beta_val*Ui, na.rm = T)
    
    #Track how many cycles it takes for the solution to converge? 
    if(track_convergence == T)
    { 
      #Lets add a convergence test (2020-04)
      #print(mloop)
      if( abs((m_new - mloop)/mloop) < 1e-15)
      {
        break
      }
    }
    
  }

  m_calc <- m_new
  intercept_calc <- y_bar - m_calc*x_bar

  #################################################
  #Aright, now do the uncertainty calcs: (Cantrell method)
  beta_bar <- sum(Wi * beta_val)/sum(Wi) #Deviation of Beta

  #Edited on 2020-04-07 - implementing York's calculation of uncertainties:
  xi = x_bar + beta_val #Calculate adjusted x values (no need for y)
  small_x_bar = sum(Wi*xi, na.rm = T)/sum(Wi, na.rm = T)
  small_ui = xi - small_x_bar
  #then calculate the uncertainties:
  sigma_m_sq = 1/(sum(Wi * (small_ui)^2))
  sigma_b_sq = (1/sum(Wi)) + sigma_m_sq*(small_x_bar^2)
  #note = the cantrell and york inplementations are equivilent! yay. 

  # #Cantrells uncertainty assessments: 
  # sigma_m_sq <- 1/sum(Wi*(beta_val - beta_bar)^2) #Sigma slope
  # sigma_b_sq <- 1/sum(Wi) + (x_bar + beta_bar)^2 * sigma_m_sq
  # #edited 2020/03/02 to add sqrts - from R package IsoplotR, york function <3
  # #Additional note 2020/04/07 - got this wrong, need to be sqrted after calcs since sigma_m_sq is used in next calc. 

  #S from Cantrell 2008 (gives the same answer as Cantrell's sheet and York 2004)
  S <- sum(Wi*(data_function$y - (m_calc*data_function$x + intercept_calc))^2) 
  gofit <- sqrt(S/(n_obs - 2)) #Form from Cantrell 2008

  #######################################
  #Calculate the correlation coefficient:
  r <- cor(data_function$x, data_function$y, use = "pairwise.complete.obs")
  r_sq <- r^2

  #yay!
  coefficient_values <- as.data.frame(cbind(m_calc, intercept_calc, r_sq, sqrt(sigma_m_sq), sqrt(sigma_b_sq), S, gofit, i))
  names(coefficient_values) <- c("slope", "intercept", "r_sq", "uncertainty_slope", "uncertainty_intercept", "S", "gofit", "converged_cycles")

  return(coefficient_values)

}
