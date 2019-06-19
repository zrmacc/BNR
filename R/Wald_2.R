# Purpose: Bilateral Wald test for bivariate normal regression via least squares.
# Updated: 19/06/18

#' Bilateral Wald Test via Least Squares.
#'
#' Performs a Wald test of the null hypothesis that a subset of the regression
#' parameters are zero for both outcomes.
#'
#' @param y1 First outcome vector.
#' @param y2 Second outcome vector.
#' @param X Model matrix.
#' @param L Logical vector, with as many entires as columns in the design
#'   matrix, indicating which columns have coefficient zero under the null.
#'
#' @importFrom stats model.matrix pchisq resid vcov
#'
#' @return A numeric vector containing the score statistic, the degrees of
#'   freedom, and a p-value.

Wald2 = function(y1,y2,X,L){
  # Test specification
  p = ncol(X);
  df = sum(L);
  if(length(L)!=p){stop("L should have one entry per column of X.")};
  if(df==0){stop("At least 1 entry of L should be TRUE.")};
  if(df==p){stop("At least 1 entry of L should be FALSE.")};

  ## Model Fitting
  M0 = fit.bnr(y1=y1,y2=y2,X=X);

  # Extract information
  I = vcov(M0,type="Information",inv=F);

  # Information keys
  key0 = c(L,L);
  key1 = c(!L,!L);

  # Partition information
  Ibb = I[key0,key0,drop=F];
  Iaa = I[key1,key1,drop=F];
  Iba = I[key0,key1,drop=F];
  # Efficient information
  V = SchurC(Ibb=Ibb,Iaa=Iaa,Iba=Iba);

  ## Test
  # Coefficients of interest
  U = coef(M0)$Point[key0];
  # Statistic
  Tw = as.numeric(matQF(X=U,A=V));
  # P value
  p = pchisq(q=Tw,df=2*df,lower.tail=F);
  # Output
  Out = c("Wald"=Tw,"df"=2*df,"p"=p);
  return(Out);
}
