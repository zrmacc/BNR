# Purpose: Master testing function for bivariate normal regression
# Updated: 19/06/19

#' Test Bivariate Normal Regression Model.
#'
#' Inference procedure for bivariate normal regression models. Performs bilateral
#' tests that a subset of the regression parameters are zero for both outcomes,
#' and unilateral tests that a subset of the regression parameters are zero for
#' the first outcome only.
#'
#' @param y1 First outcome vector.
#' @param y2 Second outcome vector.
#' @param X Model matrix.
#' @param L Logical vector, with as many entires as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#' @param test Either Score or Wald. Only Wald is available with unilateral.
#' @param both If TRUE, performs the bilateral test.
#'
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#'
#' @return A numeric vector containing the test statistic, the degrees of
#'   freedom, and a p-value.

test.bnr = function(y1,y2,X,L,test="Wald",both=TRUE){
  # Input check
  if(!is.vector(y1)){stop("A numeric vector is expected for y1.")};
  if(!is.vector(y2)){stop("A numeric vector is expected for y2.")};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  if((sum(is.na(X))>0)){stop("Missing values are not expected in the covariate matrix.")}
  if(!is.logical(L)){stop("A logical vector is expected for L.")};
  if(!(test%in%c("Score","Wald"))){stop("Please selection either: Score or Wald.")};

  # Check for score test
  score = (test=="Score");

  # Cases
  if(!both){
    Out = Wald1(y1=y1,y2=y2,X=X,L=L);
  } else if(score){
    Out = Score2(y1=y1,y2=y2,X=X,L=L);
  } else {
    Out = Wald2(y1=y1,y2=y2,X=X,L=L);
  }

  # Output
  return(Out);
}
