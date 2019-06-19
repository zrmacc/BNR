# Purpose: Score test for bivariate normal regression via least squares.
# Updated: 19/06/19

#' Bilateral Score Test via Least Squares.
#'
#' Performs a Score test of the null hypothesis that a subset of the regression
#' parameters are zero for both outcomes.
#'
#' @param y1 First outcome vector.
#' @param y2 Second outcome vector.
#' @param X Model matrix.
#' @param L Logical vector, with as many entires as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#'
#' @importFrom methods kronecker
#' @importFrom stats model.matrix pchisq resid vcov
#'
#' @return A numeric vector containing the score statistic, the degrees of
#'   freedom, and a p-value.

Score2 = function(y1,y2,X,L){
  # Test specification
  p = ncol(X);
  df = sum(L);
  if(length(L)!=p){stop("L should have one entry per column of X.")};
  if(df==0){stop("At least 1 entry of L should be TRUE.")};
  if(df==p){stop("At least 1 entry of L should be FALSE.")};

  ## Partition
  Xa = X[,L,drop=F];
  Xb = X[,!L,drop=F];

  # Fit null model
  M0 = fit.bnr(y1=y1,y2=y2,X=Xb);
  # Extract covariance
  S = vcov(M0,type="Outcome",inv=F);
  L = matInv(S);

  # Extract residuals
  e1 = resid(M0,type="Y1");
  e2 = resid(M0,type="Y2");

  ## Score
  U1 = array(0,dim=c(df,1));
  U1 = U1+L[1,1]*matIP(Xa,e1)+L[1,2]*matIP(Xa,e2);
  U2 = array(0,dim=c(df,1));
  U2 = U2+L[2,2]*matIP(Xa,e2)+L[2,1]*matIP(Xa,e1);
  U = rbind(U1,U2);

  ## Information
  # IPs
  Xa2 = matIP(Xa,Xa);
  Xb2 = matIP(Xb,Xb);
  Xab = matIP(Xa,Xb);

  # Target information
  Ibb = kronecker(L,Xa2);
  Iaa = kronecker(L,Xb2);
  Iba = kronecker(L,Xab);
  # Efficient information
  V = SchurC(Ibb=Ibb,Iaa=Iaa,Iba=Iba);

  ## Test
  # Statistic
  Ts = as.numeric(matQF(X=U,A=matInv(V)));
  # P value
  p = pchisq(q=Ts,df=2*df,lower.tail=F);
  # Output
  Out = c("Score"=Ts,"df"=2*df,"p"=p);
  return(Out);
}
