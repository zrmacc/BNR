# Purpose: Fitting function for bivariate normal regression.
# Updated: 19/06/19

#' Fit Bivariate Normal Regression Model.
#'
#' Estimation procedure for bivariate normal regression models.
#'
#' @param y1 First outcome vector.
#' @param y2 Second outcome vector.
#' @param X Target model matrix.
#' @param sig Significance level.
#'
#' @importFrom methods new
#' @importFrom stats coef pnorm qnorm resid
#' @export
#' @return An object of class 'bnr' with slots containing the estimated regression
#'  coefficients, the target-surrogate covariance matrix, the information matrices
#'  for regression parameters, and the residuals.

fit.bnr = function(y1,y2,X,sig=0.05){
  # Input check
  if(!is.vector(y1)){stop("A numeric vector is expected for y1.")};
  if(!is.vector(y2)){stop("A numeric vector is expected for y2.")};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};

  # Covariate dimension
  p = ncol(X);

  # Stage 1 regression
  m1 = fitOLS(y=y1,X=X);
  # Parameters
  beta1 = m1$Beta;
  Sigma11 = m1$V;
  # Save XtX
  XtX = m1$Ibb*Sigma11;
  rm(m1);

  # Stage 2 regression
  Z = cbind(y1,X);
  m2 = fitOLS(y=y2,X=Z);
  zeta = m2$Beta;
  Lambda22i = m2$V;
  rm(m2);

  # Recover original parameters
  delta = zeta[1];
  gamma = zeta[2:(p+1)];
  beta2  = delta*beta1+gamma;

  Sigma12 = delta*Sigma11;
  Sigma22 = Lambda22i+Sigma11*(delta)^2;

  Lambda11i = Sigma11-(Sigma12^2)/Sigma22;
  Lambda12 = -(delta/Lambda22i);

  # Final covariance
  Sigma = matrix(c(Sigma11,Sigma12,Sigma12,Sigma22),nrow=2);
  colnames(Sigma) = rownames(Sigma) = c("Y1","Y2");

  ## Information
  J = list();

  # For alpha
  J$I11 = XtX/Lambda11i;
  J$I22 = XtX/Lambda22i;
  J$I12 = XtX*Lambda12;

  # Overall
  Info = rbind(cbind(J$I11,J$I12),cbind(t(J$I12),J$I22));
  Infoi = matInv(Info);

  # Regression coefficients
  Point = c(beta1,beta2);
  # SE
  SE = sqrt(diag(Infoi));
  # CIs
  z = qnorm(p=1-(sig/2));
  L = Point-z*SE;
  U = Point+z*SE;
  P = 2*pnorm(q=abs(Point/SE),lower.tail=F);
  # Labeling
  Outcome = c(rep("Y1",p),rep("Y2",p));
  if(is.null(colnames(X))){colnames(X) = paste0("x",seq(1:ncol(X)))};
  Regressor = rep(colnames(X),times=2);
  Coeff = data.frame("Outcome"=Outcome,"Regressor"=Regressor,"Point"=Point,"SE"=SE,"L"=L,"U"=U,"p"=P);
  colnames(Info) = rownames(Info) = Coeff$Regressor;

  ## Residuals
  E = cbind(y1-MMP(X,beta1),y2-MMP(X,beta2));
  colnames(E) = c("Y1","Y2");

  ## Output
  Out = new(Class="bnr", Coefficients=Coeff,Covariance=Sigma,Information=Info,Residuals=E);
  return(Out);
};
