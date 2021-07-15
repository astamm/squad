#################### FORMATTING CODE TO BE QTS-ORIENTED #####################
library(tidyverse)
library(squat)

choix_nb.axis<-function(){
  cat("Combien d'axes souhaitez-vous conserver ?")
  nb.axis<-readline()
  return(as.numeric(nb.axis))
}

#' "Classic" quaternion to Cayley-Dickson notation
#'
#' @param q a quaternion as a tibble 1x5 with columns time, w, x, y, z
#'
#' @return a list of two tibbles containing two columns time, cplx
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' qts_q2cdq(q1)
qts_q2cdq<-function(q){
  A1<- tibble(time=q %>% pull(time), cplx = complex(real = q %>% pull(w),imaginary = q %>% pull(x)))
  A2<- tibble(time=q %>% pull(time), cplx = complex(real = q %>% pull(y),imaginary = q %>% pull(z)))
  return(list(A1,A2))
}


#' Cayley-Dickson notation to "classic" quaternion
#'
#' @param cdq a list of two tibbles with 2 columns time, cplx
#'
#' @return a quaternion as a tibble 1x5 with columns time, w, x, y, z
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' cdq<-qts_q2cdq(q1)
#' qts_cdq2q(cdq)
qts_cdq2q<-function(cdq){
  if(!identical(cdq[[1]] %>% pull(time), cdq[[2]] %>% pull(time))) stop('q1 and q2 should be same time')
  q<-tibble(time=cdq[[1]] %>% pull(time),w=Re(cdq[[1]] %>% pull(cplx)), x=Im(cdq[[1]] %>% pull(cplx)),
            y=Re(cdq[[2]] %>% pull(cplx)),z=Im(cdq[[2]] %>% pull(cplx)))
  return(q)
}

#' Conjugate of a quaternion
#'
#' @param q a quaternion as a tibble 1x5 with columns time, w, x, y, z
#'
#' @return a quaternion as a tibble 1x5 with columns time, w, x, y, z
#' @export
#'
#' @examples
#' q<-tibble(time=0,w=1,x=2,y=3,z=4)
#' qts_conjugate_q(q)
qts_conjugate_q<-function(q){
  for(i in 3:5){
    q[i]<- -q[i]
  }
  return(q)
}

#' Sum of two quaternions
#'
#' @param q1 a quaternion as a tibble 1x5 with columns time, w, x, y, z
#' @param q2 a quaternion as a tibble 1x5 with columns time, w, x, y, z
#'
#' @return a quaternion as a tibble 1x5 with columns time, w, x, y, z
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q2<-tibble(time=0, w=2,x=2,y=0,z=0)
#' qts_sum_q(q1,q2)
qts_sum_q<-function(q1,q2){
  if(!identical(q1 %>% pull(time), q2 %>% pull(time))) stop('q1 and q2 should be same time')
  return(as_tibble(mapply('+', q1[-1],q2[-1], SIMPLIFY = FALSE)) %>% add_column(time=q1 %>% pull(time), .before=1))
}

#' Size of a matrix of quaternions formatted as a list of qts
#'
#' @param M a list of quaternions time series
#'
#' @return a vector with number of line, number of columns of the matrix
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1,w=1,x=2,y=0,z=3))
#' Mat<-list(q1,q3)
#' qts_size_MatQ(Mat)
qts_size_MatQ<-function(M){
  ok<-TRUE
  temp<-dim(M[[1]])[1]
  for(i in 1:length(M)){
    if(!identical(dim(M[[i]])[1], temp)) ok<-FALSE
    temp<-dim(M[[i]])[1]
  }
  if(ok==FALSE) stop('This is not a list of qts')
  return(c(length(M), dim(M[[1]])[1]))
}

#' Sum of two quaternion matrixes
#'
#' @param A1 a list of qts
#' @param A2 a list of qts with same dimensions as A1
#'
#' @return a list of qts representing A1+A2
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1,w=1,x=2,y=0,z=3))
#' Mat<-list(q1,q3)
#' Mat2<-list(q3,q1)
#' qts_sum_MatQ(Mat,Mat2)
qts_sum_MatQ<-function(A1,A2){
  if(!identical(qts_size_MatQ(A1), qts_size_MatQ(A2))) stop('A1 and A2 should be the same size')
  m<-mapply(qts_sum_q,A1,A2, SIMPLIFY = FALSE)
  return(as.list(m, .name_repair = "minimal"))
}

#' Conjugate of a quaternion matrix
#'
#' @param X a list of qts
#'
#' @return a list of qts representing the conjugate of X
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1,w=1,x=2,y=0,z=3))
#' Mat<-list(q1,q3)
#' qts_conjugate_MatQ(Mat)
qts_conjugate_MatQ<-function(X){
  return(as.list(mapply(qts_conjugate_q, X, SIMPLIFY = FALSE)))
}

#' Tranpose a quaternion matrix
#'
#' @param X a list of qts
#'
#' @return a list of qts representing the transpose of X
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1,w=1,x=2,y=0,z=3))
#' Mat<-list(q1,q3)
#' qts_transpose_MatQ(Mat)
qts_transpose_MatQ<-function(X){
  l<-list()
  for(i in 1:dim(X[[1]])[1]){
    t<-tibble(time=integer(),w=double(),x=double(),y=double(),z=double())
    for(j in 1:length(X)){
      t<-t %>% add_row(X[[j]] %>% slice(i))
    }
    t<-t %>% mutate(time=0:(length(X)-1))
    l[[i]]<-t
  }
  return(l)
}

#' Conjugate transpose of a quaternion matrix
#'
#' @param X a list of qts
#'
#' @return a list of qts representing the conjugate transpose of X
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1,w=1,x=2,y=0,z=3))
#' Mat<-list(q1,q3)
#' qts_transconj_MatQ(Mat)
qts_transconj_MatQ<-function(X){
  return(qts_transpose_MatQ(qts_conjugate_MatQ(X)))
}

#' Product of a qts with a scalar
#'
#' @param a a scalar
#' @param qts a quaternion time series
#'
#' @return a qts representing the input qts multiplied with a
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' qts_mult(3,q1)
qts_mult<-function(a, qts){
  qts<-qts %>% mutate(w=a*w,x=a*x,y=a*y,z=a*z)
  return(qts)
}

#' Product of a quaternion matrix with a scalar
#'
#' @param a a scalar
#' @param X a list of qts
#'
#' @return a list of qts representing the input list of qts multiplied by a
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1,w=1,x=2,y=0,z=3))
#' Mat<-list(q1,q3)
#' qts_mult_MatQ(3,Mat)
qts_mult_MatQ<-function(a,X){
  return(lapply(X=X,FUN=qts_mult, a=a))
}


#' Quaternion matrix in classic notation to Cayley-Dickson notation
#'
#' @param A a list of qts
#'
#' @return a list of N lists (length(A)=N) containing two tibbles with two columns time, cplx
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1,w=1,x=2,y=0,z=3))
#' Mat<-list(q1,q3)
#' qts_MatQtoMatcdQ(Mat)
qts_MatQtoMatcdQ<-function(A){
  m<-mapply(FUN = qts_q2cdq, A, SIMPLIFY = FALSE)
  return(m)
}

#' Quaternion matrix in Cayley-Dickson notation to 'classic' notation
#'
#' @param A a list of lists containing two tibbles with two columns time, cplx
#'
#' @return a list of qts
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1,w=1,x=2,y=0,z=3))
#' Mat<-list(q1,q3)
#' qts_MatcdQtoMatQ(qts_MatQtoMatcdQ(Mat))
qts_MatcdQtoMatQ<-function(A){
  l<-list()
  for(i in 1:length(A)){
    l[[i]]<-qts_cdq2q(A[[i]])
  }
  return(l)
}

#' Adjoint complex matrix
#'
#' @description this function returns the adjoint complex matrix associated to A = A1 + A2.j 
#' (see Le Bihan and Mars, "Singular value decomposition of quaternion matrices: a new tool for vector-sensor signal processing",
#' Signal Processing 84 (2004) 1177–1199 for more information)
#' 
#'
#' @param A a list of qts
#'
#' @return a complex matrix
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1,w=1,x=2,y=3,z=4))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1,w=1,x=2,y=0,z=3))
#' Mat<-list(q1,q3)
#' qts_MatHtoC(Mat)
qts_MatHtoC<-function(A){
  
  cdA<-qts_MatQtoMatcdQ(A)
  
  N<-length(cdA)
  M<-dim(cdA[[1]][[1]])[1]
  
  A1<-t(as.matrix(cdA[[1]][[1]] %>% select(cplx)))
  A2<-t(as.matrix(cdA[[1]][[2]] %>% select(cplx)))
  
  for(i in 2:length(cdA)){
    A1<-rbind(A1,t(as.matrix(cdA[[i]][[1]] %>% select(cplx))))
    A2<-rbind(A2,t(as.matrix(cdA[[i]][[2]] %>% select(cplx))))
  }
  ChiA<-cbind(A1,A2)
  temp<-cbind(-Conj(A2), Conj(A1))
  ChiA<-rbind(ChiA,temp)

  return(ChiA)
}

#' Quaternion matrix product
#'
#' @param A a list of qts
#' @param B a list of qts
#'
#' @return a quaternion matrix in classic notation (list of qts)
#' @export
#'
#' @examples
#' q1<-tibble(time=0,w=1,x=0,y=0,z=0)
#' q1<-q1 %>% add_row(tibble(time=1:2,w=c(1,1),x=c(2,3),y=c(3,5),z=c(4,9)))
#' q3<-tibble(time=0, w=2,x=1,y=3,z=0)
#' q3<-q3 %>% add_row(tibble(time=1:2,w=c(2,1),x=c(4,3),y=c(0,5),z=c(1,0)))
#' Mat<-list(q1,q3)
#' tMat<-qts_transpose_MatQ(Mat)
#' qts_Hmatprod(Mat,tMat)
qts_Hmatprod<-function(A,B){
  
  cdA<-qts_MatQtoMatcdQ(A)
  cdB<-qts_MatQtoMatcdQ(B)
  
  A1<-t(as.matrix(cdA[[1]][[1]] %>% select(cplx)))
  A2<-t(as.matrix(cdA[[1]][[2]] %>% select(cplx)))
  B1<-t(as.matrix(cdB[[1]][[1]] %>% select(cplx)))
  B2<-t(as.matrix(cdB[[1]][[2]] %>% select(cplx)))
  
  if(length(cdA)!=1){
    for(i in 2:length(cdA)){
      A1<-rbind(A1,t(as.matrix(cdA[[i]][[1]] %>% select(cplx))))
      A2<-rbind(A2,t(as.matrix(cdA[[i]][[2]] %>% select(cplx))))
    }
  }
  if(length(cdB)!=1){
    for(j in 2:length(cdB)){
      B1<-rbind(B1,t(as.matrix(cdB[[j]][[1]] %>% select(cplx))))
      B2<-rbind(B2,t(as.matrix(cdB[[j]][[2]] %>% select(cplx))))
    }
  }
  AB1<-A1%*%B1-A2%*%Conj(B2)
  AB2<-A1%*%B2+A2%*%Conj(B1)
  
  AB<-list()
  for(k in 1:dim(AB1)[1]){
    AB[[k]]<-tibble(time=0:(dim(AB1)[2]-1), w=Re(AB1[k,]), x=Im(AB1[k,]), y=Re(AB2[k,]), z=Im(AB2[k,]))
  }
  return(AB)
}


#' Singular Value Decomposition of a quaternion matrix
#'
#' @param A a quaternion matrix as a list of qts
#'
#' @description This function computes the complete or restreint SVD of a quaternion matrix A using 
#' the complex adjoint matrix method (see Le Bihan and Mars, "Singular value decomposition of quaternion matrices: a new tool for vector-sensor signal processing",
#' Signal Processing 84 (2004) 1177–1199 for more information)
#' The user chooses the number of axis to keep for the decomposition
#'
#' @return U : a list of qts representing the quaternion matrix which contains the left eigenvectors
#' V : a list of qts representing the quaternion matrix which contains the right eigenvectors
#' D : a vector containing the real singular values
#' inertia : the percentage of inertia represented by the nb first nb.axis (nb.axis is chosen by the user)
#' svdh_A : a list of qts which is the reconstruction of the initial matrix using the SVD on the first nb.axis axis
#' nb.axis : the number of axis chosen by the user
#' 
#' @export
#'
#' @examples
#' l<-list()
#' for(i in 1:10){
#'   l[[i]]<-tibble(time=0:9, w=runif(10), x=runif(10),y=runif(10),z=runif(10))
#' }
#' qts_SVDH(l)
qts_SVDH<-function(A){
  
  cdA<-qts_MatQtoMatcdQ(A)
  
  N<-length(A)
  M<-dim(A[[1]])[1]

  #Forming the complex representation
  ChiA<-qts_MatHtoC(A)
  
  #Computing the SVD of ChiA
  svd<-svd(ChiA,nu=nrow(ChiA),nv=ncol(ChiA))
  Uc<-svd$u
  Dc<-svd$d
  Vc<-svd$v
  #allocating singular values of ChiA to those of A
  sing_values<-Dc[seq(1,length(Dc),2)]
  eigenvalues<-sing_values**2
  
  #allocating singular vectors of ChiA to those of A
  UC<-Uc[,seq(1,dim(Uc)[2],2)]
  VC<-Vc[,seq(1,dim(Vc)[2],2)]
  U1<-UC[1:N,]
  U2<- -Conj(UC[(N+1):(2*N),])
  V1<-VC[1:M,]
  V2<- -Conj(VC[(M+1):(2*M),])
  
  U<-list()
  V<-list()
  for(i in 1:N){
    U[[i]]<-tibble(time=1:N, w=Re(U1[i,]), x=Im(U1[i,]), y=Re(U2[i,]), z=Im(U2[i,]))
  }
  for(j in 1:M){
    V[[j]]<-tibble(time=1:M, w=Re(V1[j,]), x=Im(V1[j,]), y=Re(V2[j,]), z=Im(V2[j,]))
  }
  Vt<-qts_transconj_MatQ(V)

  iner<-list()
  for(i in 1:min(N,M)){
    iner[i]<-eigenvalues[i]/sum(eigenvalues)*100
  }
  p <- ggplot(as.data.frame(as.numeric(iner)), aes(x=1:min(N,M),y=as.numeric(iner))) + geom_bar(stat = "identity", fill="steelblue") +
    geom_text(aes(label=round(as.numeric(iner),2)),vjust=-0.3, color="black", size=3) + labs(y="Inertia", x="Eigenvalues") +theme_minimal()
  print(p)
  nb.axis<-choix_nb.axis()
  
  inertia<-sum(eigenvalues[1:nb.axis])/sum(eigenvalues)*100
  
  #Reconstruction of A considering nb.axis
  prod<-qts_Hmatprod(qts_transpose_MatQ(qts_transpose_MatQ(U)[1]), Vt[1])
  s<-list()
  s[[1]]<-qts_mult_MatQ(sing_values[1],prod)
  svdh_A<-s[[1]]
  if(nb.axis!=1){
    for(n in 2:nb.axis){
      prod<-qts_Hmatprod(qts_transpose_MatQ(qts_transpose_MatQ(U)[n]), Vt[n])
      s[[n]]<-qts_mult_MatQ(sing_values[n],prod)
      svdh_A<-qts_sum_MatQ(svdh_A, s[[n]])
    }
  }
  return(list('U'=U,'V'=V, 'singular_values'=sing_values,'eigenvalues'=eigenvalues, 'inertia'=inertia, 'svdh_A'=svdh_A, 'nb.axis'=nb.axis))
}
