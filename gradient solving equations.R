refr <- function(x, borders = c(2, 3), refr_level = 1 + 0i) {
  ifelse(x < borders[2] & x > borders[1], refr_level, 1 + 0i)
}




refr_matrix <- matrix(refr(x, borders = c(2, 3), refr_level = 5 + 2i), 
                      nrow = length(x), ncol = length(x), byrow = TRUE)



library(geometry)

norm2 <- function(a){
  a <- Mod(a)
  return(dot(a,a))
  
}
conjA <- function(A){
  return(Conj(t(A)))
}

grad_solve <- function(A, f, n_iter = 10e4, eps = 10e-7){
  n <- length(A[1,])
  
  if(length(A[1,])!=length(A[,1])){
    warning("Матрица должна быть квадратной")
    return(NA)
  }
  
  if(n!=length(f)){
    warning("Не совпадают размерности векторов")
    return(NA)
  }
  
  uk <- rep(0, n)
  for(i in 1:n_iter){
    
    uk0 <- uk
    rk <- A%*%uk0 - f
    
    crk <- conjA(A)%*%rk
    uk <- uk0 - crk*(norm2(crk)/norm2(A%*%crk))
    diff <- Mod(uk - uk0)
    if(max(diff)<eps){
      break
    }
  }
  print(paste("Итераций всего", i))
  return(uk)
}



a = 1
b = 4
N = 1000
h = (b - a)/ N
x = seq((a + h/2), b, h)

k = 4
G <- function(x,  k = 1) {
  return(exp(-1i * k * abs(x)) / (2 * 1i * k))
}

A <- matrix(x, nrow = length(x), ncol = length(x), byrow = FALSE)
for(i in 1:length(x)){
  A[i,] <- A[i,]-x
}

A <- refr_matrix * G(A,k)*h*(-(k**2))
"A <- G(A,k)*h*(-(k**2))"

diag(A) <- diag(A) + 1

f = cos(k*x)



u <- grad_solve(A,f, n_iter = 1000)


head(u, n = 10)

resid <- A%*%u-f
head(resid, n = 10)

plot(x, Re(u))
abline(v = c(2, 3), col = "red")


v <- solve(A,f)
A%*%v-f
head(v)
head(u)
