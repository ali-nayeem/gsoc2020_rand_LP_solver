library("volesti")
library("linprog")

RCP_LP <- function(A, b, c, maxiter=500, sample=100, maximum = FALSE, verbose = TRUE ) {
  
  iterVec <- c(0)
  scoreVec <- c(0)
  
  if (maximum) {
    c <- c*-1
  }

  H <- Hpolytope$new(A,b)
  ind <- rep(0, 500)
  perf <- rep(0, 500)
  for (i in 1:maxiter)
  {
    
    points2 <- sample_points(H, sample, random_walk = list("walk" = "RDHR", "walk_length" = 3)) #, 'CDHR' "RDHR"
    plot(points2[1,], points2[2,])
    points3 <- t(points2)%*%c
    min <- which.min(points3[,1])
    xmin <- points2[,min]
    # print(min)
    # print(xmin)
    # print(points3[min])
    iterVec <- append(iterVec, i)
    scoreVec <- append(scoreVec, points3[min]) 
    score <- points3[min]
    if (maximum) {
      score <- score*-1
    }
    if (verbose) {
      print(sprintf(fmt = "Iteration: %d, score: %g", i, score))
    }
    if (sd(points3) < 0.00001) {
      break
    }

    H$A <- rbind(H$A, c) #*-1
    H$b <- append(H$b, points3[min]) #*-1
    
  }
  if (maximum) {
    scoreVec <- scoreVec*-1
  }
  return(list("score"=score, "point"=xmin, "iterV"=iterVec, "scoreV"=scoreVec))
}


# maximize 2x1 + 4x2 + 3x3
# subject to :
#   3x1 + 4x2 + 2x3 =< 60
#   2x1 + x2 + 2x3 =< 40
#   x1 + 3x2 + 2x3 =< 80
#   x1, x2, x3 >= 0

cvec <- c(2, 4, 3) 
bvec <- c(60, 40, 80, 0, 0, 0) 
Amat <- rbind( c( 3, 4, 2 ),
               c( 2, 1, 2 ),
               c( 1, 3, 2 ),
               c( -1, 0, 0 ),
               c( 0, -1, 0 ),
               c( 0, 0, -1 ))
res <- RCP_LP(Amat, bvec, cvec, maximum = TRUE)
plot(res$iterV, res$scoreV, type="s", xlab = 'Iteration', ylab = 'Score')
print(str(res))

# bench = readMps('afiro.mps', solve=FALSE)
# res <- RCP_LP(bench$Amat, bench$bvec, bench$cvec)
# plot(res$iterV, res$scoreV, type="s", xlab = 'Iteration', ylab = 'Score')
# print(str(res))