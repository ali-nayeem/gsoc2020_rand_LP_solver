library(volesti)
library(ggplot2)
pdf("plots.pdf")

for (seed in c(0, 1, 2, 3))
{

  P <- gen_rand_vpoly(2, 6, seed=seed)
  points1 <- sample_points(P, 500)
  g  <- plot(ggplot(data.frame(x=points1[1,], y=points1[2,])) + geom_point(aes(x=x, y=y))  + ggtitle(sprintf("Sampling from a random V-polytope: seed=%s", seed)))
  
  H <- gen_rand_hpoly(2, 6)
  points2 <- sample_points(H, 500)
  g  <- plot(ggplot(data.frame(x=points2[1,], y=points2[2,])) + geom_point(aes(x=x, y=y))  + ggtitle(sprintf("Sampling from a random H-polytope: seed=%s", seed)))
  
}
dev.off()