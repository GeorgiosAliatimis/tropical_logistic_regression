clust <- function(omega){
  e = length(omega)
  m = (1+sqrt(1+8*e))/2 # so that m choose 2 = e
  dist_mat = matrix(0,m,m)
  dist_mat[lower.tri(dist_mat)]= omega 
  dist_mat = dist_mat + t(dist_mat)
  u = hclust(as.dist(dist_mat),method="complete")
  u = as.phylo(u)
  # u = phangorn::upgma(dist_mat)
  C =cophenetic(u)
  C[lower.tri(C)]
}



