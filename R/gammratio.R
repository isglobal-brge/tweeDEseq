
gammratio <- function(n,d) {
    aux <- gamma(d)
    if(n==0){
        aux
    } else {
      #..#for(i in 1:n)
        for(i in seq_len(n)) {
            aux <- aux*(n+d-i)/i
        }
    }
    aux
}
