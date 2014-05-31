countDna <- function(file) {
    data <- readLines(file)
    
    res <- split(seq(nchar(data)), unlist(strsplit(data, '')))
    a <- length(res[['A']])
    g <- length(res[['C']])
    c <- length(res[['G']])
    t <- length(res[['T']])
    print(c(a, g, c, t))
}