args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}else if (length(args)==1) {
  args[2] = "estimated_lambda"
}
x <- as.matrix(read.csv(args[1], header = F))
x <- sort(x)
l <- length(x)
alpha <- 0.025
x <- x[round(alpha*l):round((1-alpha)*l)]
dmode <- function(x) {
      den <- density(x, kernel=c("gaussian"))
        ( den$x[den$y==max(den$y)] )
    }
write(as.numeric(dmode(x)), file = args[2], append = FALSE, sep = " ")