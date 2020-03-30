\name{cv.grbLMM}
\alias{cv.grbLMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
k-fold cross validation of grbLMM()
}
\description{
Performs k-fold cross validation of grbLMM() by splitting the data into subsets and executing grbLMM() with different training and test data k times.
}
\usage{
cv.grbLMM(k, y, X, Z, id, m.stop = 500, ny = .1, cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
Number of folds for the cross validation
}
\item{y}{
The dependent variable
}
\item{X}{
Numeric design matrix for the fixed effects (no intercept column)
}
\item{Z}{
Numeric design matrix for the random effects, same shape as X with intercepts in first column
}
\item{id}{
Grouping factor
}
\item{m.stop}{
Given number of maximum iterations
}
\item{ny}{
Learning rate, i.e. step length, of the algorithm
}
\item{cores}{
Number of cores for parallel computing. Requires the parallel package.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(666, kind = "L'Ecuyer-CMRG")

# basic example using Orthodont data
y_Orth = Orthodont$distance
X_Orth = cbind(Orthodont$Sex, Orthodont$age)
Z_Orth = matrix(1, length(Orthodont$distance), 1)
id_Orth = Orthodont$Subject

# perform 5-fold cross validation with a total of 200 iterations using 5 cores parallel
model.bst = cv.grbLMM(k = 5, y_Orth, X_Orth, Z_Orth, id_Orth, m.stop = 200, cores = 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
