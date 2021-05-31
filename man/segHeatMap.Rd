\name{segHeatMap}
\alias{segHeatMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{segHeatMap
%%  ~~function to do ... ~~
}
\description{This function is almost the same as the function segment. The only difference is that it returns a submatrix of the HiC data in which the upper triangular part is the original HiC data but the lower triangular part is constructed from the detected change-points and their hierarchical orders. User can draw a heatmap by the return matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
segHeatMap(argv, s, e)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{argv}{Full path and name of the text file "arguments_HiCKey.txt" which contains arguments setting
%%     ~~Describe \code{argv} here~~
}
  \item{s}{Start index of the submatrix. Its default value is 0. Note that the index is 0 based.
%%     ~~Describe \code{s} here~~
}
  \item{e}{End index of the submatrix. Its default value is -1, which means the end of the whole HiC matrix. Note that the last index of the submatrix is actually end - 1 following C++ convention.
%%     ~~Describe \code{e} here~~
}
}
\details{Please refer to https://github.com/YingruWuGit/HiCKey for arguments setting in "arguments_HiCKey.txt".
%%  ~~ If necessary, more details than the description above ~~
}
\value{It returns a Rcpp::NumericMatrix from s to e. Its upper triangular part is the original HiC data but the lower triangular part is constructed from the detected change-points and their hierarchical orders. It also generates an output file containing change-point locations, orders and p-values.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Yingru Wu
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
X = segHeatMap("C:/Users/Yingru/Source/Repos/HiCKey/examples/arguments_HiCKey.txt")
heatmap(X, scale = "none", Rowv = NA, Colv = NA, col = gray.colors(50, start = 1, end = 0, gamma = 0.15))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
