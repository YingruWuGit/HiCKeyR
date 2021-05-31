\name{segment}
\alias{segment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{segment
%%  ~~function to do ... ~~
}
\description{This function and the corresponding class were written in C++ and compiled under C++14 standard. It uses change-point testing to decipher hierarchical organization of TADs in HiC data and outputs change-point locations, orders and p-values.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
segment(argv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{argv}{Full path and name of the text file "arguments_HiCKey.txt" which contains arguments setting.
%%     ~~Describe \code{argv} here~~
}
}
\details{Please refer to https://github.com/YingruWuGit/HiCKey for arguments setting in "arguments_HiCKey.txt".
%%  ~~ If necessary, more details than the description above ~~
}
\value{It does not return a value, but generates an output file containing change-point locations, orders and p-values.
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
segment("C:/Users/Yingru/Source/Repos/HiCKey/examples/arguments_HiCKey.txt")
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
