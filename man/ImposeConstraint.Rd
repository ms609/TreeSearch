% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImposeConstraint.R
\name{ImposeConstraint}
\alias{ImposeConstraint}
\title{Force a tree to match a constraint}
\usage{
ImposeConstraint(tree, constraint)
}
\arguments{
\item{tree}{A tree of class \code{\link{phylo}}.}

\item{constraint}{An object of class \code{phyDat}; returned trees will be
perfectly compatible with each character in \code{constraint}.
See \href{https://ms609.github.io/TreeSearch/articles/inapplicable.html}{vignette}
for further examples.}
}
\value{
\code{ImposeConstraint()} returns a tree of class \code{phylo}, consistent
with \code{constraint}.
}
\description{
Modify a tree such that it matches a specified constraint.
This is at present a somewhat crude implementation that attempts to retain
much of the strucure of \code{tree} whilst guaranteeing compatibility with
each entry in \code{constraint}.
}
\examples{
library('TreeTools')
tips <- letters[1:9]
tree <- as.phylo(1, 9, tips)
plot(tree)

constraint <- StringToPhyDat('0000?1111 000111111 0000??110', tips, FALSE)
plot(ImposeConstraint(tree, constraint))
}
\author{
\href{https://smithlabdurham.github.io/}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
