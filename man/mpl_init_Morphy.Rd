% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpl_morphyex.R
\name{mpl_init_Morphy}
\alias{mpl_init_Morphy}
\title{Sets up the dimensions of the dataset.}
\usage{
mpl_init_Morphy(numtaxa, numchars, morphyobj)
}
\arguments{
\item{numtaxa}{The number of taxa (or tips/terminals).}

\item{numchars}{The number of characters (i.e. transformation series) in the
data set.}

\item{morphyobj}{An instance of the Morphy object.}
}
\value{
Morphy error code.
}
\description{
Provides initial dimensions for the dataset, which will
constrain any input matrix supplied to Morphy.
}
\seealso{
Other Morphy API functions: 
\code{\link{GapHandler}()},
\code{\link{MorphyErrorCheck}()},
\code{\link{MorphyWeights}()},
\code{\link{PhyDat2Morphy}()},
\code{\link{SingleCharMorphy}()},
\code{\link{UnloadMorphy}()},
\code{\link{is.morphyPtr}()},
\code{\link{mpl_apply_tipdata}()},
\code{\link{mpl_attach_rawdata}()},
\code{\link{mpl_attach_symbols}()},
\code{\link{mpl_delete_Morphy}()},
\code{\link{mpl_first_down_recon}()},
\code{\link{mpl_first_up_recon}()},
\code{\link{mpl_get_charac_weight}()},
\code{\link{mpl_get_gaphandl}()},
\code{\link{mpl_get_num_charac}()},
\code{\link{mpl_get_num_internal_nodes}()},
\code{\link{mpl_get_numtaxa}()},
\code{\link{mpl_get_symbols}()},
\code{\link{mpl_new_Morphy}()},
\code{\link{mpl_second_down_recon}()},
\code{\link{mpl_second_up_recon}()},
\code{\link{mpl_set_charac_weight}()},
\code{\link{mpl_set_num_internal_nodes}()},
\code{\link{mpl_set_parsim_t}()},
\code{\link{mpl_translate_error}()},
\code{\link{mpl_update_lower_root}()},
\code{\link{mpl_update_tip}()},
\code{\link{summary.morphyPtr}()}
}
\author{
Martin Brazeau
}
\concept{Morphy API functions}
\keyword{internal}
