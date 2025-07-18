% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpl_morphyex.R
\name{mpl_set_parsim_t}
\alias{mpl_set_parsim_t}
\title{Sets a character's parsimony function type}
\usage{
mpl_set_parsim_t(char_id, tname = "typename", morphyobj)
}
\arguments{
\item{char_id}{The number of the character (transformation series) as defined
in the input matrix.  The first character is numbered 1 (one).}

\item{tname}{The parsimony function type as defined in morphydefs.h}

\item{morphyobj}{An instance of the Morphy object.}
}
\value{
A Morphy error code.
}
\description{
Set the parsimony function type to one defined in the
morphydefs.h header file. Setting the character to type NONE_T will also
cause it to be excluded from any further calculations.
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
\code{\link{mpl_init_Morphy}()},
\code{\link{mpl_new_Morphy}()},
\code{\link{mpl_second_down_recon}()},
\code{\link{mpl_second_up_recon}()},
\code{\link{mpl_set_charac_weight}()},
\code{\link{mpl_set_num_internal_nodes}()},
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
