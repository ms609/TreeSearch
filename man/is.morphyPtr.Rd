% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpl_morphy_objects.R
\name{is.morphyPtr}
\alias{is.morphyPtr}
\title{Is an object a valid Morphy object?}
\usage{
is.morphyPtr(morphyObj)
}
\arguments{
\item{morphyObj}{Object of class \code{morphy}, perhaps created with
\code{\link[=PhyDat2Morphy]{PhyDat2Morphy()}}.}
}
\value{
\code{is.morphyPtr()} returns \code{TRUE} if \code{morphyObj} is a valid morphy
pointer, \code{FALSE} otherwise.
}
\description{
Is an object a valid Morphy object?
}
\seealso{
Other Morphy API functions: 
\code{\link{GapHandler}()},
\code{\link{MorphyErrorCheck}()},
\code{\link{MorphyWeights}()},
\code{\link{PhyDat2Morphy}()},
\code{\link{SingleCharMorphy}()},
\code{\link{UnloadMorphy}()},
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
\code{\link{mpl_set_parsim_t}()},
\code{\link{mpl_translate_error}()},
\code{\link{mpl_update_lower_root}()},
\code{\link{mpl_update_tip}()},
\code{\link{summary.morphyPtr}()}
}
\author{
\href{https://smithlabdurham.github.io/}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{Morphy API functions}
