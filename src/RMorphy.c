#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <float.h>
#include <R.h>
#include <Rinternals.h>
#include "mpl.h"
#include "RMorphyUtils.h"
#include "RMorphy.h"

SEXP _R_wrap_mpl_new_Morphy(void) {
    Morphy new_Morphy = mpl_new_Morphy();
    SEXP result = PROTECT(R_MakeExternalPtr(new_Morphy, R_NilValue, R_NilValue));
    R_PreserveObject(result);
    R_RegisterCFinalizerEx(result, _finalize_Morphy, TRUE);
    UNPROTECT(1);
    
    return result;
}

void _finalize_Morphy (SEXP MorphyHandl) {
    Morphy handl = R_ExternalPtrAddr(MorphyHandl);
    if (handl == NULL) return;
    mpl_delete_Morphy(handl);
    R_ClearExternalPtr(MorphyHandl);
    R_ReleaseObject(MorphyHandl);
}

SEXP _R_wrap_mpl_delete_Morphy(SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));
    
    Morphy handl = R_ExternalPtrAddr(MorphyHandl);
    if (handl == NULL) {
        INTEGER(Rret)[0] = NA_INTEGER;
    } else {
        INTEGER(Rret)[0] = mpl_delete_Morphy(handl);
        R_ReleaseObject(MorphyHandl);
        R_ClearExternalPtr(MorphyHandl);
    }
    UNPROTECT(1);
    
    return Rret;
}

SEXP _R_wrap_mpl_init_Morphy(SEXP Rntax, SEXP Rnchar, SEXP MorphHandl) {
    int ret = 0;
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    Morphy handl = R_ExternalPtrAddr(MorphHandl);
    int ntax = INTEGER(Rntax)[0];
    int nchar = INTEGER(Rnchar)[0];

    ret = mpl_init_Morphy(ntax, nchar, handl);

    INTEGER(Rret)[0] = ret;
    UNPROTECT(1);
    return Rret;
}
    
SEXP _R_wrap_mpl_get_numtaxa(SEXP MorphHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = mpl_get_numtaxa(R_ExternalPtrAddr(MorphHandl));

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_get_num_charac(SEXP MorphHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));
    
    INTEGER(Rret)[0] = mpl_get_num_charac(R_ExternalPtrAddr(MorphHandl));

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_attach_symbols(SEXP Rsymbols, SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    int Mret = 0;
    const char *Msymbols = CHAR(asChar(Rsymbols));

    Mret = mpl_attach_symbols(Msymbols, R_ExternalPtrAddr(MorphyHandl));
    INTEGER(Rret)[0] = Mret;

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_get_symbols(SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(STRSXP, 1));

    char* symbols = mpl_get_symbols(R_ExternalPtrAddr(MorphyHandl));

    Rret = mkString(symbols);
    
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_attach_rawdata(SEXP Rmatrix, SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    int Mret = 0;
    const char *Mmatrix = CHAR(asChar(Rmatrix));

    Mret = mpl_attach_rawdata(Mmatrix, R_ExternalPtrAddr(MorphyHandl));
    
    INTEGER(Rret)[0] = Mret;
    UNPROTECT(1);

    return Rret;
}

SEXP _R_wrap_mpl_set_parsim_t(SEXP RcharID, SEXP Rchtype, SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));
    MPLchtype chtype;
    int Mret = 0;

    const char* chtypename = CHAR(asChar(Rchtype));
    
    chtype = _R_mpl_str2chtype(chtypename);
    Mret = mpl_set_parsim_t(INTEGER(RcharID)[0], chtype,
                            R_ExternalPtrAddr(MorphyHandl));

    INTEGER(Rret)[0] = Mret;

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_set_gaphandl(SEXP Rgaptype, SEXP MorphyHandl) {
	SEXP Rret = PROTECT(allocVector(INTSXP, 1));
	MPLgap_t gaptype;
	int Mret = 0;

	const char* gaptypename = CHAR(asChar(Rgaptype));
	
	gaptype = _R_mpl_str2gaptype(gaptypename);
	Mret = mpl_set_gaphandl(gaptype, R_ExternalPtrAddr(MorphyHandl));

	INTEGER(Rret)[0] = Mret;

	UNPROTECT(1);
	return Rret;
}

SEXP _R_wrap_mpl_get_gaphandl(SEXP MorphyHandl) {
	SEXP Rret = PROTECT(allocVector(INTSXP, 1));
	int Mret = 0;

	Mret = mpl_query_gaphandl(R_ExternalPtrAddr(MorphyHandl));

	INTEGER(Rret)[0] = Mret;

	UNPROTECT(1);
	return Rret;
}

SEXP _R_wrap_mpl_set_num_internal_nodes(SEXP Rnnodes, SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = mpl_set_num_internal_nodes
                        (INTEGER(Rnnodes)[0], 
                         R_ExternalPtrAddr(MorphyHandl));

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_get_num_internal_nodes(SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = mpl_get_num_internal_nodes
                        (R_ExternalPtrAddr(MorphyHandl));

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_apply_tipdata(SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = mpl_apply_tipdata(R_ExternalPtrAddr(MorphyHandl));
    
    UNPROTECT(1);    
    return Rret;
}


SEXP _R_wrap_mpl_set_charac_weight(SEXP RcharID, SEXP Rweight, 
                                   SEXP MorphyHandl) {
  SEXP Rret = PROTECT(allocVector(INTSXP, 1));
  INTEGER(Rret)[0] = mpl_set_charac_weight(INTEGER(RcharID)[0], REAL(Rweight)[0],
                                           R_ExternalPtrAddr(MorphyHandl));
  UNPROTECT(1);
  return Rret;
}

SEXP _R_wrap_mpl_get_charac_weight(SEXP RcharID, SEXP MorphyHandl) {
  SEXP Rret, Wapprox, Wexact;
  PROTECT(Rret = allocVector(VECSXP, 2));
  PROTECT(Wapprox = allocVector(INTSXP, 1));
  PROTECT(Wexact  = allocVector(REALSXP, 1));
  
  INTEGER(Wapprox)[0] = mpl_get_charac_weight(REAL(Wexact), INTEGER(RcharID)[0],
                                              R_ExternalPtrAddr(MorphyHandl));
  SET_VECTOR_ELT(Rret, 0, Wapprox);
  SET_VECTOR_ELT(Rret, 1, Wexact);
  UNPROTECT(3);
  return Rret;
}

SEXP _R_wrap_mpl_first_down_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id,
                                  SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_first_down_recon(INTEGER(Rnode_id)[0], INTEGER(Rleft_id)[0],
                         INTEGER(Rright_id)[0],
                         R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_first_up_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id,
                                SEXP Ranc_id, SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_first_up_recon(INTEGER(Rnode_id)[0], INTEGER(Rleft_id)[0],
                       INTEGER(Rright_id)[0], INTEGER(Ranc_id)[0],
                       R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_second_down_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id,
                                   SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_second_down_recon(INTEGER(Rnode_id)[0], INTEGER(Rleft_id)[0],
                       	  INTEGER(Rright_id)[0],
                          R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_second_up_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id,
                                 SEXP Ranc_id, SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_second_up_recon(INTEGER(Rnode_id)[0], INTEGER(Rleft_id)[0],
                       	INTEGER(Rright_id)[0], INTEGER(Ranc_id)[0],
                        R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_update_tip(SEXP tip_id, SEXP anc_id, SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_update_tip(INTEGER(tip_id)[0], INTEGER(anc_id)[0],
                        R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_update_lower_root(SEXP lower_id, SEXP upper_id,
                                   SEXP MorphyHandl) {
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_update_lower_root(INTEGER(lower_id)[0], INTEGER(upper_id)[0],
                        R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

void morphy_length (const int *ancestor, const int *left, const int *right, 
                    Morphy handl, int *score) {
  int i;
  const int n_taxa = mpl_get_numtaxa(handl); 
  const int n_internal = mpl_get_num_internal_nodes(handl);
  const int root_node = n_taxa;
  const int max_node = n_taxa + n_internal;
  
  for (i = max_node - 1; i >= n_taxa; i--) { /* First Downpass */
    *score += mpl_first_down_recon(i,  left[i - n_taxa], right[i - n_taxa],
                                   handl);
  }
 /* We could use a spare internal node with index = max_node as a dummy root node.
  * Instead: just pass the root node as its own ancestor. */
  mpl_update_lower_root(root_node, root_node, handl); 
 
  for (i = root_node; i != max_node; i++) { /* First uppass: internal nodes */
    *score += mpl_first_up_recon(i, left[i - n_taxa], right[i - n_taxa],
                                 ancestor[i], handl);
  }
  for (i = 0; i != n_taxa; i++) { /* First uppass: update tips */
    mpl_update_tip(i, ancestor[i], handl);
  }
  
  for (i = max_node - 1; i >= n_taxa; i--) { /* Second Downpass */
    *score += mpl_second_down_recon(i, left[i - n_taxa], right[i - n_taxa],
                                    handl);
  }
 
}

SEXP MORPHYLENGTH(SEXP R_ancestors, SEXP R_left, SEXP R_right, SEXP MorphyHandl) {
  Morphy handl = R_ExternalPtrAddr(MorphyHandl);
  /* R_descendants and R_ancestors have already had one subtracted to convert 
   * them to an index */
  const int
    /* INTEGER gives pointer to first element of an R vector*/
    *ancestor = INTEGER(R_ancestors),
    *left = INTEGER(R_left),
    *right = INTEGER(R_right)
  ;
            
  /* Declare and protect result, to return to R */
  SEXP Rres = PROTECT(allocVector(INTSXP, 1));
  
  /* Initialize return variables */
  int *score;
  score = INTEGER(Rres);
  *score = 0;
  morphy_length(ancestor, left, right, handl, score); /* Updates score */
   
  UNPROTECT(1);
  return Rres;
}
