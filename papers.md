# Papers Referenced in TreeSearch Optimization Work

Tracking citations for write-up purposes.

## Search Algorithm Infrastructure

| Key | Citation | Relevance |
|-----|----------|-----------|
| Goloboff1999 | Goloboff, P.A. (1999). Analyzing large data sets in reasonable times: solutions for composite optima. *Cladistics* 15(4): 415–428. doi:10.1006/clad.1999.0122 | Foundational: sectorial searches (XSS/RSS/CSS), drift search, outer-cycle interleaving pattern (§2.3) |
| Nixon1999 | Nixon, K.C. (1999). The Parsimony Ratchet, a new method for rapid parsimony analysis. *Cladistics* 15(4): 407–414. doi:10.1111/j.1096-0031.1999.tb00277.x | Parsimony ratchet: weight perturbation to escape local optima |
| Goloboff2016 | Goloboff, P.A. & Catalano, S.A. (2016). TNT version 1.5, including a full implementation of phylogenetic morphometrics. *Cladistics* 32(3): 221–238. doi:10.1111/cla.12160 | TNT reference: xmult, driven search pipeline |

## Wagner Tree Construction

| Key | Citation | Relevance |
|-----|----------|-----------|
| Goloboff2014 | Goloboff, P.A. (2014). Extended implied weighting. *Cladistics* 30(3): 260–272. doi:10.1111/cla.12047 | §3.3: biased taxon-addition order for Wagner trees; stochastic sampling from informativeness-weighted distribution (T-188) |
| Kluge1969 | Kluge, A.G. & Farris, J.S. (1969). Quantitative phyletics and the evolution of anurans. *Systematic Zoology* 18(1): 1–32. | Original Wagner tree construction |

## NNI Perturbation

| Key | Citation | Relevance |
|-----|----------|-----------|
| Nguyen2015 | Nguyen, L.-T., Schmidt, H.A., von Haeseler, A. & Minh, B.Q. (2015). IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. *Molecular Biology and Evolution* 32(1): 268–274. doi:10.1093/molbev/msu300 | Stochastic NNI-perturbation strategy (doRandomNNIs); topology-space escape mechanism (T-186) |

## Parsimony Scoring

| Key | Citation | Relevance |
|-----|----------|-----------|
| Fitch1971 | Fitch, W.M. (1971). Toward defining the course of evolution: Minimum change for a specific tree topology. *Systematic Zoology* 20(4): 406–416. | Standard Fitch parsimony algorithm |
| Brazeau2019 | Brazeau, M.D., Guillerme, T. & Smith, M.R. (2019). An algorithm for morphological phylogenetic analysis with inapplicable data. *Systematic Biology* 68(4): 619–631. doi:10.1093/sysbio/syy083 | Three-pass inapplicable algorithm (NA scoring) |

## Inapplicable-Handling Alternatives

| Key | Citation | Relevance |
|-----|----------|-----------|
| Hopkins2021 | Hopkins, M.J. & St. John, K. (2021). A new approach to inapplicable characters in phylogenetics. *Systematic Biology* 70(4): 764–781. | HSJ dissimilarity-metric scoring (T-116–T-118) |
| Goloboff2021 | Goloboff, P.A., Torres, A. & Arias, J.S. (2021). Weighted parsimony outperforms other methods of phylogenetic inference under models appropriate for morphology. *Cladistics* 37(6): 569–588. | X-transformation recoding for step-matrix inapplicable handling (T-122) |

## Implied Weighting

| Key | Citation | Relevance |
|-----|----------|-----------|
| Goloboff1993 | Goloboff, P.A. (1993). Estimating character weights during tree search. *Cladistics* 9(1): 83–91. | Implied weighting: k/(e+k) |
| Goloboff2018 | Goloboff, P.A., Torres, A. & Arias, J.S. (2018). Weighted parsimony outperforms other methods of phylogenetic inference under models appropriate for morphology. *Cladistics* 35(4): 407–437. doi:10.1111/cla.12prior | Recommended k≈10; IW outperforms EW on morphological data |

## Resampling

| Key | Citation | Relevance |
|-----|----------|-----------|
| Goloboff2003 | Goloboff, P.A. et al. (2003). Improvements to resampling measures of group support. *Cladistics* 19(4): 324–332. doi:10.1016/S0748-3007(03)00060-4 | Symmetric resampling (jackknife/bootstrap frequencies) |

## Profile Parsimony

| Key | Citation | Relevance |
|-----|----------|-----------|
| Faith2001 | Faith, D.P. & Trueman, J.W.H. (2001). Towards an inclusive philosophy for phylogenetic inference. *Systematic Biology* 50(3): 331–350. doi:10.1080/10635150118627 | Profile parsimony: information-based character weights |

