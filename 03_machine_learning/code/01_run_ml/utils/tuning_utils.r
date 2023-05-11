make_tuner <- function(inner,
                       measure,
                       learner,
                       ps,
                       term_evals,
                       n_evals_afs,
                       method_at,
                       fselector,
                       method_afs
                       ) {

    # Autotuner
    at <- auto_tuner(
            method = method_at,
            learner = learner,
            search_space = ps,
            resampling = inner,
            measure = measure,
            term_evals = term_evals
    )
    # Autotuner features
    if (fselector == TRUE) {
            at <- auto_fselector(
                        method = method_afs,
                        learner = at,
                        resampling = inner,
                        measure = measure,
                        term_evals = n_evals_afs
    )
    }

    return(at)
}