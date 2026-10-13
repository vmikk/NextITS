## The OptimOTU rank loop, with a stop rank for taxonomic splitting
##
## Sourced by optimotu_cluster.R after optimotu is loaded. This is optimotu::optimotu() from optimotu
## 0.10.0.9001 (R/optimotu.R; commit b10b06b, MIT, (c) 2025 Brendan Furneaux)
## with three changes, and no other differences in the clustering logic:
##
## 1. `stop_rank`. Every step at rank k works inside one group of rank k-1, and a
##    sequence clustered de novo gets all of its lower ranks in that same
##    seq_cluster() call. So once the loop has passed rank R, the sequences that
##    are classified or attached at R only ever interact with sequences sharing
##    their rank-R name. With `stop_rank = R` they are returned as `pending`
##    (seq_id plus their resolved ranks down to R) instead of being clustered
##    further; the loop still runs to the last rank, but only to finish naming
##    the de novo pseudotaxa, which are complete and returned as `final`.
##    Clustering each rank-R taxon of `pending` separately, with its resolved
##    ranks given as certain calls, then reproduces an unsplit run.
##
##    Truncating `ranks` instead would not: calc_subtaxon_thresholds() makes the
##    de novo thresholds monotonic with cummax() over all lower ranks.
##
## 2. `closed_ref_from`. Closed-reference attachment runs only at this rank and
##    below; above it, sequences not classified at a rank go straight to de novo
##    clustering among themselves. Attachment compares every unknown with every
##    known sequence of the parent taxon, so at high ranks (an unknown phylum
##    within all Fungi) it is the most expensive step by far. NULL keeps
##    optimotu()'s behaviour (attachment at every rank).
##
## 3. Deterministic closed-reference ties. closed_ref_cluster() picks a random
##    reference when a query is equally close to several. The draw depends on the
##    state of the RNG and on the row order seq_search() returns, so it differs
##    between a split and an unsplit run. Here the tie goes to the smallest
##    reference ID.

## closed_ref_cluster(), with the tie rule above
closed_ref_cluster_det <- function(query, ref, threshold, ...) {
  last_out <- data.frame(seq_id = character(0), ref_id = character(0))
  out <- list(last_out)
  while (optimotu:::sequence_size(ref) > 0 && optimotu:::sequence_size(query) > 0) {
    result <- optimotu::seq_search(query = query, ref = ref, threshold = threshold, ...)
    if (nrow(result) > 0) {
      result <- result[result$dist <= threshold, c("seq_id", "ref_id", "dist"), drop = FALSE]
      if (nrow(result) > 0) {
        result <- result[order(result$seq_id, result$dist, result$ref_id), , drop = FALSE]
        result <- result[!duplicated(result$seq_id), c("seq_id", "ref_id"), drop = FALSE]
      }
    }
    if (nrow(last_out) == 0) {
      last_out <- result <- result[, c("seq_id", "ref_id")]
    } else {
      result <- result[, c("seq_id", "ref_id")]
      names(result)[2] <- "temp"
      result <- merge(result, last_out, by.x = "temp", by.y = "seq_id")
      last_out <- result <- result[, c("seq_id", "ref_id")]
    }
    out <- c(out, list(result))
    ref <- optimotu:::select_sequence(query, result$seq_id)
    query <- optimotu:::select_sequence(query, result$seq_id, negate = TRUE)
  }
  do.call(rbind, out)
}

## optimotu::optimotu() for a data.frame of sequences (seq_id, seq), whose inputs
## the caller has already validated
##
## @return list(final = data.frame(seq_id, <ranks>), pending = NULL or
##   data.frame(seq_id, <ranks up to stop_rank>))
cluster_by_rank <- function(seqs, tax_prob, threshold_optima, ranks, measure = NULL,
                            prob_thresh = NULL, dist_config, clust_config,
                            parallel_config, verbose = FALSE, stop_rank = NULL,
                            closed_ref_from = NULL) {
  if (!is.null(stop_rank) && !stop_rank %in% ranks[-c(1L, length(ranks))]) {
    fail("The split rank must be one of: ",
         paste(ranks[-c(1L, length(ranks))], collapse = ", "))
  }
  if (!is.null(closed_ref_from) && !closed_ref_from %in% ranks[-1L]) {
    fail("closed_ref_from must be one of: ", paste(ranks[-1L], collapse = ", "))
  }
  first_closed_ref <- if (is.null(closed_ref_from)) 2L else match(closed_ref_from, ranks)
  verbose <- as.integer(verbose)
  n_ranks <- length(ranks)
  pending <- NULL

  taxon_table <- tax_prob[tax_prob$rank == ranks[1], c("seq_id", "taxon")]
  names(taxon_table)[2] <- ranks[1]

  pseudotaxon_table <- data.frame(seq_id = character(0))
  pseudotaxon_table[[ranks[1]]] <- character(0)

  for (r in seq_len(n_ranks - 1)) {
    .rank <- ranks[r + 1]
    .parent_rank <- ranks[r]

    .parent_taxa <- taxon_table
    .parent_pseudotaxa <- pseudotaxon_table

    if (is.null(prob_thresh)) {
      known_taxa <- tax_prob$rank == .rank & !is.na(tax_prob$taxon)
    } else {
      known_taxa <- tax_prob$rank == .rank &
        !is.na(tax_prob$taxon) &
        tax_prob$prob >= prob_thresh
    }
    known_taxon_table <- tax_prob[known_taxa, c("seq_id", "taxon")]
    names(known_taxon_table)[2] <- .rank
    known_taxon_table <- merge(.parent_taxa, known_taxon_table, by = "seq_id", all.x = TRUE)

    preclosed_taxon_table <- split(known_taxon_table, known_taxon_table[[.parent_rank]])
    to_closedref_cluster <- vapply(
      preclosed_taxon_table,
      function(x) any(is.na(x[[.rank]])) & !all(is.na(x[[.rank]])),
      logical(1)
    )
    preclosed_taxon_table <- preclosed_taxon_table[to_closedref_cluster]
    # above closed_ref_from, unknowns are only clustered de novo
    if (r + 1L < first_closed_ref) preclosed_taxon_table <- preclosed_taxon_table[0]

    thresholds <- optimotu::calc_taxon_thresholds(
      rank = .parent_rank,
      taxon_table = known_taxon_table,
      optima = threshold_optima,
      ranks = ranks,
      measure = measure
    )

    clusters_closed_ref <- lapply(
      preclosed_taxon_table,
      function(parent_group) {
        if (verbose) {
          cat("Closed reference clustering within", .parent_rank,
              parent_group[[.parent_rank]][1], "\n")
        }
        unknowns <- is.na(parent_group[[.rank]])
        taxon <- parent_group[[.parent_rank]][1]
        clusters <- closed_ref_cluster_det(
          query = optimotu:::select_sequence(seqs, parent_group$seq_id[unknowns]),
          ref = optimotu:::select_sequence(seqs, parent_group$seq_id[!unknowns]),
          threshold = optimotu::threshold_as_dist(thresholds[taxon]),
          dist_config = dist_config,
          parallel_config = parallel_config,
          verbose = max(verbose - 1L, 0L)
        )
        new_taxa <- merge(clusters, parent_group[, c("seq_id", .rank)],
                          by.x = "ref_id", by.y = "seq_id")
        new_taxa[, c("seq_id", .rank)]
      }
    )
    clusters_closed_ref <- do.call(rbind, clusters_closed_ref)

    known_taxon_table[
      match(clusters_closed_ref$seq_id, known_taxon_table$seq_id),
      .rank
    ] <- clusters_closed_ref[[.rank]]

    predenovo_taxon_table <- known_taxon_table[is.na(known_taxon_table[[.rank]]), ]

    denovo_thresholds <- optimotu::calc_subtaxon_thresholds(
      rank = .parent_rank,
      taxon_table = predenovo_taxon_table,
      optima = threshold_optima,
      ranks = ranks,
      measure = measure
    )

    predenovo_taxon_table <- split(predenovo_taxon_table, predenovo_taxon_table[[.parent_rank]])

    clusters_denovo <- lapply(
      predenovo_taxon_table,
      function(parent_group) {
        .parent_taxon <- parent_group[[.parent_rank]][1]
        if (verbose) {
          cat("De novo clustering within", .parent_rank, .parent_taxon, "\n")
        }
        if (nrow(parent_group) > 1) {
          clusters <- optimotu::seq_cluster(
            seq = optimotu:::select_sequence(seqs, parent_group$seq_id),
            threshold_config = optimotu::threshold_set(
              if (.parent_taxon %in% names(denovo_thresholds)) {
                denovo_thresholds[[.parent_taxon]]
              } else {
                denovo_thresholds[["_NA_"]]
              }
            ),
            dist_config = dist_config,
            clust_config = clust_config,
            parallel_config = parallel_config,
            output_type = "matrix",
            verbose = max(0L, verbose - 1L)
          )
          clusters <- t(clusters)
          clusters <- as.data.frame(clusters, stringsAsFactors = FALSE)
          for (new_rank in names(clusters)) {
            parent_group[[new_rank]] <- clusters[[new_rank]]
          }
        } else {
          for (new_rank in c(.rank, optimotu:::subranks(.rank, ranks))) {
            parent_group[[new_rank]] <- 0L
          }
        }
        parent_group
      }
    )
    clusters_denovo <- do.call(rbind, clusters_denovo)

    if (is.null(clusters_denovo) || nrow(clusters_denovo) == 0) {
      clusters_denovo <- data.frame(seq_id = character(0))
      for (new_rank in optimotu:::superranks(.rank, ranks)) {
        clusters_denovo[[new_rank]] <- character(0)
      }
      for (new_rank in c(.rank, optimotu:::subranks(.rank, ranks))) {
        clusters_denovo[[new_rank]] <- integer(0)
      }
    }

    taxon_table <- known_taxon_table[!is.na(known_taxon_table[[.rank]]), , drop = FALSE]

    pseudotaxon_table <- rbind(clusters_denovo, .parent_pseudotaxa)
    pseudotaxon_table <- sort_by(pseudotaxon_table, ~seq_id)
    pseudotaxon_table[[.rank]] <- paste(pseudotaxon_table[[.parent_rank]],
                                        pseudotaxon_table[[.rank]])
    pseudotaxon_table[[.rank]] <- ordered(
      pseudotaxon_table[[.rank]],
      levels = unique(pseudotaxon_table[[.rank]]),
      labels = optimotu:::make_seq_names(
        n = length(unique(pseudotaxon_table[[.rank]])),
        prefix = paste0("pseudo", .rank, "_")
      )
    ) |>
      as.character()

    # hand the sequences resolved at the stop rank over to per-taxon runs; the
    # remaining iterations only name the pseudotaxa below it
    if (!is.null(stop_rank) && identical(.rank, stop_rank)) {
      pending <- taxon_table
      rownames(pending) <- NULL
      taxon_table <- taxon_table[0, , drop = FALSE]
    }
  }
  out <- sort_by(rbind(taxon_table, pseudotaxon_table), ~seq_id)
  rownames(out) <- NULL
  list(final = out, pending = pending)
}

## seq_cluster()'s data.frame method re-dispatches by evaluating a call to an
## unexported method in the caller's frame, so the loop has to run with the
## optimotu namespace as its enclosure, as optimotu() itself does
.loop_env <- new.env(parent = asNamespace("optimotu"))
.loop_env$fail <- fail
environment(closed_ref_cluster_det) <- .loop_env
.loop_env$closed_ref_cluster_det <- closed_ref_cluster_det
environment(cluster_by_rank) <- .loop_env
