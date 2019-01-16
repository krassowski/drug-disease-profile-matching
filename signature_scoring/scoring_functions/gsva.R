gsva.with_probabilities <- function(
  expression, expression_classes, gene_sets, procedure,
  method='gsva', mx.diff=TRUE, permutations=1000, adjust.method='fdr',
  fdr.threshold=0.01, include_control=T, cores=1, progress=T,
  limit_to_gene_sets=F,
  ...
) {
  stopifnot(procedure %in% c('gene_permutation', 'bayes'))
  cat('GSVA with probabilities starting')
  design <- cbind(condition = expression_classes != 'normal')
  phenoData <- Biobase::AnnotatedDataFrame(data = as.data.frame(as.table(design)))
  row.names(phenoData) <- colnames(expression)
  expression_set <- Biobase::ExpressionSet(assayData = data.matrix(expression), phenoData = phenoData)

  expressions <- Biobase::exprs(expression_set)
  genesInExpressionData <- rownames(expressions)
  # also, see GSVA::computeGeneSetsOverlap()
  if (limit_to_gene_sets != FALSE) {
      gene_sets <- gene_sets[sapply(names(gene_sets), function(name) name %in% limit_to_gene_sets)]
  }

  subset <- gene_sets[
    sapply(
        gene_sets,
        function(genes) length(genes[genes %in% genesInExpressionData]) > 1
    )
  ]

  cat(paste('Method calculation', cores, 'cores'))

  result <- GSVA::gsva(expressions, subset, method = method, verbose = F, parallel.sz = cores, mx.diff = mx.diff, ...)

  if (procedure == 'gene_permutation')
  {
    cat('Permuting')

    genes_order <- rownames(expressions)
    result <- as.data.frame(result)

    if (include_control) {
        result$difference <- result$condition - result$control
    }
    else {
        result$difference <- result$condition
    }

    if (progress) {
        apply = pbapply::pblapply
        if (cores > 1) {
            apply = pbmcapply::pbmclapply
        }
    } else {
        apply = lapply
        if (cores > 1) {
            apply = parallel::mclapply
        }
    }
    options(mc.cores = cores)

    i_permutations <- c(1:permutations)
    permutation_effect_sizes <- apply(i_permutations, function(i) {
      rownames(expressions) <- sample(genes_order)  # a random permutation
      random_result <- GSVA::gsva(expressions, subset, method = method, verbose = F, parallel.sz = cores, mx.diff = mx.diff, ...)
      if (include_control) {
        random_effect_size <- random_result[,'condition'] - random_result[,'control']
      }
      else {
        random_effect_size <- random_result[,'condition']
      }
      random_effect_size
    })

    random_effect_sizes <- as.data.frame(permutation_effect_sizes, col.names = i_permutations)

    result$p_value <- mapply(
      function(effect_size, gene_set_name) {
        if (effect_size > 0) {
          is_random_permutation_more_extreme_cnt <- sum(random_effect_sizes[gene_set_name,] > effect_size)
        }
        else if (effect_size == 0) {
             is_random_permutation_more_extreme_cnt <- permutations
        }
        else {
          is_random_permutation_more_extreme_cnt <- sum(random_effect_sizes[gene_set_name,] < effect_size)
        }
        is_random_permutation_more_extreme_cnt / permutations
      },
      result$difference,
      rownames(result)
    )

    result$fdr <- p.adjust(result$p_value, method = adjust.method)

    colnames(result)[colnames(result) == 'fdr'] <- 'fdr_q-val'
    colnames(result)[colnames(result) == 'difference'] <- 'nes'
  }
  else {
    # result of the gsva is then used to create a table
    design <- cbind(all = 1, condition = expression_classes != 'normal')
    fit <- limma::lmFit(result, design)
    fit <- limma::eBayes(fit)
    # consider adding p.value=p_value_cutoff
    result <- limma::topTable(fit, coef = 'condition', number = Inf, adjust.method = adjust.method)

    colnames(result)[colnames(result) == 'adj.P.Val'] <- 'fdr_q-val'
    colnames(result)[colnames(result) == 'logFC'] <- 'nes'
  }
  result
}
