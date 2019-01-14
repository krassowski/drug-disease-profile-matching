gsva.with_probabilities <- function(
  expression, expression_classes, gene_sets, procedure,
  method='gsva', mx.diff=TRUE, permutations=1000, adjust.method='fdr',
  ...
) {
  stopifnot(procedure %in% c('gene_permutation', 'bayes'))

  design <- cbind(condition = expression_classes != 'normal')
  phenoData <- Biobase::AnnotatedDataFrame(data = as.data.frame(as.table(design)))
  row.names(phenoData) <- colnames(expression)
  expression_set <- Biobase::ExpressionSet(assayData = data.matrix(expression), phenoData = phenoData)

  # transform to named list from GeneSetCollection class object
  geneSets <- GSEABase::geneIds(gene_sets)

  expressions <- Biobase::exprs(expression_set)
  genesInExpressionData <- rownames(expressions)
  # also, see GSVA::computeGeneSetsOverlap()
  overlaps <- sapply(geneSets, function(genes) genes[genes %in% genesInExpressionData])
  subset <- overlaps[sapply(overlaps, function(genes) length(genes) > 1)]

  result <- GSVA::gsva(expressions, subset, method = method, verbose = F, parallel.sz = 1, mx.diff = mx.diff, ...)

  if (procedure == 'gene_permutation')
  {
    genes_order <- rownames(expressions)

    result <- as.data.frame(result)
    result$difference <- result$condition - result$control

    n_permutations <- permutations
    permutations <- c(1:n_permutations)
    permutation_effect_sizes <- lapply(permutations, function(i) {{
      rownames(expressions) <- sample(genes_order)  # a random permutation
      random_result <- GSVA::gsva(expressions, subset, method = method, verbose = F, parallel.sz = 1, mx.diff = mx.diff, ...)
      random_effect_size <- random_result[,'condition'] - random_result[,'control']
      random_effect_size
    }})

    random_effect_sizes <- as.data.frame(permutation_effect_sizes, col.names = permutations)

    result$p_value <- mapply(
      function(effect_size, gene_set_name) {{
        if (effect_size >= 0) {{
          is_random_permutation_more_extreme <- random_effect_sizes[gene_set_name,] > effect_size
        }}
        else {{
          is_random_permutation_more_extreme <- random_effect_sizes[gene_set_name,] < effect_size
        }}
        sum(is_random_permutation_more_extreme) / n_permutations
      }},
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
