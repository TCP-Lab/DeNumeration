#' Convert a DESeq2 "result" object to an enumeration object
#' 
#' @param data A DESeq2 object from the 'results' function
#' @param alpha The alpha p-value threshold above which a gene is considered
#'    to be not differentially expressed. E.g. 0.05 would cause any gene with a
#'    p-value of 0.06 to be not differentially expressed.
#' @param fold_change The absolute fold change value a gene must have to be
#'    considered to as differentially expressed. E.g. 2 would consider a gene
#'    with fold change 1 to be not deferentially expressed.
#' @param zero_lab Label to give in the enumeration to non-differentially expressed genes.
#' @param pos_lab Label to give in the enumeration to upregulated genes.
#' @param neg_lab Label to give in the enumeration to downregulated expressed genes.
#' @param colname The column name to give to the variable with the labels in the resulting frame
#' 
#' @return A data.frame with one row per input gene and a column with the proper
#'    enumeration labels based on the filter in the column "colname", by default
#'    "class".
#' 
#' @author Hedmad
#'    
#' @export
to_enumeration_vector <- function(
        data, alpha = 0.01, fold_change = 1,
        zero_lab = "zero", pos_lab = "positive", neg_lab = "negative",
        colname = "class"
) {
    # This makes just a single dframe with one col: the enumeration
    new_data <- data[, "baseMean", drop=FALSE]
    # This is just to not lose the rownames
    
    res <- factor(
        rep(zero_lab, nrow(data)),
        ordered = TRUE,
        # since ordered is true, the order of the labels here is important
        levels = c(pos_lab, zero_lab, neg_lab)
    )
    
    res[data$log2FoldChange > fold_change & data$padj < alpha] <- pos_lab
    res[data$log2FoldChange < - fold_change & data$padj < alpha] <- neg_lab
    
    new_data[[colname]] <- res
    new_data$baseMean <- NULL
    
    return(new_data)
}

#' Apply `fun` to all results out of a deseq object.
#' 
#' @param deseq_object A Deseq Dataset with results (after the DESeq call)
#' @param results A list of result names (as from `DESeq2::listNames`), for each
#'    result on which the function should be applied to.
#' @param fun A function to be applied
#' @param ... Further arguments to be passed to `fun`
#' 
#' @returns A list with one frame per item in `results`, with the output of calling
#'    `fun` on that result.
#'    
#' @export
apply_to_results <- function(deseq_object, results, fun, ...) {
    res <- list()
    for (result in results) {
        print(result)
        res[[result]] <- fun(DESeq2::results(deseq_object, name = result), ...)
    }
    
    res
}


#' Plot an enumeration frame from a list of enumeration data
#' 
#' This function takes care of computing frequency tables from input enumerations
#' and plotting the resulting frequency graph.
#' 
#' @param enum_data A list with one slot per enumeration, with items returned by
#'    `to_enumeration_vector`.
#' @param order_by Order of the rows in the resulting plot. Either "frequency", to
#'    sort by frequency, or "label" to sort by categorical labels (up, zero, down).
#' @param zero_lab Label in the enumeration for non-differentially expressed genes.
#' @param pos_lab Label in the enumeration for upregulated genes.
#' @param neg_lab Label in the enumeration for downregulated expressed genes.
#' @param title Title of the plot.
#' @param top_n If not NULL, an integer. Keep only `top_n` rows by frequency.
#' @param exclude_all_negative If TRUE, the case where all variables are zero is
#'    omitted from the plot.
#' @param category_renames A list of `old_label = new_label` with the new names
#'    of the enumeration categories. The old labels are the column names of the
#'    input enumeration data.
#' 
#' @param positive_symbol The symbol in the plot to give to the upregulated case.
#' @param negative_symbol The symbol in the plot to give to the negative case.
#' @param positive_color The color in the plot to give to the upregulated case.
#' @param negative_color The color in the plot to give to the negative case.
#' @param symbol_size The size of the positive and negative legend symbols.
#' @param symbol_fontface The font face of the positive and negative symbols.
#' 
#' @param labels_y_nudge By default, the symbols are centered on the correct
#'    position. However, based on the used font, the symbols might be misaligned
#'    in the Y axis. The "nudge" to give to each symbol is plot size and
#'    font dependent, so it's difficult to compute beforehand.
#'    Tweak this parameter to move the symbols up (e.g. +0.02) or down (e.g. -0.02)
#'    so that they are perfectly aligned for your specific usecase.
#' 
#' @return A `ggproto` object with the resulting plot.
#' 
#' @author Hedmad
#' 
#' @export
plot_enumeration_frame <- function(
        enum_data, order_by = c("frequency", "label"),
        zero_lab = "zero", pos_lab = "positive", neg_lab = "negative",
        title = NULL, top_n = NULL, labels_y_nudge=0,
        positive_symbol = "+", negative_symbol = "-",
        symbol_size = 10,
        positive_color = "red", negative_color = "blue",
        symbol_fontface = "bold",
        exclude_all_negative = FALSE,
        category_renames = NULL
) {
    # Argument pre-parsing
    order_by <- order_by[1]
    # we need to preserve the ordered factors, so I take out one of them here
    ord_factor_levels <- enum_data[[1]] %>% levels()
    
    # create the frequency matrix and recast it as d-frame
    frequencies <- enum_data %>% as.data.frame() %>% stats::ftable() %>% as.data.frame()
    all_factors <- utils::head(colnames(frequencies), -1)
    
    # return to ordered factors
    for (lab in all_factors) {
        frequencies[[lab]] <- factor(frequencies[[lab]], ordered = TRUE, levels = ord_factor_levels)
    }
    
    # Delete the "all negatives" case, if asked
    if (exclude_all_negative) {
        truth_key <- apply(frequencies[, all_factors], 1, \(x) {all(x == zero_lab)})
        # Stop if we detect more than one all negative case (something is terribly wrong)
        assertthat::assert_that(sum(truth_key) == 1, msg = "More than one all zero case found! Aborting.")
        frequencies <- frequencies[!truth_key, ] # we negate, as we want to exclude just that case
    }
    
    # add a static label of the combination of factors
    frequencies$label <- if (length(all_factors) == 1) {
        frequencies[[all_factors[1]]]
    } else {
        apply( frequencies[, all_factors], 1, paste, collapse="_")
    }
    
    
    # If we just want to show the top N, we dramatically subset here.
    # Drumroll...
    if (!is.null(top_n)) {
        frequencies <- frequencies %>% dplyr::slice_max(Freq, n = top_n)
        print(frequencies)
    }
    
    # Now we need to sort the frame by both "freq" and by the labels
    frequencies$order <- if (order_by == "frequency") {
        order(frequencies$Freq)
    } else if (order_by == "label") {
        # factors are a bit less easy
        do.call(order, as.list(frequencies[, all_factors, drop = FALSE]))
    } else {
        stop(paste0("Invalid argument 'order': '", order_by, "'. Possibilities: 'frequency', 'label'"))
    }
    
    bar_plot <- frequencies %>% ggplot2::ggplot(ggplot2::aes(x = Freq, y = label)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::ylab("") + ggplot2::xlab("Gene Count") +
        ggplot2::scale_y_discrete(
            limits = frequencies$label[frequencies$order],
            labels=NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::geom_text(stat="identity", ggplot2::aes(label = Freq), hjust = -0.1) +
        ggplot2::coord_cartesian(xlim = c(0, max(frequencies$Freq) * 1.10))+
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.ticks.length.y = ggplot2::unit(0, "pt")
        )
    
    melted_legend <- reshape2::melt(
        frequencies[, c("label", all_factors)],
        id.vars = "label",
        measure_vars=all_factors
    )
    
    # Convert from the three labels to a point that we can display
    melted_legend$pos <- ""
    melted_legend$pos[melted_legend$value == pos_lab] <- positive_symbol
    melted_legend$neg <- ""
    melted_legend$neg[melted_legend$value == neg_lab] <- negative_symbol
    
    ## This is to appease the NOTEs section of devtools::check()
    # See https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    # Fuck you CRAN, and fuck you R.
    
    Freq <- label <- margin <- neg <- pos <- variable <- NULL
    
    labels <- if (!is.null(category_renames)) {
        new_labs <- c()
        for (label in levels(melted_legend$variable)) {
            if (label %in% names(category_renames)) {
                new_labs <- c(new_labs, category_renames[[label]])
            } else {
                new_labs <- c(new_labs, label)
            }
        }
        new_labs
    } else {
        levels(melted_legend$variable)
    }
    
    label_plot <- melted_legend %>% ggplot2::ggplot(ggplot2::aes(x = variable, y = label)) +
        ggplot2::geom_text(
            ggplot2::aes(label=pos),
            size = symbol_size, colour = negative_color, fontface=symbol_fontface,
            vjust=0.5, hjust=0.5,
            nudge_y=labels_y_nudge
        ) +
        ggplot2::geom_text(
            ggplot2::aes(label=neg),
            size = symbol_size, colour = positive_color, fontface=symbol_fontface,
            vjust=0.5, hjust=0.5,
            nudge_y=labels_y_nudge
        ) +
        ggplot2::scale_x_discrete(
            labels=labels
        ) +
        ggplot2::scale_y_discrete(
            limits=melted_legend$label[frequencies$order],
            labels=NULL
        ) +
        ggplot2::ylab(NULL) +
        ggplot2::xlab(NULL) +
        ggplot2::theme_minimal()+
        ggplot2::theme(
            plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"),
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)
        )
    
    # We have the plots, we just have to tape them together.
    label_plot + patchwork::plot_spacer() + bar_plot +
        patchwork::plot_layout(
            ncol=3, nrow=1, widths=c(0.2, -0.04, 0.9),
            guides = "collect"
        ) +
        patchwork::plot_annotation(
            title = title
        )
}


#' D-enumerate a DESeq2 Dataset object with included results to an enumerator
#' 
#' With this function you can quickly move from a DESeq2 dataset with
#' computed results (i.e. after calling `DESeq()`) and a complicated
#' experimental design to something simpler that can be better observed
#' with functions such as `plot_enumeration_frame`
#' 
#' @param computed_deseq_object DESeq2 Dataset object with included results
#' @param results Either a vector of named results (as in from `DESeq2::resultsNames`)
#'    or the special cases "nointercept" for all results without the intercept
#'    and "all" for all results from the object.
#' @param new_labels A list of `(old_label = new_label)` to replace the often
#'    lengthy result labels in the DESeq2 object with more understandable ones.
#'    The old labels are the same ones as specified in the `results` parameter.
#' @param ... Further arguments to be passed to `to_enumeration_vector`.
#'    Specifically, consider tweaking `alpha` and `fold_change`.
#' 
#' @returns A list of with one slot per result, with the computed enumerations.
#' 
#' @export
denumerate <- function(
        computed_deseq_object, results = "nointercept",
        new_labels = NULL, ...
) {
    # Parse the "results" parameter
    possible_results <- DESeq2::resultsNames(computed_deseq_object)
    if (results == "nointercept") {
        assertthat::assert_that(possible_results[1] == "Intercept")
        results <- possible_results[-1]
    } else if (results == "all") {
        results <- possible_results
    } else {
        assertthat::assert_that(
            all(results %in% possible_results)
        )
    }
    
    # Execute the enumerations and rename the internal frames, so we can merge
    enumerations <- apply_to_results(
        computed_deseq_object,
        results,
        to_enumeration_vector,
        ...
    )
    # I merge just with cbind(), so I want to be extra sure the rows are all
    # in the same order
    test_row_order <- enumerations[[1]] %>% rownames()
    for (name in names(enumerations)) {
        colnames(enumerations[[name]]) <- name
        assertthat::are_equal(rownames(enumerations[[name]]), test_row_order)
    }
    enum_frame <- Reduce(cbind, enumerations)
    
    if (!is.null(new_labels)) {
        enum_frame <- tibble::as_tibble(enum_frame) %>% dplyr::rename_with(\(x) {
            .fn <- \(x) {
                if (x %in% names(new_labels)) {
                    return(new_labels[[x]])
                } else {
                    return(x)
                }
            }
            return(sapply(x, .fn))
        }) %>% as.data.frame()
    }
    
    # Restore the row names
    row.names(enum_frame) <- test_row_order
    
    enum_frame
}
