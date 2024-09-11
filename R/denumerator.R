requireNamespace("ggplot2")
requireNamespace("assertthat")
requireNamespace("reshape2")
requireNamespace("patchwork")

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

apply_to_results <- function(deseq_object, results, fun, ..., pass_names = FALSE) {
    res <- list()
    for (result in results) {
        print(result)
        res[[result]] <- fun(DESeq2::results(deseq_object, name = result), ...)
    }
    
    res
}

plot_enumeration_frame <- function(
        enum_data, order_by = c("frequency", "label"),
        zero_lab = "zero", pos_lab = "positive", neg_lab = "negative",
        title = NULL, top_n = NULL, labels_y_nudge=0,
        positive_symbol = "+", negative_symbol = "-",
        symbol_size = 10,
        positive_color = "red", negative_color = "blue",
        symbol_fontface = "bold"
) {
    # Argument pre-parsing
    order_by <- order_by[1]
    # we need to preserve the ordered factors, so I take out one of them here
    ord_factor_levels <- enum_data[[1]] %>% levels()
    
    # create the frequency matrix and recast it as d-frame
    frequencies <- enum_data %>% as.data.frame() %>% ftable() %>% as.data.frame()
    all_factors <- head(colnames(frequencies), -1)
    
    # return to ordered factors
    for (lab in all_factors) {
        frequencies[[lab]] <- factor(frequencies[[lab]], ordered = TRUE, levels = ord_factor_levels)
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
        frequencies <- frequencies %>% slice_max(Freq, n = top_n)
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
    
    bar_plot <- frequencies %>% ggplot(aes(x = Freq, y = label)) +
        geom_bar(stat = "identity") +
        ylab("") + xlab("Gene Count") +
        scale_y_discrete(
            limits = frequencies$label[frequencies$order],
            labels=NULL
        ) +
        theme_bw() +
        geom_text(stat="identity", aes(label = Freq), hjust = -0.1) +
        coord_cartesian(xlim = c(0, max(frequencies$Freq) * 1.10))+
        theme_minimal() +
        theme(
            plot.margin = margin(0, 0, 0, 0, "pt"),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.length.y = unit(0, "pt")
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
    
    label_plot <- melted_legend %>% ggplot(aes(x = variable, y = label)) +
        geom_text(
            aes(label=pos),
            size = symbol_size, colour = negative_color, fontface=symbol_fontface,
            vjust=0.5, hjust=0.5,
            nudge_y=labels_y_nudge
        ) +
        geom_text(
            aes(label=neg),
            size = symbol_size, colour = positive_color, fontface=symbol_fontface,
            vjust=0.5, hjust=0.5,
            nudge_y=labels_y_nudge
        ) +
        scale_y_discrete(
            limits=melted_legend$label[frequencies$order],
            labels=NULL
        ) +
        ylab(NULL) +
        xlab(NULL) +
        theme_minimal()+
        theme(
            plot.margin = margin(0, 0, 0, 0, "pt"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        )
    
    label_plot + patchwork::plot_spacer() + bar_plot +
        patchwork::plot_layout(
            ncol=3, nrow=1, widths=c(0.2, -0.04, 0.9),
            guides = "collect"
        ) +
        patchwork::plot_annotation(
            title = title
        )
}

denumerate <- function(
        computed_deseq_object, results = "nointercept",
        new_labels = NULL,
        ...
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
        pass_names = TRUE
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
        enum_frame <- as.tibble(enum_frame) %>% rename_with(\(x) {
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
    
    plot_enumeration_frame(enum_frame, ...)
}
