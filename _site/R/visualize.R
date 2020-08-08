if(getRversion() >= "2.15.1") utils::globalVariables(c("Group"))

#' plots the groups discovered by lcm_tree on the tree for observations.
#'
#' @param x Output from `moretrees()`
#' An object of class "moretrees_result".
#' @param group.text.size Text size for the group labels
#' @param group.text.offset Offset of the group label from the
#' leaves of the tree
#' @param legend.text.size Text size for legend
#' @param layout Layout for the tree, most likely "rectangular" (the default)
#' or "slanted", but see the `layout` option of `ggtree()` for more
#' possibilities
#' @param horizontal If TRUE (the default), the tree will be plotted with
#' the root node at the top and all other nodes below. If FALSE, the tree
#' will be plotted with the root node to the left and all other nodes to
#' the right.
#' @param ... Not used.
#' @return A plot showing the groups discovered by MOReTreeS on the original
#' leaf_ids tree.
#' @examples
#' # See vignette
#'
#' @importFrom ggtree %<+%
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
plot.lcm_tree <- function(x,
                          group.text.size = 4,
                          group.text.offset = 0.1,
                          legend.text.size = 10,
                          layout = "rectangular",
                          horizontal = TRUE,
                          ...) {
# the following is from moretrees:
  # Get data.frame with groups info
  leaves <- names(igraph::V(x$mytree)[igraph::degree(x$mytree, mode = "out") == 0])
  groups.df <- data.frame(leaves = leaves, Group = as.factor(x$prob_est$group))
  cols_g <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(x$prob_est$group)))
  #cols_g <- RColorBrewer::brewer.pal(length(levels(groups.df$Group)), "Set3")

  # Plot
  p <- ggtree::ggtree(x$mytree, ladderize = FALSE, layout = layout)  %<+% groups.df
  if (horizontal) {
    p <- p + ggplot2::coord_flip() + ggplot2::scale_x_reverse()
    group.text.offset <- - group.text.offset
  }
  p <- p + ggtree::geom_tippoint(ggplot2::aes(color = Group),
                                 shape = 15, size = 4) +
    ggtree::geom_tiplab(ggplot2::aes(label = Group), size = group.text.size,
                        offset = group.text.offset,
                        vjust = 0.5,
                        hjust = 0.5) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = legend.text.size),
                   legend.title = ggplot2::element_text(size = legend.text.size))
  if (length(cols_g) == length(levels(groups.df$Group))) {
    p <- p + ggplot2::scale_color_manual(values = cols_g)
  }

  # Return
  return(p)
}
