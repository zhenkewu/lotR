if(getRversion() >= "2.15.1") utils::globalVariables(c("Group","Group_nm","name"))

#' plots the groups discovered by lcm_tree on the tree for observations.
#'
#' @param x Output from `moretrees()`
#' An object of class "moretrees_result".
#' @param layout layout for the tree, by default it is `"circular"`;
#' other layouts such as `"rectangular"` can run, but currently has issues.
#' @param colnames_offset_y default is `0`; see `gheatmap` from `ggtree`.
#' @param heatmap_width default is `1`; see `gheatmap` from `ggtree`.
#' @param ... Not used.
#' @return A plot showing the groups discovered by MOReTreeS on the original
#' leaf_ids tree.
#' @examples
#' # See vignette
#'
#' @import ggtree
#' @importFrom igraph E
#' @importFrom ape as.phylo
#' @importFrom stats median
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom grDevices colorRampPalette topo.colors
#' @importFrom RColorBrewer brewer.pal
#' @export
plot.lcm_tree <- function(x,
                          # group.text.size = 4,
                          # group.text.offset = 0.1,
                          # legend.text.size = 10,
                          layout = "circular",
                          colnames_offset_y=0,
                          heatmap_width = 1,
                          ...) {

  classprob <- data.frame(x$prob_est$pi)
  classprob$id <- unique(names(x$dsgn$leaf_ids))
  colnames(classprob) <- c(paste("class",1:x$mod$hyper_fixed$K,
                                 sep=""),"id")
  classprob <- classprob[,c(x$mod$hyper_fixed$K+1,1:x$mod$hyper_fixed$K)]
  classprob$Group <- as.factor(x$prob_est$group)
  classprob$Group_nm <- as.character(classprob$Group)
  col_g_uniq <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(x$prob_est$group)))
  rownames(classprob) <- classprob$leaf_nm
  classprob$name <- classprob$id

  thetree_here <- as.phylo(x$mytree) # add the edge weights to the transformed object;
  thetree_here$edge.length <- E(x$mytree)$weight
  thetree_edge_scaled <- thetree_here
  thetree_edge_scaled$edge.length <- thetree_here$edge.length*0.5/median(thetree_here$edge.length)
  p <- ggtree(thetree_edge_scaled,layout=layout) %<+% classprob +
    geom_tippoint(aes(color=Group))+
    scale_color_manual(values=col_g_uniq)+
    geom_tiplab2(aes(label=Group_nm,color=Group), align=T, linetype=NA,
                 size=3, offset=1, hjust=0.5)+
    geom_tiplab2(aes(label=name,color=Group), align=T, linetype=NA,
                 size=3, offset=7, hjust=0.5)

  heatmapData <- classprob[,1+(1:x$mod$hyper_fixed$K)]
  rownames(heatmapData) <- classprob$name

  heatmap.colours <- c("white","grey","seagreen3","darkgreen",
                       "green","brown","tan","red","orange",
                       "pink","magenta","purple","blue","skyblue3",
                       "blue","skyblue2")
   p <- gheatmap(p, heatmapData, offset = 12,color=NULL,
           colnames_position="top",
           colnames_angle=90, colnames_offset_y = 0,
           hjust=0, font.size=heatmap_width,width=heatmap_width)+
    scale_fill_gradientn(colours=topo.colors(100),na.value = "transparent",
                         breaks=c(0,0.5,1),labels=c(0,0.5,1),
                         limits=c(0,1))
   p

  # Return
  return(p)
}



###### @param group.text.size Text size for the group labels
##### @param group.text.offset Offset of the group label from the
##### leaves of the tree
##### @param legend.text.size Text size for legend
##### @param layout Layout for the tree, most likely "rectangular" (the default)
##### or "slanted", but see the `layout` option of `ggtree()` for more
##### possibilities
##### @param horizontal If TRUE (the default), the tree will be plotted with
##### the root node at the top and all other nodes below. If FALSE, the tree
##### will be plotted with the root node to the left and all other nodes to
##### the right.



# plot.lcm_tree <- function(x,
#                           group.text.size = 4,
#                           group.text.offset = 0.1,
#                           legend.text.size = 10,
#                           layout = "rectangular",
#                           horizontal = TRUE,
#                           ...) {
# # the following is from moretrees:
#   # Get data.frame with groups info
#   leaves <- names(igraph::V(x$mytree)[igraph::degree(x$mytree, mode = "out") == 0])
#   groups.df <- data.frame(leaves = leaves, Group = as.factor(x$prob_est$group))
#   cols_g <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(x$prob_est$group)))
#   #cols_g <- RColorBrewer::brewer.pal(length(levels(groups.df$Group)), "Set3")
#
#   # Plot
#   p <- ggtree::ggtree(x$mytree, ladderize = FALSE, layout = layout)  %<+% groups.df
#   if (horizontal) {
#     p <- p + ggplot2::coord_flip() + ggplot2::scale_x_reverse()
#     group.text.offset <- - group.text.offset
#   }
#   p <- p + ggtree::geom_tippoint(ggplot2::aes(color = Group),
#                                  shape = 15, size = 4) +
#     ggtree::geom_tiplab(ggplot2::aes(label = Group), size = group.text.size,
#                         offset = group.text.offset,
#                         vjust = 0.5,
#                         hjust = 0.5) +
#     ggplot2::theme(legend.text = ggplot2::element_text(size = legend.text.size),
#                    legend.title = ggplot2::element_text(size = legend.text.size))
#   if (length(cols_g) == length(levels(groups.df$Group))) {
#     p <- p + ggplot2::scale_color_manual(values = cols_g)
#   }
#
#   # Return
#   return(p)
# }




