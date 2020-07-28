###############################################################################
# Example production code for fitting latent class models with tree-structured
# shrinkage
# Zhenke Wu| zhenkewu@gmail.com
# First Version: July 10, 2020
###############################################################################

################################################################################
#                   Explore the tree encoding
#                   Jun 21st, 2020
#                   Zhenke Wu | zhenkewu@gmail.com
# NB:
# 1. Check the order of the internal and leaf nodes, which come first
# 2. ancestor (VI updates need this), leaf.descendants (VI update needs this),
#    descendants.
# 3. Visualize the tree (interactively so you can collapse upon clicking
#    on the nodes); also add estimates and some summaries to each leaf node.
# 4. adhoc collapsed groups - for computing separate LCM fits in each group
#    without letting them borrowing information. This means we need quick LCM
#    algorithm for each group.
# 5. need to create cross-validation folds within each leaf (but this is
#    infeasible for small leaves)
# 6. edge lengths - this is new to this kind of tree
# 7. Q: how to decide a reasonable set of initial values? random restart
#      to prevent stuck in local modes; criteria for stopping (nested -
#        hyperparameter update, VI parameter update, max??)
# 8. Q: Results to store [need to compile]
# 9. Q: how much difference does the edge-weight vs no edge weight make?
# 10.Q: How important is the tree information in shaping the posterior estimates;
#      what if the tree information is not that accurate? This is likely true
#      because the tree is estimated from ML with uncertainty.
# 11: some leaves may not have instances of one class! And some response probability
#     can be close to zero or 1 - this means collapsing?
#  12. may need to convert to igraph object. # leaves - 133;
#  13. some edge lengths are zeros. Just directly set these node-specific
#     additive component to be zero - because there is no information
#     for estimating s_u or gamma_u or alpha_u
#       # A tibble: 4 x 4
## # parent  node branch.length label
## # <int> <int>         <dbl> <chr>
## #   1    227    86             0 ST127_#_52
## # 2    259   119             0 ST963_#_4
## # 3    263   128             0 ST393_#_17
## # 4    231   232             0 NA
#   14. for initialization, need to use LCM with VI implementation (Grimmer, program
#      by my own or use existing package.)

###############################################################################
library(treeio)
library(ape)
library(stringr)
library(tidytree)
library(igraph)
library(compositions) # needed for alr? perhaps not; check.
library(foreach) # what is the difference with 'future'.
library(lotR)

#plotting:
library(magrittr)
library(RColorBrewer)
library(colorRamps)
library(dplyr)
library(ggtree)

rm(list=ls())
work_dir <- "/Users/zhenkewu/Dropbox/ZW/professional/=research_agenda/latent_variable/latent_class_model_with_phylogenetic_tree/tree_lcm_code"
res_dir <- "inst/example_figure"

thetree0 <- read.nexus(file.path(work_dir,"data/101816_STrep_ST4more_pav25Gub.filtered_polymorphic_sites_fixed.phylip_phyml_tree.nex"))
thetree0$tip.label <- str_extract_all(thetree0$tip.label,"ST[0-9]+$",simplify=TRUE)
#thetree$tip.label[109] <- "JJ1887"
thetree <- drop.tip(thetree0,109) # remove a reference tip used in constructing
# the phylogenetic tree using core genome. This ends up with 133 leaves.

# tree_df <- as.data.frame(cbind(thetree$edge,thetree$edge.length))
# colnames(tree_df) <- c("node_id1","node_id2","length")

# read in data (Mobile Genetic Elements):
dat_mge0 <- read.csv(file.path(work_dir,"data/HostElementData_hcl_mobile_20200424.csv"))
dat_mge0$human   <- c(rowSums(dat_mge0[,c("human_h1","human_h2","human_m1","human_ot")]))
dat_mge0$chicken <- rowSums(dat_mge0[,c("chicken_h1","chicken_h2","chicken_m1","chicken_ot")])
dat_mge0$turkey  <- rowSums(dat_mge0[,c("turkey_h1","turkey_h2","turkey_m1","turkey_ot")])

## check that the tips are a subset of the MLST in the dat_mge:
dat_mge <- dat_mge0

# the ones on the left are from data, the ones on the right are from the tree
# but how does the tree call the first node, in what order?
match_ind <- match(dat_mge$MLST,thetree$tip.label)
# not every one in the data is matched to the tree. Now let's see how many
# isolates are there in these unmatched STs; also some

sum(is.na(dat_mge0$MLST))      # 13 isolates have no sequence type info.
dat_mge$MLST[is.na(match_ind)] # we remove these isolates from main analysis because they are not mapped in the tree.
ct_ST <- table(dat_mge$MLST[!is.na(match_ind)])
hist(ct_ST)
sort(ct_ST) # 133 MLSTs that appear in the
# tree, a minimum of 4 and a maximum of 223 isolates per leaf.

ct_ST_human_or_not0 <- table(dat_mge[!is.na(match_ind),c("MLST","human")])

# reorder:
ct_ST_human_or_not <- ct_ST_human_or_not0[match(thetree$tip.label,
                                                rownames(ct_ST_human_or_not0)),]

tip_labels_without_ct <- thetree$tip.label
# match the leaves to the ST counts in the data:
tip_labels_with_ct <- paste(thetree$tip.label,
                            ct_ST[match(thetree$tip.label,names(ct_ST))],sep="_#_")
thetree$tip.label  <- tip_labels_with_ct

dat_mge$ct_MLST    <- dat_mge$MLST
dat_mge$ct_MLST[!is.na(match_ind)] <- paste(dat_mge$MLST[!is.na(match_ind)],
                                            ct_ST[match(dat_mge$MLST[!is.na(match_ind)],names(ct_ST))],sep="_#_")

# plot the tree (phylo object) for easy reference:
pdf(file.path(res_dir,"figs/thetree.pdf"),width=30,height=40)
plot(thetree,use.edge.length = TRUE,align.tip.label = TRUE)
nodelabels()
tiplabels()
#edgelabels()
dev.off()

# can we plot histogram of elements.
ind_EL <- grepl("^EL[0-9]+",colnames(dat_mge)) # get elements column.
EL_nm  <- colnames(dat_mge)[ind_EL]            # elements name.
# mean abundance by ST.
form_EL <- as.formula(paste0("cbind( ",paste(EL_nm,collapse=" , "),")~ MLST"))
abundance_by_ST0 <- aggregate(form_EL,dat_mge[!is.na(match_ind),],mean)
abundance_by_ST  <- as.matrix(abundance_by_ST0[,-1])

#rownames(abundance_by_ST) <- abundance_by_ST0$MLST
rownames(abundance_by_ST) <- paste(abundance_by_ST0$MLST,
                                   ct_ST[match(abundance_by_ST0$MLST,names(ct_ST))],
                                   sep="_#_")

# split into clades:
clade1 <- extract.clade(thetree,253)
edge_id1    <- match(clade1$edge.length,thetree$edge.length)
tip_id1     <- match(clade1$tip.label,thetree$tip.label)

clade2 <- extract.clade(thetree,209)
edge_id2    <- match(clade2$edge.length,thetree$edge.length)
tip_id2     <- match(clade2$tip.label,thetree$tip.label)

clade3 <- extract.clade(thetree,137)
edge_id3    <- match(clade3$edge.length,thetree$edge.length)
tip_id3     <- match(clade3$tip.label,thetree$tip.label)

edgecol <- rep("black",length(thetree$edge.length))
edgecol[edge_id1] <- "purple"
edgecol[edge_id2] <- "blue"
edgecol[edge_id3] <- "orange"

tipcol <- rep("black",length(thetree$edge.length))
tipcol[tip_id1] <- "purple"
tipcol[tip_id2] <- "blue"
tipcol[tip_id3] <- "orange"

plot(thetree,edge.color=edgecol,tip.color=tipcol)

library(phytools)
pdf(file.path(res_dir,"figs/thetree_EL_abundance.pdf"),width=30,height=40)
par(xpd=NA,oma=c(1,10,1,25),mar=c(1,3,1,3))
phylo.heatmap(thetree,abundance_by_ST,grid=TRUE,fsize = 1.5)
#split=c(0.7,0.3))

denom <- 2
rect(xleft=2,
     ybottom=seq(0,1,length=nrow(abundance_by_ST))-.5/(nrow(abundance_by_ST)),
     xright =2+ct_ST_human_or_not[,1]/max(ct_ST)/denom,
     ytop   =seq(0,1,length=nrow(abundance_by_ST))+.5/(nrow(abundance_by_ST)),
     col="gray",lwd=5,border=NA)
rect(xleft=2+ct_ST_human_or_not[,1]/max(ct_ST)/denom,
     ybottom=seq(0,1,length=nrow(abundance_by_ST))-.5/(nrow(abundance_by_ST)),
     xright =2+ct_ST_human_or_not[,1]/max(ct_ST)/denom+ct_ST_human_or_not[,2]/max(ct_ST)/denom,
     ytop   =seq(0,1,length=nrow(abundance_by_ST))+.5/(nrow(abundance_by_ST)),
     col="dodgerblue2",lwd=5,border=NA)
arrows(x0=2,
       y0=-0.01,
       x1 =2+(50+max(rowSums(ct_ST_human_or_not)))/max(ct_ST)/denom,
       y1   =-0.01,
       col="gray")
arrows(x0=2,
       y0=1.01,
       x1 =2+(50+max(rowSums(ct_ST_human_or_not)))/max(ct_ST)/denom,
       y1   =1.01,
       col="gray")
points(2+seq(0,250,length=5)/max(ct_ST)/denom,rep(-0.01,5),
       pch="|",col="gray")

text(2+seq(0,250,by=50)/max(ct_ST)/denom,rep(-0.01,5),seq(0,250,by=50),cex=2)
text(2+seq(0,250,by=50)/max(ct_ST)/denom,rep(1.01,5),seq(0,250,by=50),cex=2)
text(2+150/max(ct_ST)/denom,-0.03,"# of isolates",cex=3)
text(2+150/max(ct_ST)/denom,1.03,"# of isolates",cex=3)
segments(x0=2+seq(0,250,by=50)/max(ct_ST)/denom,
         y0=-0.01,
         x1=2+seq(0,250,by=50)/max(ct_ST)/denom,
         y1=1.01,
         col="gray")
dev.off()

##
## exploring the lengths:
## potentially remove ST770 (see how tau distance is relevant here?)
summary(thetree$edge.length)

which(as_tibble(thetree)$branch.length>0.14)
#[1] 133 135
as_tibble(thetree)[c(133,135),]
# A tibble: 2 x 4
#parent  node branch.length label
#<int> <int>         <dbl> <chr>
#  1    134   133         0.147 ST770_#_4
#2    134   135         0.147 NA
which(as_tibble(thetree)$branch.length==0)
#[1]  86 119 128 232
as_tibble(thetree)[c(86, 119, 128, 232),]

# other functions can find siblings, offsprings, etc.

# how to modify edge, node attributes? (igraph seems to have the attributes)
#?? currently don't know how to do this. This is important because
# 1) want to specify the ancestors.children etc.
# 2) want to specify the random variables in the tree-shrinkage prior; may need
#  multiple of them for different classes and dimensions in latent class models
# 3) need this to store estimates of the quantities in 2)

thetree_igraph <- as.igraph(thetree) # 133 tips and 132 internal nodes.
# This tree has a tip named JJ1887 which is a reference ST; which we removed.
E(thetree_igraph)$weight <- thetree$edge.length

# find ancestor (it is easy to do so using the igraph package.)
all_vertex_id <- sort(unique(unlist(as_tibble(thetree)[,c("parent","node")])))
leaf_id <- as_tibble(thetree)$node[!is.na(as_tibble(thetree)$label)]
ancestor_list <- descendent_list <-
  leaf_descendent_list <- vector("list",length = length(all_vertex_id))
for (i in all_vertex_id){
  ancestor_list[[i]] <- unique(unlist(ancestor(as_tibble(thetree),i)[,c("parent","node")]))
  descendent_list[[i]] <- getDescendants(thetree,i)
  leaf_descendent_list[[i]] <- intersect(leaf_id,descendent_list[[i]])
}

## the labels of the tip labels are:
thetree$tip.label

nodes <- names(igraph::V(thetree_igraph))

leaves <- names(igraph::V(thetree_igraph)[igraph::degree(thetree_igraph, mode = "out") == 0])
nodes <- c(nodes[!(nodes %in% leaves)], leaves)

library(ggtree)
library(ggplot2)
library(ggraph)

pdf(file.path(res_dir,"figs/thetree_ggtree.pdf"),width=30,height=40)
ggtree(thetree_igraph, ladderize = F) +
  geom_tiplab(geom = "label") + geom_nodelab(geom = "label") +
  theme(plot.margin = unit(c(0, 1.5, 0, 0.2), "cm")) +
  coord_cartesian(clip = "off") +
  scale_y_reverse()
dev.off()

# another way of visualization.
ggraph(thetree_igraph, layout = 'dendrogram', circular = FALSE) +
  geom_edge_diagonal() +
  geom_node_point() +
  theme_void()


# specify the groups of hyperparameters (just example, two levels):
# V(thetree_igraph)$levels <- rep(1,length(V(thetree_igraph)))
# V(thetree_igraph)$levels[match(names(igraph::V(thetree_igraph)[igraph::degree(thetree_igraph, mode = "out") == 0]),
#                               names(igraph::V(thetree_igraph)))] <- 2
# V(thetree_igraph)$levels[1] <- 3 # separate rootnode.

#V(thetree_igraph)$levels <- rep(2,length(V(thetree_igraph)))
# V(thetree_igraph)$levels[c(1)] <- 1

# # seven different
# V(thetree_igraph)$levels <- rep(7,length(V(thetree_igraph)))
# V(thetree_igraph)$levels[c(1)] <- 1
# V(thetree_igraph)$levels[c(2)] <- 2
# V(thetree_igraph)$levels[c(3)] <- 3
# V(thetree_igraph)$levels[c(4)] <- 4
# V(thetree_igraph)$levels[c(76)] <- 5
# V(thetree_igraph)$levels[c(120)] <- 6

# overall LCM:
# V(thetree_igraph)$levels       <- rep(2,length(V(thetree_igraph)))
# V(thetree_igraph)$levels[c(1)] <- 1

#V(thetree_igraph)$levels <- rep(1,length(V(thetree_igraph)))

# # generate a designed tree:
# dat_thetree <- design_tree(dat_mge[!is.na(match_ind),ind_EL], # these are not ordered yet.
#                            dat_mge[!is.na(match_ind),"ct_MLST"],
#                            thetree_igraph, weighted_edge = TRUE) # <-- this is actually not used in the lcm_tree.
#


# Y        <- dat_mge[!is.na(match_ind),ind_EL]   <---- real data

# works, our works best!
# tau <- c(0.5, 0.3, 0.2)
# theta <- rbind(rep(rep(c(0.8, 0.2), each = 1),9),
#                rep(rep(c(0.2, 0.8), each = 1),9),
#                rep(rep(c(0.9, 0.9), each = 1),9))
#
# # #
# # # SIMULATION
# # #
# tau   <- c(0.5, 0.3, 0.15,0.05)
# # theta <- rbind(rep(rep(c(0.01, 0.01), each = 1),9),
# #                rep(rep(c(0.01, 0.99), each = 1),9),
# #                rep(rep(c(0.99, 0.01), each = 1),9),
# #                rep(rep(c(0.5, 0.5), each = 1),9)) # works
#
#
# theta <- rbind(rep(rep(c(0.01, 0.01), each = 1),9),
#                rep(rep(c(0.01, 0.99), each = 1),9),
#                rep(rep(c(0.99, 0.01), each = 1),9),
#                c(0.1, 0.7,
#                  0.1, 0.7,
#                  0.1, 0.7,
#                  0.1, 0.7,
#                  0.1, 0.7,
#                  0.1, 0.7,
#                  0.1, 0.7,
#                  0.1, 0.7,
#                  0.1, 0.7))
# #rep(rep(c(0.4, 0.7), each = 1),9))
#
# image(t(theta))
# #                #c(rep(rep(c(0.05, 0.95), each = 1),9)[-18],0.1))
# Y <- BayesLCA::rlca(2663, itemprob = theta, classprob = tau)
#
# K     <- length(tau)
# p <- length(V(mytree))
# J <- ncol(Y)
#
# # #works well: well separated.
# # tau <- c(0.5, 0.3, 0.2)
# # theta <- rbind(rep(rep(c(0.8, 0.2), each = 1),9), rep(rep(c(0.2, 0.8), each = 1),9),
# #                rep(rep(c(0.9, 0.1,0.1,0.1,0.1,0.9,0.1,0.1,0.9), each = 1)))
# # Y <- BayesLCA::rlca(2663, itemprob = theta, classprob = tau)
#
# ## not well separated.
# # tau <- c(0.5, 0.3, 0.2)
# # theta <- rbind(rep(rep(c(0.8, 0.2), each = 1),9), rep(rep(c(0.2, 0.8), each = 1),9),
# #                rep(rep(c(0.05, 0.95,0.5), each = 1),6))
# # Y <- BayesLCA::rlca(2663, itemprob = theta, classprob = tau)
#
# dsgn0    <- design_tree(Y,outcomes,mytree)
# outcomes <- dat_mge[!is.na(match_ind),"ct_MLST"]
# mytree   <- thetree_igraph
#
# #set.seed(10) # this fix the initial mu_gamma mualpha if NULL; this also matters, because the initialization may cause some
# # blca item probabilities to be zero or one, throwing the error log(1-p) NaN.
# nrestarts <- 1
# # doParallel::registerDoParallel(cores = nrestarts)
# log_dir <- "restart_logs"
# dir.create(log_dir)
# # # mod0     <- lcm_tree(Y,outcomes,mytree,
# # #                      rootnode      = "Node1", # <-- may be redundant?
# # #                      weighted_edge = FALSE,
# # #                      hyper_fixed   = list(K=3,a=c(1,1,1,1,1,1,1),b=c(1,1,1,1,1,1,1),tau_update_levels = c(1:7),s_u_zeroset = (1:265)[-c(1,2,3,4,76,120)],s_u_oneset = c(1,2,3,4,76,120)),
# # #                      #hyper_fixed   = list(K=3,a=c(100000,1),b=c(1,1),tau_update_levels = c(1,2),s_u_zeroset = (1:265)[-1],s_u_oneset = c(1)),
# # #                      #hyper_fixed   = list(K=3,a=c(1,1),b=c(1,1),tau_update_levels = c(1,2),s_u_zeroset = (1:265)[-c(1,2,3,4,76,120)]),
# # #                      #hyper_fixed   = list(K=3,a=c(1,1),b=c(1,1),s_u_zeroset = 2:265),
# # #                      #hyperparams_init = list(tau_1=c(2.25^2,1e-10),tau_2=c(2.25^2,1e-10)),
# # #                      # hyper_fixed   = list(K=3,a=c(999,.0001),b=c(.0001,9999),tau_update_levels = 1),
# # #                      hyperparams_init = list(tau_1=rep(10,7),tau_2=rep(10,7)),
# # #                      nrestarts     = nrestarts,
# # #                      print_freq = 1,update_hyper_freq = 5, max_iter = 500,
# # #                      allow_continue = FALSE)#,
# # # #log_restarts = TRUE, log_dir = log_dir) # logs.
# #
# #
# #set.seed(345083)
# mod0     <- lcm_tree(Y,outcomes,mytree,
#                      rootnode      = "Node1", # <-- may be redundant?
#                      weighted_edge = FALSE,
#                      hyper_fixed   = list(K=K,a=c(1),b=c(1),tau_update_levels = c(1),
#                                           #s_u_zeroset = NULL,s_u_oneset = NULL),
#                                           #s_u_zeroset = 999,s_u_oneset = 999),
#                                           s_u_zeroset = (1:265)[-c(1)],s_u_oneset = c(1)),
#                      #s_u_zeroset = (1:265)[-c(1,2,3,4,76,120)],s_u_oneset = c(1,2,3,4,76,120)),
#                      hyperparams_init = list(tau_1=matrix(9/4,nrow=1,ncol=K-1),tau_2=array(9/4,c(1,J,K))),
#                      # hyperparams_init = list(tau_1=c(mean(c(alr(tau))^2)),tau_2=c(mean(logit(t(theta))^2))),
#                      # vi_params_init = list(mu_alpha = split_along_dim(t(replicate(p,alr(tau))),1),
#                      #                        mu_gamma = split_along_dim(array(logit(t(theta)),c(J,K,p)),3)
#                      #                        ),
#                      nrestarts     = nrestarts,
#                      # random_init   = FALSE,
#                      # random_init_vals = list(tau_lims = c(0.5, 1.5),
#                      #                         psi_sd_frac = 0.5,
#                      #                         phi_sd_frac = 0.5,
#                      #                         mu_gamma_sd_frac = .5,
#                      #                         mu_alpha_sd_frac = .5,
#                      #                         u_sd_frac = 0.2),
#                      print_freq = 10,update_hyper_freq = 10, max_iter = 5000,
#                      allow_continue = FALSE)
#
# # ,log_restarts = TRUE, log_dir = log_dir
# #                      )
# # unlink(log_dir, recursive = T)
# # closeAllConnections()
#
# # mod0     <- lcm_tree(Y,outcomes,mytree,
# #                      rootnode      = "Node1", # <-- may be redundant?
# #                      weighted_edge = FALSE,
# #                      hyper_fixed   = list(K=K,a=c(1,1),b=c(1,1),tau_update_levels = c(99),
# #                                           #s_u_zeroset = 999,s_u_oneset = 999),
# #                                           s_u_zeroset = (1:265)[-c(1)],s_u_oneset = c(1)),
# #                      #s_u_zeroset = (1:265)[-c(1,2,3,4,76,120)],s_u_oneset = c(1,2,3,4,76,120)),
# #                      hyperparams_init = list(tau_1=c(mean(c(alr(tau))^2),9/4),tau_2=c(mean(logit(t(theta))^2),9/4)),
# #                      # vi_params_init = list(mu_alpha = split_along_dim(t(replicate(p,alr(tau))),1),
# #                      #                       mu_gamma = split_along_dim(array(logit(t(theta)),c(J,K,p)),3)
# #                      #                       ),
# #                      nrestarts     = nrestarts,
# #                      # random_init   = FALSE,
# #                      # random_init_vals = list(tau_lims = c(0.5, 1.5),
# #                      #                         psi_sd_frac = 0.5,
# #                      #                         phi_sd_frac = 0.5,
# #                      #                         mu_gamma_sd_frac = .5,
# #                      #                         mu_alpha_sd_frac = .5,
# #                      #                         u_sd_frac = 0.2),
# #                      print_freq = 1,update_hyper_freq = 10, max_iter = 5000,
# #                      allow_continue = FALSE)#,
# #log_restarts = TRUE, log_dir = log_dir) # logs.
# #
# #
# # unlink(log_dir, recursive = T)
# # closeAllConnections()
#
# # summarize posterior results:
# plot(mod0)
# # pdf("0712_non_weighted.pdf",width=30,height=10)
# # dev.off()
#
# #image(mod0$prob_est$pi_collapsed)
# mod0$prob_est$pi_collapsed
#
#
# mexpit(c(mod0$prob_est$eta_est[1,],0))
#
# colMeans(mod0$mod$vi_params$rmat)
#
#
# rowMeans(apply(cbind(MASS::mvrnorm(10000,mod0$prob_est$eta_est[1,],mod0$prob_est$eta_var_est[,,1]),0),1,mexpit))
#
#
# #colnames(dsgn0$Y) <- c("EL2","EL3", "EL40", "EL35", "EL44", "EL45", "EL46", "EL12", "EL18","EL19", "EL36", "EL37", "EL38", "EL39" ,"EL41", "EL42", "EL43", "EL50")
# # blca:
# # VB:
# fullLCM_vb <- BayesLCA::blca(dsgn0$Y,K,method="vb",verbose=FALSE)
# fullLCM_vb
# # EM:
# fullLCM_em <- BayesLCA::blca(dsgn0$Y,K,method="em",verbose=FALSE)
# fullLCM_em
#
# ## poLCA:
# # form = as.formula("cbind(EL2,EL3, EL40, EL35, EL44, EL45, EL46, EL12, EL18,
# #                   EL19, EL36, EL37, EL38, EL39 ,EL41, EL42, EL43, EL50)~1")
# # for (j in 1:ncol(dsgn0$Y)){dsgn0$Y[,j] <- as.integer(as.factor(dsgn0$Y[,j]))}
# # res2 <- poLCA::poLCA(form, as.data.frame(dsgn0$Y),nclass=K)
#
# proposed <- mod0$prob_est$theta_collapsed#[,c(1,3,2)]
#
# par(mfcol=c(3,3))
# image(t(theta),main="truth")
# # differences
# image((proposed-t(theta)),main="diff")
# hist((proposed-t(theta)),xlim=c(-1,1))
#
# image(proposed,main="proposed")
# image(t(fullLCM_vb$itemprob),main="BayesLCA package: vb")
# image(t(fullLCM_em$itemprob),main="BayesLCA package: em")
#
# # pi estimate:
# barplot(rbind(tau,
#               mod0$prob_est$pi_collapsed,
#               fullLCM_vb$classprob,
#               fullLCM_em$classprob),main="class assignment", beside=TRUE)#,legend.text = c("truth","proposed","vb","em"))
#
#
#
#
# # # class assignments:
# # table(apply(mod0$mod$vi_params$rmat,1,which.max))
# # table(apply(fullLCM_vb$Z,1,which.max))
# # table(apply(fullLCM_em$Z,1,which.max))
# #
# # image(fullLCM_vb$Z[1:10,])
# # image(fullLCM_em$Z[1:10,])
# # image(mod0$mod$vi_params$rmat[1:10,])
# #
# #
#
#
# #mod1 <- continue_lcm_tree(mod0$old_mod,print_freq = 5,update_hyper_freq = 20)
# #
# # # with Random Restarts, with edge weights:
# # nrestarts <- 3
# # doParallel::registerDoParallel(cores = nrestarts)
# # set.seed(345083)
# # log_dir <- "restart_logs"
# # dir.create(log_dir)
# #
# # mod1     <- lcm_tree(Y,outcomes,mytree,
# #                      rootnode = "Node1", # <-- may be redundant?
# #                      weighted_edge = TRUE,
# #                      nrestarts = nrestarts,
# #                      print_freq = 1,update_hyper_freq = 10,
# #                      log_restarts = TRUE, log_dir = log_dir)
# # unlink(log_dir, recursive = T)
# # # Print the final ELBO for all restarts
# # print(c(mod1$mod$hyperparams$ELBO,
# # sapply(mod1$mod_restarts, function(mod) mod$hyperparams$ELBO)), digits = 8)
# #
# #
# # # use the existing mod1 fit; so the initialization is better:
# # nrestarts <- 1
# # # doParallel::registerDoParallel(cores = nrestarts)
# # # set.seed(345083)
# # # log_dir <- "restart_logs"
# # # dir.create(log_dir)
# # mod2     <- lcm_tree(Y,outcomes,mytree,
# #                      rootnode = "Node1", # <-- may be redundant?
# #                      weighted_edge = TRUE,
# #                      vi_params_init = mod1$mod$vi_params,
# #                      hyperparams_init = mod1$mod$hyperparams,
# #                      nrestarts = nrestarts,
# #                      print_freq = 1,update_hyper_freq = 5,max_iter = 11#,
# #                      #log_restarts = TRUE, log_dir = log_dir
# #                      )
# # # unlink(log_dir, recursive = T)
