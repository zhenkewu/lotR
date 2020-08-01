###############################################################################
# Example production code for fitting latent class models with tree-structured
# shrinkage
# Zhenke Wu| zhenkewu@gmail.com
# First Version: July 10, 2020
# 2nd version: July 31, 2020
###############################################################################
library(treeio)
library(ape)
library(stringr)
library(tidytree)
library(igraph)
library(foreach) # what is the difference with 'future'?
library(lotR)

#plotting:
library(magrittr)
library(RColorBrewer)
library(colorRamps)
library(dplyr)
library(ggtree)

rm(list=ls())
work_dir <- "/Users/zhenkewu/Dropbox/ZW/professional/=research_agenda/latent_variable/latent_class_model_with_phylogenetic_tree/tree_lcm_code"
res_dir  <- "inst/example_figure"

thetree0 <- read.nexus(file.path(work_dir,"data/101816_STrep_ST4more_pav25Gub.filtered_polymorphic_sites_fixed.phylip_phyml_tree.nex"))
thetree0$tip.label <- str_extract_all(thetree0$tip.label,"ST[0-9]+$",simplify=TRUE)
#thetree$tip.label[109] <- "JJ1887"
thetree <- drop.tip(thetree0,109) # remove a reference tip used in constructing
# the phylogenetic tree using core genome; we now have pL = 133 leaves.

# tree_df <- as.data.frame(cbind(thetree$edge,thetree$edge.length))
# colnames(tree_df) <- c("node_id1","node_id2","length")

# read in data (Mobile Genetic Elements):
dat_mge0         <- read.csv(file.path(work_dir,"data/HostElementData_hcl_mobile_20200424.csv"))
dat_mge0$human   <- c(rowSums(dat_mge0[,c("human_h1","human_h2","human_m1","human_ot")]))
dat_mge0$chicken <- rowSums(dat_mge0[,c("chicken_h1","chicken_h2","chicken_m1","chicken_ot")])
dat_mge0$turkey  <- rowSums(dat_mge0[,c("turkey_h1","turkey_h2","turkey_m1","turkey_ot")])

## check that the tips (read from the tree) are a subset of the MLST in the dat_mge:
## NB: add this check in the function that designs the tree.
dat_mge <- dat_mge0

# the ones on the left are from data, the ones on the right are from the tree
# but how does the tree call the first node, in what order?
match_ind <- match(dat_mge$MLST,thetree$tip.label)
# not every one in the data is matched to the tree (the tree
# was built for leaves with more than 4 obs; so some ST with fewer observations
# are not included in the tree estimation). Now let's see how many
# isolates are there in these unmatched STs

sum(is.na(dat_mge0$MLST))      # 13 isolates have no sequence type info.
dat_mge$MLST[is.na(match_ind)] # remove these isolates from main analysis because they are not mapped to the tree.
ct_ST <- table(dat_mge$MLST[!is.na(match_ind)])
hist(ct_ST)
sort(ct_ST) # 133 MLSTs that appear in the tree, a minimum of 4 and a maximum of 223 isolates per leaf.

ct_ST_human_or_not0 <- table(dat_mge[!is.na(match_ind),c("MLST","human")])
##        human
## MLST   0  1
## ST10   64 31
## ST101  96  3
## ST1011  4  0
## ST106   4  0
## ST1079  9  0
## ST1141  7  0


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


###############################################################################
#  Setting up simulation and recovery of truths to illustrate 'lotR'
###############################################################################

# V(thetree_igraph)$levels       <- rep(2,length(V(thetree_igraph)))
# V(thetree_igraph)$levels[c(1)] <- 1
V(thetree_igraph)$levels <- rep(1,length(V(thetree_igraph)))
V(thetree_igraph)$levels[match(names(igraph::V(thetree_igraph)[igraph::degree(thetree_igraph, mode = "out") == 0]),
                               names(igraph::V(thetree_igraph)))] <- 2
# V(thetree_igraph)$levels[1] <- 3 # separate rootnode.
tau   <- c(0.5,0.3,0.15,0.05)

theta <- rbind(rep(rep(c(0.9, 0.9), each = 1),9),
               rep(rep(c(0.01, 0.01), each = 1),9),
               rep(rep(c(0.5, 0.5), each = 1),9),
               rep(rep(c(0.1, 0.7), each = 1),9)
)
image(t(theta))

#
# simulated data
#
fold <- 1
prop_obs <- 0.0
simu <- lotR_blca(2663*fold, itemprob = theta, classprob = tau)
Y <- simu$x
Z <- simu$z
Z_obs  <- as.matrix(cbind(1:nrow(Y),Z))
Z_obs[-sort(sample(1:nrow(Y),floor(prop_obs*nrow(Y)),replace=FALSE)),2] <- NA

curr_leaves <- rep(dat_mge[!is.na(match_ind),"ct_MLST"],fold)

#
# real data:
#
# Y <- rbind(dat_mge[!is.na(match_ind),ind_EL]
# )
# curr_leaves <- c(dat_mge[!is.na(match_ind),"ct_MLST"]
# )

K     <- length(tau)
p     <- length(V(thetree_igraph))
J     <- ncol(Y)

dsgn0  <- design_tree(Y,curr_leaves,thetree_igraph,rootnode = "Node1",weighted_edge = FALSE,Z_obs)

nrestarts <- 1
doParallel::registerDoParallel(cores = nrestarts)
#log_dir <- file.path("/Users/zhenkewu/Dropbox/ZW/professional/=research_agenda/latent_variable/latent_class_model_with_phylogenetic_tree/tree_lcm_code",
#"restart_logs")
log_dir <- tempdir()
dir.create(log_dir)
set.seed(345083)
par(mfrow=c(3,3));plot(0,0)

# profvis::profvis({
mod0     <- lcm_tree(Y,curr_leaves,thetree_igraph,
                     rootnode      = "Node1", # <-- may be redundant?
                     weighted_edge = TRUE,
                     Z_obs         = Z_obs,
                     parallel      = TRUE,
                     hyper_fixed   = list(K=K,a=c(1,1),b=c(1,1),
                                          tau_update_levels = c(1,2),
                                          s_u_zeroset = NULL,s_u_oneset = NULL),
                     #s_u_zeroset = (1:265)[-c(1,2,3,4,76,120)],s_u_oneset = c(1,2,3,4,76,120)),
                     #s_u_zeroset = (1:265)[-c(1)],s_u_oneset = c(1)),
                     # hyperparams_init = list(tau_1=rep(900/4,2),
                     #                         tau_2=rep(900/4,2)),
                     # hyperparams_init = list(tau_1=matrix(900/4,nrow=2,ncol=K-1),
                     #                         tau_2=array(900/4,c(2,J,K))),
                     vi_params_init = list(prob=rep(0.05,p)),
                     random_init = !TRUE,
                     nrestarts     = nrestarts,
                     quiet         = !TRUE,
                     plot_fig      = TRUE,
                     shared_tau    = FALSE,
                     print_freq    = 10,update_hyper_freq = 50, max_iter = 500,
                     tol           = 1e-8,
                     tol_hyper     = 1e-4,
                     allow_continue = FALSE,
                     log_restarts =!TRUE,log_dir = log_dir)
# })

# unlink(log_dir, recursive = T)
# closeAllConnections()

# summarize posterior results:
plot(mod0)

#image(mod0$prob_est$pi_collapsed)
image(mod0$prob_est$pi_collapsed)
which(mod0$mod$vi_params$prob>0.5)

#colnames(dsgn0$Y) <- c("EL2","EL3", "EL40", "EL35", "EL44", "EL45", "EL46", "EL12", "EL18","EL19", "EL36", "EL37", "EL38", "EL39" ,"EL41", "EL42", "EL43", "EL50")
# blca:
# VB:

fullLCM_vb <- BayesLCA::blca(dsgn0$Y,K,method="vb",verbose=FALSE)
fullLCM_vb
image(t(fullLCM_vb$itemprob),main="BayesLCA package: vb")

# EM:
fullLCM_em <- BayesLCA::blca(dsgn0$Y,K,method="em",verbose=FALSE)
fullLCM_em
image(t(fullLCM_em$itemprob),main="BayesLCA package: em")


barplot(rbind(c(tau,rep(0,K-length(tau))),
              sort(mod0$prob_est$pi_collapsed,decreasing = TRUE),
              fullLCM_vb$classprob,#,
              fullLCM_em$classprob
),
main="class probabilities", beside=TRUE)#,legend.text = c("truth","proposed","vb","em"))


## poLCA:
# form = as.formula("cbind(EL2,EL3, EL40, EL35, EL44, EL45, EL46, EL12, EL18,
#                   EL19, EL36, EL37, EL38, EL39 ,EL41, EL42, EL43, EL50)~1")
# for (j in 1:ncol(dsgn0$Y)){dsgn0$Y[,j] <- as.integer(as.factor(dsgn0$Y[,j]))}
# res2 <- poLCA::poLCA(form, as.data.frame(dsgn0$Y),nclass=K)

proposed <- mod0$prob_est$theta_collapsed

png("inst/example_figure/comparison_with_std.png",width=10,height=9,units = "in",res=300)
par(mfcol=c(3,3))

image(t(theta),main="truth")

# differences
image((proposed-t(theta)),main="diff")
hist((proposed-t(theta)),xlim=c(-1,1))

image(proposed,main="proposed")
image(t(fullLCM_vb$itemprob),main="BayesLCA package: vb")
image(t(fullLCM_em$itemprob),main="BayesLCA package: em")

# pi estimate:
barplot(rbind(c(tau,rep(0,K-length(tau))),
              # mod0$prob_est$pi_collapsed,
              sort(mod0$prob_est$pi_collapsed,decreasing = TRUE),
              fullLCM_vb$classprob,
              fullLCM_em$classprob),
        main="class probabilities", beside=TRUE,
        legend.text = c("truth","proposed","vb","em"))
dev.off()
