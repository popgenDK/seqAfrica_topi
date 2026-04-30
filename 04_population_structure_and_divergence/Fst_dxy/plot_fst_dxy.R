library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)


tab <- read.table(paste0("topi.allsites.2dsfs.tsv"), header = TRUE)
country <- read.table("/projects/seqafrica/people/mdn487/Topi/topiBased/Topi_topi_variable_nomultiallelics_noindels_polymorphic_goodSites.fam")[, 1]

#mirror matrix around diagonal for plotting
tab_mirror <- tab
tab_mirror$A <- tab$B
tab_mirror$B <- tab$A
tab <- rbind(tab, tab_mirror)

# To not let very negative values mess up scale
tab$fst[tab$fst < -0.03] <- NA

# plotorder <- readLines(paste0(spec, ".plotorder"))

#print(unique(tab$A)[!(unique(tab$A) %in% plotorder)]) # test to see if any missing in plotorder file

# tab$A <- factor(tab$A, levels = plotorder)
# tab$B <- factor(tab$B, levels = plotorder)


p1 <- ggplot(tab, aes(A, B, fill = fst)) +
geom_tile(color = "#f8f8f8",
        lwd = 0.1,
        linetype = 1) + 
scale_fill_gradientn(colors = rev(viridis::viridis(9))) +
# scale_fill_gradient2(low = "blue", mid = "green", high = "yellow", midpoint = 0.325) +
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
ggtitle(paste("Pairwise Fst for topi"))

p2 <- ggplot(tab, aes(A, B, fill = dxy)) +
geom_tile(color = "#f8f8f8",
        lwd = 0.1,
        linetype = 1) + 
scale_fill_gradientn(colors = rev(viridis::viridis(9))) +
# scale_fill_gradient2(low = "blue", mid = "green", high = "yellow", midpoint = 0.325) +
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
ggtitle(paste("Pairwise dXY for topi"))


last_in_loc <- cumsum(rle(sort(country))$lengths)

p1 <- p1 + geom_vline(xintercept = last_in_loc + 0.5, color = "white", size = 1.5) +
geom_hline(yintercept = last_in_loc + 0.5, color = "white", size = 1.5)
p2 <- p2 + geom_vline(xintercept = last_in_loc + 0.5, color = "white", size = 1.5) +
geom_hline(yintercept = last_in_loc + 0.5, color = "white", size = 1.5)

png(paste0("topi.allsites.fst_dxy.png"), width=2640, height=1080)
grid.arrange(p1, p2, ncol = 2)
dev.off()

