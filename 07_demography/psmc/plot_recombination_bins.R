########### PLot recombination in each bin from Genis
psmc.result.all <-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
  X<-scan(file=file,what="",sep="\n",quiet=TRUE)
  
  START<-grep("^RD",X)
  END<-grep("^//",X)
  
  X<-X[START[i.iteration+1]:END[i.iteration+1]]
  
  TR<-grep("^TR",X,value=TRUE)
  RS<-grep("^RS",X,value=TRUE)
  MM <- grep("^MM", X, value=T)
  
  write(TR,"temp.psmc.result")
  theta0 <- as.numeric(read.table("temp.psmc.result")[1,2])
  rho0 <- as.numeric(read.table("temp.psmc.result")[1,2])
  Cpi <- as.numeric(strsplit(strsplit(MM,split=",")[[1]][1], ":")[[1]][2])
  Csigma <- 1 / (rho0 * Cpi)
  N0<-theta0/4/mu/s # this is the same as theta0/(4*mu)/s, don't panic next time
  
  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  
  Generation<-as.numeric(2*N0*a[,3]) # a[,3] is t_k
  Ne<-as.numeric(N0*a[,4]) #a[,4] is lambda_k
  Nrecomb <- as.numeric(a[,5])
  
  n.points<-length(Ne)
  YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  NRecombinations <- c(as.numeric(rbind(Nrecomb[-n.points],Nrecomb[-n.points])),
                       Nrecomb[n.points])
  
  file.remove("temp.psmc.result")
  data.frame(YearsAgo, Ne, NRecombinations)
}

files <-list.files("/path/psmc/results/PSMC_OUT",".psmc",full.names=T)

mu <- 1e-8
g <- 7.5
psmcFiles <- files
res <- lapply(psmcFiles, psmc.result,mu=mu,g=g)
names(res) <- gsub(".psmc", "",basename(psmcFiles))

inds <- names(res)
pop = as.data.frame(read.table("/path/topi/topi68inds.poplist", col.names=c("Pop","Ind","Color")))
res2 <- lapply(res, function(x) x[-((nrow(x) - 8):nrow(x)),])

# ---- Palette and group mapping ----
ind_to_pop <- setNames(pop$Pop, pop$Ind)
final_palette <- setNames(pop$Color, pop$Pop)

res_recomb <- lapply(files, psmc.result.all, mu = mu, g = g)
names(res_recomb) <- gsub(".psmc", "", basename(files))
res_recomb <- lapply(res_recomb, function(x) x[-((nrow(x) - 8):nrow(x)),])

### Loop over populations ### 
for (pop_name in unique(pop$Pop)) {
  # Get individuals in this population
  inds_in_pop <- pop$Ind[pop$Pop == pop_name]
  inds_in_pop <- intersect(inds_in_pop, names(res_recomb))
  
  if (length(inds_in_pop) == 0) {
    message("Skipping population ", pop_name, " — no matching individuals.")
    next
  }
  # Extract data for plotting
  pop_curves <- res_recomb[inds_in_pop]
  
  
  # Get Y-axis range across this population
  y_vals <- unlist(lapply(pop_curves, function(x) x$NRecombinations))
  y_vals <- y_vals[is.finite(y_vals) & y_vals > 0]
  
  if (length(y_vals) == 0) {
    warning("Skipping population ", pop_name, ": no valid recombination values.")
    next
  }

  # Determine Y-axis log-scaled limits
  y_min <- 10^(floor(log10(min(y_vals))))
  y_max <- 10^(ceiling(log10(max(y_vals))))

  png(file = paste0("recomb_", pop_name, ".png"), width = 3200, height = 2400, res = 300)
  par(mar = c(5, 5, 4, 1) + 0.1)

  plot(x = pop_curves[[1]]$YearsAgo, y = pop_curves[[1]]$NRecombinations, type = 'l', log = 'xy', col = final_palette[pop_name],
     lwd = 4, xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', ylim = c(y_min, y_max), cex.lab = 1.5, cex.axis = 1.5, main = paste("Population:", pop_name))

  # Add custom y-axis ticks (log scale)
  y_ticks <- 10^(floor(log10(y_min)):ceiling(log10(y_max)))
  axis(2, at = y_ticks, labels = y_ticks, las = 1, cex.axis = 1.2)

  # Custom x-axis (optional, if you want ticks)
  x_ticks <- 10^(seq(3, 6, by = 0.5))
  axis(1, at = x_ticks, labels = formatC(x_ticks, format = "e"), cex.axis = 1.2)

  # Plot all individuals
  for (ind in inds_in_pop[-1]) {
    lines(x = pop_curves[[ind]]$YearsAgo,
          y = pop_curves[[ind]]$NRecombinations,
          col = final_palette[pop_name],
          lwd = 2)
  }

  title(xlab = sprintf("Years ago (mu=%.2e, g=%.1f)", mu, g),
        cex.lab = 2, line = 3.5)
  title(ylab = "Number of recombinations in bin",
        cex.lab = 2, line = 3.5)
  
  dev.off()
}
