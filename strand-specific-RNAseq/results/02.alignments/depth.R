rm(list=ls())

# Load datasets
fwd.df <- read.delim('fwd.depth.txt', sep = '\t', header = FALSE)
names(fwd.df) <- c('chr', 'start', 'end', 'read.counts', 'base.number', 'window.size', 'average.coverage')
head(fwd.df)

rev.df <- read.delim('rev.depth.txt', sep = '\t', header = FALSE)
names(rev.df) <- c('chr', 'start', 'end', 'read.counts', 'base.number', 'window.size', 'average.coverage')
head(rev.df)


# Load packages
library(ggplot2)
library(ggpubr)


# plot
# 4    chrIV   74
# 12  chrXII   58
# 8   chrVII   52

if(F){
  ### chrIV 
  fwd.chrI <- fwd.df[which(fwd.df$chr == "chrIV"), ]
  rev.chrI <- rev.df[which(rev.df$chr == "chrIV"), ]
  p1 <- ggplot(data = fwd.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="red", alpha=0.5, position="identity") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), 
                                                              breaks = c(0, 250, 500, 750), 
                                                              limits = c(0, 750)) + 
    theme_classic() + 
    xlab("Chromosome IV (Bp)") +
    # ylab("Read Counts/Per 1kb") + 
    # xlab("") + 
    ylab("") +
    # theme(
    #   axis.title.x = element_text(vjust=2)
    # ) + 
    annotate("text", x = 1350000, y = 600, label = "Positive Strand", size=5)
  p1
  
  
  p2 <- ggplot(data = rev.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="blue", alpha=0.5, position="identity") + 
    theme_classic() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), 
                                                              breaks = c(0, 250, 500, 750)) +  
    scale_y_reverse(breaks = c(0, 250, 500, 750),
                    limits = c(750, 0)) + 
    geom_hline(yintercept = 0) + 
    xlab("") + 
    ylab("") +
    annotate("text", x = 1350000, y = 600, label = "Negative Strand", size=5)
  p2
  
  
  figure <- ggarrange(p1, p2,
                      ncol = 1, nrow = 2)
  figure
  
  
  anno.fig <- annotate_figure(
    figure,
    left = text_grob("Read Depth (Per 1kb)",
                     color = "black", rot = 90)
  )
  anno.fig
  
  #save 
  ggsave('chrIV-readdepth.pdf', anno.fig, height = 6, width = 10)
}


if(F){
  ### chrXII 
  fwd.chrI <- fwd.df[which(fwd.df$chr == "chrXII"), ]
  rev.chrI <- rev.df[which(rev.df$chr == "chrXII"), ]
  p1 <- ggplot(data = fwd.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="red", alpha=0.5, position="identity") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), 
                                                              breaks = c(0, 250, 500, 750), 
                                                              limits = c(0, 750)) + 
    theme_classic() + 
    xlab("Chromosome XII (Bp)") +
    # ylab("Read Counts/Per 1kb") + 
    # xlab("") + 
    ylab("") +
    # theme(
    #   axis.title.x = element_text(vjust=2)
    # ) + 
    annotate("text", x = 860000, y = 600, label = "Positive Strand", size=5)
  p1
  
  
  p2 <- ggplot(data = rev.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="blue", alpha=0.5, position="identity") + 
    theme_classic() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), 
                                                              breaks = c(0, 250, 500, 750)) +  
    scale_y_reverse(breaks = c(0, 250, 500, 750),
                    limits = c(750, 0)) + 
    geom_hline(yintercept = 0) + 
    xlab("") + 
    ylab("") +
    annotate("text", x = 860000, y = 600, label = "Negative Strand", size=5)
  p2
  
  
  figure <- ggarrange(p1, p2,
                      ncol = 1, nrow = 2)
  figure
  
  
  anno.fig <- annotate_figure(
    figure,
    left = text_grob("Read Depth (Per 1kb)",
                     color = "black", rot = 90)
  )
  anno.fig
  
  #save 
  ggsave('chrXII-readdepth.pdf', anno.fig, height = 6, width = 10)
}


if(F){
  ### chrVII 
  fwd.chrI <- fwd.df[which(fwd.df$chr == "chrVII"), ]
  rev.chrI <- rev.df[which(rev.df$chr == "chrVII"), ]
  p1 <- ggplot(data = fwd.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="red", alpha=0.5, position="identity") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), 
                                                              breaks = c(0, 250, 500, 750), 
                                                              limits = c(0, 750)) + 
    theme_classic() + 
    xlab("Chromosome VII (Bp)") +
    # ylab("Read Counts/Per 1kb") + 
    # xlab("") + 
    ylab("") +
    # theme(
    #   axis.title.x = element_text(vjust=2)
    # ) + 
    annotate("text", x = 860000, y = 600, label = "Positive Strand", size=5)
  p1
  
  
  p2 <- ggplot(data = rev.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="blue", alpha=0.5, position="identity") + 
    theme_classic() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), 
                                                              breaks = c(0, 250, 500, 750)) +  
    scale_y_reverse(breaks = c(0, 250, 500, 750),
                    limits = c(750, 0)) + 
    geom_hline(yintercept = 0) + 
    xlab("") + 
    ylab("") +
    annotate("text", x = 860000, y = 700, label = "Negative Strand", size=5)
  p2
  
  
  figure <- ggarrange(p1, p2,
                      ncol = 1, nrow = 2)
  figure
  
  
  anno.fig <- annotate_figure(
    figure,
    left = text_grob("Read Depth (Per 1kb)",
                     color = "black", rot = 90)
  )
  anno.fig
  
  #save 
  ggsave('chrVII-readdepth.pdf', anno.fig, height = 6, width = 10)
}


if(F){
  ### chrII; chrM 
  fwd.chrI <- fwd.df[which(fwd.df$chr == "chrM"), ]
  rev.chrI <- rev.df[which(rev.df$chr == "chrM"), ]
  p1 <- ggplot(data = fwd.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="red", alpha=0.5, position="identity") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), 
                                                              breaks = c(0, 250, 500, 750), 
                                                              limits = c(0, 750)) + 
    theme_classic() + 
    xlab("Mitochondria (Bp)") +
    # ylab("Read Counts/Per 1kb") + 
    # xlab("") + 
    ylab("") +
    # theme(
    #   axis.title.x = element_text(vjust=2)
    # ) + 
    annotate("text", x = 70000, y = 400, label = "Positive Strand", size=5)
  p1
  
  
  p2 <- ggplot(data = rev.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="blue", alpha=0.5, position="identity") + 
    theme_classic() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), 
                                                              breaks = c(0, 250, 500, 750)) +  
    scale_y_reverse(breaks = c(0, 250, 500, 750),
                    limits = c(750, 0)) + 
    geom_hline(yintercept = 0) + 
    xlab("") + 
    ylab("") +
    annotate("text", x = 70000, y = 400, label = "Negative Strand", size=5)
  p2
  
  
  figure <- ggarrange(p1, p2,
                      ncol = 1, nrow = 2)
  figure
  
  
  anno.fig <- annotate_figure(
    figure,
    left = text_grob("Read Depth (Per 1kb)",
                     color = "black", rot = 90)
  )
  anno.fig
  
  #save 
  ggsave('Mitochondria-readdepth.pdf', anno.fig, height = 6, width = 10)
}


if(F){
  ### chromosome I
  fwd.chrI <- fwd.df[which(fwd.df$chr == "chrI"), ]
  rev.chrI <- rev.df[which(rev.df$chr == "chrI"), ]
  p1 <- ggplot(data = fwd.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="red", alpha=0.5, position="identity") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
    theme_classic() + 
    xlab("Chromosome I") +
    # ylab("Read Counts/Per 1kb") + 
    # xlab("") + 
    ylab("") +
    # theme(
    #   axis.title.x = element_text(vjust=2)
    # ) + 
    annotate("text", x = 190000, y = 9500, label = "Positive Strand", size=5)
  p1
  
  
  p2 <- ggplot(data = rev.chrI, aes(x=(start+end)/2, y = read.counts)) + 
    # geom_line() + 
    geom_area(colour = "white", fill="blue", alpha=0.5, position="identity") + 
    theme_classic() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
    scale_y_reverse() + 
    geom_hline(yintercept = 0) + 
    xlab("") + 
    ylab("") +
    annotate("text", x = 190000, y = 34000, label = "Negative Strand", size=5)
  p2
  
  
  figure <- ggarrange(p1, p2,
                      ncol = 1, nrow = 2)
  figure
  
  
  anno.fig <- annotate_figure(
    figure,
    left = text_grob("Read Depth (Per 1kb)",
                     color = "black", rot = 90)
  )
  anno.fig
  
  #save 
  ggsave('chrI-readdepth.pdf', anno.fig, height = 6, width = 10)
}
