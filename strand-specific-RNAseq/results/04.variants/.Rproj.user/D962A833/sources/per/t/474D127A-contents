library(ggplot2)

dat <- read.table('ST.txt', sep = '\t')
names(dat) <- c('Type', 'Count')
dat$group <- factor(rep(c('a', 'c', 'g', 't'), each = 3))

p1 <- ggplot(data = dat, aes(x = Type, y = Count)) + geom_bar(stat = 'identity', aes(fill = group)) + 
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.position="none")
p1

# ggsave('ST-count.pdf', height = 6, width = 10)
