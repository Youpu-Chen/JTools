#-------------------------------------------------------------------------------
# this script is used to run Tetmer on the k-mer spectrum
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# load packages
#-------------------------------------------------------------------------------
rm(list = ls())
library(Tetmer)
library(argparser)


# Create a parser
p <- arg_parser("Tetmer plot")

# Add command line arguments
p <- add_argument(p, "--input", help="the input name of kmer spectrum", type="character")
p <- add_argument(p, "--output", help="the output name of Tetmer plot", type = "character")
p <- add_argument(p, "--kmer", help="number of decimal places", type = "numeric")

# Parse the command line arguments
argv <- parse_args(p)


Tetmer_plot_title <- strsplit(argv$input, split = "[.]")[[1]][1]
mySpec = read.spectrum(argv$input, Tetmer_plot_title, argv$kmer)
# tetmer(mySpec)
pdf(file=argv$output, width = 8, height = 6)
plot(mySpec, xlim=c(0, 150), ylim=c(0, 10000000))
dev.off()