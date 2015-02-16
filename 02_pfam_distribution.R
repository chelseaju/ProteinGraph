
##	Plot histogram for pfam

library(ggplot2)
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 2){
        stop(paste("Invalid Arguments\n",
        "Usage: R --no-save --slave < 02_pfam_distribution.R --args inputfile outputdir\n",
        "\t iputfile = input filename \n",
        "\t outputdir = output directory\n",
        sep=""));
}

file <- options[1];
dir <- options[2];

data <- read.table(file, header = T);
colnames(data) <- c("Name", "ProteinCount");

png(paste(dir, "/", "pfam_histogram_small.png", sep=""), height = 800, width = 800)
small_data <- data[data$ProteinCount < 50,];
qplot(ProteinCount, data=small_data, geom="histogram", binwidth = 2);
dev.off();


png(paste(dir, "/", "pfam_histogram_mid.png", sep=""), height = 800, width = 800)
mid_data <- data[data$ProteinCount > 25 & data$ProteinCount < 300,];
qplot(ProteinCount, data=mid_data, geom="histogram", binwidth = 2);
dev.off();


png(paste(dir, "/", "pfam_histogram_big.png", sep=""), height = 800, width = 800)
big_data <- data[data$ProteinCount >= 300,];
qplot(ProteinCount, data=big_data, geom="histogram", binwidth = 2);
dev.off();



















