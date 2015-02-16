
##	Plot histogram for pfam

library(ggplot2)
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 2){
        stop(paste("Invalid Arguments\n",
        "Usage: R --no-save --slave < 02_pdb_distribution.R --args inputfile outputdir\n",
        "\t iputfile = input filename \n",
        "\t outputdir = output directory\n",
        sep=""));
}

file <- options[1];
dir <- options[2];

data <- read.table(file, header = T);
colnames(data) <- c("Name", "Pfam_Count");

png(paste(dir, "/", "pdb_histogram.png", sep=""), height = 800, width = 800)
qplot(Pfam_Count, data=data, geom="histogram", binwidth = 2);
dev.off();

