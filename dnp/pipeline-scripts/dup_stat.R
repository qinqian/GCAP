args <- commandArgs(trailingOnly = TRUE)

read_no_hr <- function(x) read.table(x, as.is=T, row.names=NULL, header =F)
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

output <- args[1]
id <- args[2]
type <- args[3]
print(args)

dup_bed <- list()

for (i in args[4:length(args)]) {
  dup_bed[[i]] <- read_no_hr(i)
  dup_bed[[i]] <- table(dup_bed[[i]]$V4)
  dup_bed[[i]] <- as.data.frame(dup_bed[[i]])
}

if (is.installed("ggplot2") && type == "barplot") {
  suppressMessages(require(ggplot2))
  for (i in seq(along.with = dup_bed)) {
    pdf(output, width=16, height=8)
      a <- ggplot(dup_bed[[i]], aes(x=Var1, y=Freq)) + geom_bar(stat = "identity", position = "stack",width=0.5)
      a <- a + labs(x = "Tags Number in the same location", y = expression(log10(Number)), fill = NULL)
      a <- a + theme(
           axis.text.x = element_text(size=7,  colour="grey20", angle=60, vjust=2, hjust=2),
           axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"))
      print(a)
  }
  } else if (type == "scatter") {
  pdf(output, width=16, height=8)
  ## barplot(test.table, col="blue", log="y", ylim=c(1, max(test.table)))
  for (i in seq(along.with = dup_bed)) {
    par(cex=1.3, mar=c(2,5,2,2))
    if (i == 1) {
        plot(as.numeric(dup_bed[[i]]$Var1), dup_bed[[i]]$Freq, log="y", col=i, main = "Distribution of tags", type = "h", pch=18, lwd=1.2,
             xlab = "Duplicates number", ylab = "Location Frequency") } else {
        lines(as.numeric(dup_bed[[i]]$Var1), dup_bed[[i]]$Freq, col = i)
        }
      }
  }

dev.off()