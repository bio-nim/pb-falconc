#!/usr/bin/env Rscript

require("ggplot2")
require("optparse")
option_list = list(make_option(c("-f", "--file"), type="character", default=NULL, help="raptor m4 pile", metavar="character"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

outfn<-paste(opt$file, ".rap-plot.pdf", sep="")
print(outfn)

dat<-read.table(opt$file,  colClasses=c("character", "character", "integer", "numeric", "integer", "integer", "integer", "integer",  "integer", "integer", "integer", "integer", "character", "character", "character", "character", "character" ))
print(head(dat))
dat<-dat[order(dat$V7),]

dat$o <- 1:length(dat$V1)
loc<-dat[dat$V13 == 'u',]
loc3<-loc[loc$V7 != loc$V8,]
loc2<-loc[loc$V6 != 0,]

pname<-names(table(dat$V1))

tp<-ggplot(dat, aes(x=V6, xend=V7, y=o, yend=o, colour=V13))+geom_segment(size=1.5)+geom_point(data=loc3, aes(x=V7, y=o), colour="black",size=2)+geom_point(data=loc2, aes(x=V6, y=o), colour="black", size=2)+labs(x="A-read pos", y="end sorted B-reads")+theme_classic(16)+theme(legend.position="top")+scale_color_brewer(name="overlap type", palette="Dark2")+ggtitle(paste ("pread:", pname))

ggsave(filename=outfn, plot=tp, width=20)
