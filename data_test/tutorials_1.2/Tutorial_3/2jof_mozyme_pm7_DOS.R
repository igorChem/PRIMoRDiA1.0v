#!/usr/bin/env Rscript
require(ggplot2) 
dos = read.table( '2jof_mozyme_pm7.DOS',header=T)
attach(dos) 
p <-ggplot(dos, aes( x=Energy) )+
	geom_density(fill='blue',bw=1) 
png('2jof_mozyme_pm7.DOS.png',width = 5, height = 3.5, units ='in', res = 400)
p
dev.off()