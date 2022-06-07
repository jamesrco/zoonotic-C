library(plotrix)
library(RColorBrewer)

# read in biomass, export data

rawCdat <- read.csv2("/Users/jamesrco/Dropbox/Research Materials & Journal Articles/Ocean NCS & CDR/EDF Bezos Earth Fund project/Reports/BEF State of science report figures/Open_ocean_NCS_C_data.csv", header = T, sep = ",")

# convert data types; append a field for a numerical reference # in the final figures

rawCdat[,6:14] <- lapply(rawCdat[,6:14],as.numeric)
rawCdat$refNum <- as.numeric(as.factor(rawCdat$Reference))

# subsetting

Biomass.today <- rawCdat[rawCdat$Biomass.or.export=="Biomass",]

# factor variable by subcategory, for the range plot

yrange.Biomass <- as.factor(Biomass.today$Subcategory)

# make a range plot with the biomass data, by subcategory; using plotrix since we need an axis break

gap.plot(Biomass.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
         as.numeric(yrange.Biomass),
         pch=0, col = 0,
         gap=c(3.25,23.75), gap.axis="x",
         ylim=c(0.5,4.5),
         ylab = "Category",
         xlab = "Biomass (Pg C)",
     xlim=c(min(Biomass.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.),
            max(c(Biomass.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,Biomass.today$Current.hi.bound..Pg.C.or.Pg.C.yr.1.), na.rm = T)),
     xtics = c(seq(0,25,1)),breakcol="black")

# create range rectangles using min/max bounds from uncertainties or values

rectBounds <- data.frame(c(min(unlist(c(Biomass.today[Biomass.today$Subcategory=="All fish",9:10])), na.rm = T),
                min(unlist(c(Biomass.today[Biomass.today$Subcategory=="Mesopelagic species",9:10])), na.rm = T),
                min(unlist(c(Biomass.today[Biomass.today$Subcategory=="Whales",9:10])), na.rm = T),
                min(unlist(c(Biomass.today[Biomass.today$Subcategory=="z - All ocean animals (upper limit)",9:10])), na.rm = T)),
                c(sort(unique(as.numeric(yrange.Biomass)))-.25),
                c(max(unlist(c(Biomass.today[Biomass.today$Subcategory=="All fish",c(9,11)])), na.rm = T),
                  max(unlist(c(Biomass.today[Biomass.today$Subcategory=="Mesopelagic species",c(9,11)])), na.rm = T),
                  max(unlist(c(Biomass.today[Biomass.today$Subcategory=="Whales",c(9,11)])), na.rm = T),
                  max(unlist(c(Biomass.today[Biomass.today$Subcategory=="z - All ocean animals (upper limit)",c(9,11)])), na.rm = T)),
                c(sort(unique(as.numeric(yrange.Biomass)))+.25))

gap.plot(rect(rectBounds[,1],rectBounds[,2],rectBounds[,3],rectBounds[,4],col="lightgrey",lty = 0), add = T, gap=c(3.25,23.75), gap.axis="x",ylim=c(0.5,4.5),
         pch = 0, col = 0,
         xlim=c(min(Biomass.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.),
                max(c(Biomass.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,Biomass.today$Current.hi.bound..Pg.C.or.Pg.C.yr.1.), na.rm = T)))

axis.break(axis = 1, breakpos = 3.25) # add the axis break symbol

# append some vertical lines that show where the reference values are

segments(Biomass.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
         as.numeric(yrange.Biomass)-.25,
         y1 = as.numeric(yrange.Biomass)+.25)

# text with our reference #s

text(Biomass.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
     as.numeric(yrange.Biomass),
     Biomass.today$refNum, cex = 0.6)

# plot the flux data

# subset first

Export.today <- rawCdat[rawCdat$Biomass.or.export=="Export",]

# factor variable by subcategory, for the range plot

xrange.Export <- as.factor(Export.today$Subcategory)

# make a plot, log scale first, then a smaller on regular scale for use as an inset

# current values 

plot(as.numeric(xrange.Export), Export.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1., type = "p", col = 0, xlim=c(0.5,5.5), ylim=c(250,0.000001), log="y", xlab = "Flux", ylab = "Export from surface ocean (Pg C yr-1)", yaxt = "n", yaxs="i")
axty <- c(axTicks(2),1e-06)
myTickNumbers <- as.numeric(sapply(log10(axty), format, scientifc = FALSE))
myTickLabels <- sapply(myTickNumbers,function(i)
  as.expression(bquote(10^ .(i))))
myTickLabels[c(1:3)] <- c(100,10,1)

axis(2, at = axty, labels = myTickLabels)

# range rectangles

# current estimates

rectBounds.Export <- data.frame(c(sort(unique(as.numeric(xrange.Export)))-.35),
                       rep(0, length(levels(xrange.Export))),
                       c(sort(unique(as.numeric(xrange.Export)))-.05),
                       c(max(unlist(c(Export.today[Export.today$Subcategory=="aa - Carbon cycle fluxes",c(9,11)])), na.rm = T),
                         max(unlist(c(Export.today[Export.today$Subcategory=="All fish",c(9,11)])), na.rm = T),
                         max(unlist(c(Export.today[Export.today$Subcategory=="Epipelagic species",c(9,11)])), na.rm = T),
                         max(unlist(c(Export.today[Export.today$Subcategory=="Mesopelagic species",c(9,11)])), na.rm = T),
                         max(unlist(c(Export.today[Export.today$Subcategory=="Whales",c(9,11)])), na.rm = T)))
                       
rect(rectBounds.Export[,1],rep(0.0000001, length(levels(xrange.Export))),rectBounds.Export[,3],rectBounds.Export[,4],col="lightgrey",lty = 0)

# pre-whaling estimates for whales and for all fish (from Bianchi et al. 2021)

rectBounds.Export.bestcase <- data.frame(c(sort(unique(as.numeric(xrange.Export)))+.05),
                                rep(0, length(levels(xrange.Export))),
                                c(sort(unique(as.numeric(xrange.Export)))+.35),
                                c(max(unlist(c(Export.today[Export.today$Subcategory=="aa - Carbon cycle fluxes",c(6,8)])), na.rm = T),
                                  max(unlist(c(Export.today[Export.today$Subcategory=="All fish",c(6,8)])), na.rm = T),
                                  max(unlist(c(Export.today[Export.today$Subcategory=="Epipelagic species",c(12,14)])), na.rm = T),
                                  max(unlist(c(Export.today[Export.today$Subcategory=="Mesopelagic species",c(12,14)])), na.rm = T),
                                  max(unlist(c(Export.today[Export.today$Subcategory=="Whales",c(6,8)])), na.rm = T)))

rect(rectBounds.Export.bestcase[,1],rep(0.0000001, length(levels(xrange.Export))),rectBounds.Export.bestcase[,3],rectBounds.Export.bestcase[,4],col="lightgrey",lty = 0)

# for whales, separate rectangles for the with and without indirect fertilization component 

rect(rectBounds.Export[5,1],
     max(unlist(c(Export.today[Export.today$Subcategory=="Whales" & Export.today$Export.type=="Whale fall",c(9,11)])), na.rm = T),
     rectBounds.Export[5,3],
     max(unlist(c(Export.today[Export.today$Subcategory=="Whales" & Export.today$Export.type!="Whale fall",c(9,11)])), na.rm = T),
     col="darkgrey",lty = 0)

rect(rectBounds.Export.bestcase[5,1],
     min(unlist(c(Export.today[Export.today$Subcategory=="Whales" & Export.today$Export.type=="Whale fall",c(6,8)])), na.rm = T),
     rectBounds.Export.bestcase[5,3],
     max(unlist(c(Export.today[Export.today$Subcategory=="Whales" & Export.today$Export.type!="Whale fall",c(6,8)])), na.rm = T),
     col="darkgrey",lty = 0)

# superimpose lines

segments(as.numeric(xrange.Export)-.35,
         Export.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
         x1 = as.numeric(xrange.Export)-.05)

segments(as.numeric(xrange.Export)[c(1,4,7)]+.05,
         Export.today$Pre.whaling.pre.industrial.fishing..Pg.C.or.Pg.C.yr.1.[c(1,4,7)],
         x1 = as.numeric(xrange.Export)[c(1,4,7)]+.35)

# superimpose text refs

text(as.numeric(xrange.Export)-0.2,
     Export.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
     Export.today$refNum, cex = 0.6)

text(as.numeric(xrange.Export)[c(1,4,7)]+0.2,
     Export.today$Pre.whaling.pre.industrial.fishing..Pg.C.or.Pg.C.yr.1.[c(1,4,7)],
     Export.today$refNum[c(1,4,7)], cex = 0.6)

# separate text for the overall export range estimates

text(1-0.38,
     10^-2.5,
     paste(max(Export.today$refNum),", ",max(Export.today$refNum)+1, sep = ""), cex = 0.6, srt = 90)

# non-log plot

plot(as.numeric(xrange.Export), Export.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1., type = "p", col = 0, xlim=c(0.5,5.5), ylim=c(17,0),xlab = "Flux", ylab = "Export from surface ocean (Pg C yr-1)", yaxs="i")

# range rectangles

rectBounds.Export <- data.frame(c(sort(unique(as.numeric(xrange.Export)))-.35),
                                rep(0, length(levels(xrange.Export))),
                                c(sort(unique(as.numeric(xrange.Export)))-.05),
                                c(max(unlist(c(Export.today[Export.today$Subcategory=="aa - Carbon cycle fluxes",c(9,11)])), na.rm = T),
                                  max(unlist(c(Export.today[Export.today$Subcategory=="All fish",c(9,11)])), na.rm = T),
                                  max(unlist(c(Export.today[Export.today$Subcategory=="Epipelagic species",c(9,11)])), na.rm = T),
                                  max(unlist(c(Export.today[Export.today$Subcategory=="Mesopelagic species",c(9,11)])), na.rm = T),
                                  max(unlist(c(Export.today[Export.today$Subcategory=="Whales",c(9,11)])), na.rm = T)))

rect(rectBounds.Export[,1],rep(0, length(levels(xrange.Export))),rectBounds.Export[,3],rectBounds.Export[,4],col="lightgrey",lty = 0)

rectBounds.Export.bestcase <- data.frame(c(sort(unique(as.numeric(xrange.Export)))+.05),
                                         rep(0, length(levels(xrange.Export))),
                                         c(sort(unique(as.numeric(xrange.Export)))+.35),
                                         c(max(unlist(c(Export.today[Export.today$Subcategory=="aa - Carbon cycle fluxes",c(6,8)])), na.rm = T),
                                           max(unlist(c(Export.today[Export.today$Subcategory=="All fish",c(6,8)])), na.rm = T),
                                           max(unlist(c(Export.today[Export.today$Subcategory=="Epipelagic species",c(12,14)])), na.rm = T),
                                           max(unlist(c(Export.today[Export.today$Subcategory=="Mesopelagic species",c(12,14)])), na.rm = T),
                                           max(unlist(c(Export.today[Export.today$Subcategory=="Whales",c(6,8)])), na.rm = T)))

rect(rectBounds.Export.bestcase[,1],rep(0, length(levels(xrange.Export))),rectBounds.Export.bestcase[,3],rectBounds.Export.bestcase[,4],col="lightgrey",lty = 0)

# for whales, separate rectangles for the with and without indirect fertilization component 

rect(rectBounds.Export[5,1],
     max(unlist(c(Export.today[Export.today$Subcategory=="Whales" & Export.today$Export.type=="Whale fall",c(9,11)])), na.rm = T),
     rectBounds.Export[5,3],
     max(unlist(c(Export.today[Export.today$Subcategory=="Whales" & Export.today$Export.type!="Whale fall",c(9,11)])), na.rm = T),
     col="darkgrey",lty = 0)

rect(rectBounds.Export.bestcase[5,1],
     min(unlist(c(Export.today[Export.today$Subcategory=="Whales" & Export.today$Export.type=="Whale fall",c(6,8)])), na.rm = T),
     rectBounds.Export.bestcase[5,3],
     max(unlist(c(Export.today[Export.today$Subcategory=="Whales" & Export.today$Export.type!="Whale fall",c(6,8)])), na.rm = T),
     col="darkgrey",lty = 0)

# superimpose lines

segments(as.numeric(xrange.Export)-.35,
         Export.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
         x1 = as.numeric(xrange.Export)-.05)

segments(as.numeric(xrange.Export)[c(1,4,7)]+.05,
         Export.today$Pre.whaling.pre.industrial.fishing..Pg.C.or.Pg.C.yr.1.[c(1,4,7)],
         x1 = as.numeric(xrange.Export)[c(1,4,7)]+.35)

# superimpose text refs

text(as.numeric(xrange.Export)-0.2,
     Export.today$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
     Export.today$refNum, cex = 0.6)

text(as.numeric(xrange.Export)[c(1,4,7)]+0.2,
     Export.today$Pre.whaling.pre.industrial.fishing..Pg.C.or.Pg.C.yr.1.[c(1,4,7)],
     Export.today$refNum[c(1,4,7)], cex = 0.6)

# separate text for the overall export range estimates

text(1-0.38,
     0.5*max(Export.today[Export.today$Subcategory=="aa - Carbon cycle fluxes",c("Current.hi.bound..Pg.C.or.Pg.C.yr.1.")]),
     paste(max(Export.today$refNum),", ",max(Export.today$refNum)+1, sep = ""), cex = 0.6, srt = 90)

# plot of fishing emissions 

FishEmisstoAtm <- rawCdat[rawCdat$Biomass.or.export=="Flux from ocean to atmosphere",]

# log plot

plot(1, FishEmisstoAtm$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1., type = "p", col = 0, xlim=c(0.5,5.5), ylim=c(0.000001,250), log="y", xlab = "Flux", ylab = "Flux to atmosphere (Pg C yr-1)", yaxt = "n", yaxs="i")
axty <- c(1e-06,axTicks(2))
myTickNumbers <- as.numeric(sapply(log10(axty), format, scientifc = FALSE))
myTickLabels <- sapply(myTickNumbers,function(i)
  as.expression(bquote(10^ .(i))))
myTickLabels[c(7:9)] <- c(1,10,100)

axis(2, at = axty, labels = myTickLabels)

rectBounds.FishEmisstoAtm <- data.frame(c(1-.35),
                                rep(0, 1),
                                c(1-.05),
                                c(max(unlist(c(FishEmisstoAtm[c(9,11)])), na.rm = T)))

rect(rectBounds.FishEmisstoAtm[,1],rep(0.000001, 1),rectBounds.FishEmisstoAtm[,3],rectBounds.FishEmisstoAtm[,4],col="lightgrey",lty = 0)

# superimpose lines

segments(1-.35,
         FishEmisstoAtm$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
         x1 = 1-.05)

# superimpose text refs

text(1-0.2,
     FishEmisstoAtm$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
     FishEmisstoAtm$refNum, cex = 0.6)

# non-log plot

plot(1, FishEmisstoAtm$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1., type = "p", col = 0, xlim=c(0.5,5.5), ylim=c(0,17),xlab = "Flux", ylab = "Flux to atmosphere (Pg C yr-1)", yaxs="i")

rectBounds.FishEmisstoAtm <- data.frame(c(1-.35),
                                        rep(0, 1),
                                        c(1-.05),
                                        c(max(unlist(c(FishEmisstoAtm[c(9,11)])), na.rm = T)))

rect(rectBounds.FishEmisstoAtm[,1],rep(0,1),rectBounds.FishEmisstoAtm[,3],rectBounds.FishEmisstoAtm[,4],col="lightgrey",lty = 0)

# superimpose lines

segments(1-.35,
         FishEmisstoAtm$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
         x1 = 1-.05)

# superimpose text refs

text(1-0.2,
     FishEmisstoAtm$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
     FishEmisstoAtm$refNum, cex = 0.6)

# plot of trawling sediment remineralization

TrawlSedRemin <- rawCdat[rawCdat$Biomass.or.export=="Flux from sediment",]

# log plot

plot(rep(1,2), TrawlSedRemin$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1., type = "p", col = 0, xlim=c(0.5,5.5), ylim=c(0.000001,250), log="y", xlab = "Flux", ylab = "Seafloor sediment remineralized (Pg C yr-1)", yaxt = "n", yaxs="i")
axty <- c(1e-06,axTicks(2))
myTickNumbers <- as.numeric(sapply(log10(axty), format, scientifc = FALSE))
myTickLabels <- sapply(myTickNumbers,function(i)
  as.expression(bquote(10^ .(i))))
myTickLabels[c(7:9)] <- c(1,10,100)

axis(2, at = axty, labels = myTickLabels)

rectBounds.TrawlSedRemin<- data.frame(c(1-.35),
                                        rep(0, 1),
                                        c(1-.05),
                                        c(max(unlist(c(TrawlSedRemin[c(9,11)])), na.rm = T)))

rect(rectBounds.TrawlSedRemin[,1],rep(0.000001, 1),rectBounds.TrawlSedRemin[,3],rectBounds.TrawlSedRemin[,4],col="lightgrey",lty = 0)

# superimpose lines

segments(1-.35,
         TrawlSedRemin$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
         x1 = 1-.05)

# superimpose text refs

text(1-0.2,
     TrawlSedRemin$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
     TrawlSedRemin$refNum, cex = 0.6)

# non-log plot

plot(rep(1,2), TrawlSedRemin$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1., type = "p", col = 0, xlim=c(0.5,5.5), ylim=c(0,17),xlab = "Flux", ylab = "Seafloor sediment remineralized (Pg C yr-1)", yaxs="i")

rectBounds.TrawlSedRemin <- data.frame(c(1-.35),
                                        rep(0, 1),
                                        c(1-.05),
                                        c(max(unlist(c(TrawlSedRemin[c(9,11)])), na.rm = T)))

rect(rectBounds.TrawlSedRemin[,1],rep(0,1),rectBounds.TrawlSedRemin[,3],rectBounds.TrawlSedRemin[,4],col="lightgrey",lty = 0)

# superimpose lines

segments(1-.35,
         TrawlSedRemin$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
         x1 = 1-.05)

# superimpose text refs

text(1-0.2,
     TrawlSedRemin$Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.,
     TrawlSedRemin$refNum, cex = 0.6)

# plot of overall ocean C storage

# subset to relevant
bulkOceanC.rawDat <- rawCdat[rawCdat$Overall.category=="Bulk C storage reservoirs" | (rawCdat$Overall.category=="Fish" & rawCdat$Biomass.or.export=="Biomass"),]

# create target DF for plot data

bulkOceanC.plotDat <- data.frame(matrix(ncol = 4, nrow = 9))
colnames(bulkOceanC.plotDat) <- c("Reservoir","Storage (Pg C)","Reference","refNum")
bulkOceanC.bulkCats <- c("Living ocean biota","DOC","DIC","Anthropogenic DIC","Sediment C_org (top 1 m)","Sediment C_org (1-5 cm)",
                         "Sediment C_org (5-30 cm)","Sediment C_org (30-100 cm)","Sediment C_org (all reactive sediments)")
bulkOceanC.plotDat$Reservoir <- bulkOceanC.bulkCats
bulkOceanC.plotDat$Reference <- as.character(bulkOceanC.plotDat$Reference)

# populate plot data DF from raw data

bulkOceanC.plotDat[bulkOceanC.plotDat$Reservoir=="Living ocean biota",c(2:4)] <-
  c(max(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Fish" & bulkOceanC.rawDat$Biomass.or.export=="Biomass",c("Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.","Current.hi.bound..Pg.C.or.Pg.C.yr.1.")], na.rm = T),
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Fish" & bulkOceanC.rawDat$Biomass.or.export=="Biomass",c("Reference")], collapse = "; "),
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Fish" & bulkOceanC.rawDat$Biomass.or.export=="Biomass",c("refNum")], collapse = ", "))

bulkOceanC.plotDat[bulkOceanC.plotDat$Reservoir=="DOC",c(2:4)] <-
  c(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Bulk C storage reservoirs" & bulkOceanC.rawDat$Subcategory=="Ocean DOC",c("Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.")],
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Bulk C storage reservoirs" & bulkOceanC.rawDat$Subcategory=="Ocean DOC",c("Reference")], collapse = "; "),
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Bulk C storage reservoirs" & bulkOceanC.rawDat$Subcategory=="Ocean DOC",c("refNum")], collapse = ", "))

bulkOceanC.plotDat[bulkOceanC.plotDat$Reservoir=="DIC",c(2:4)] <-
  c(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Bulk C storage reservoirs" & bulkOceanC.rawDat$Subcategory=="Ocean DIC" & bulkOceanC.rawDat$Biomass.or.export!="Anthro DIC",c("Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.")],
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Bulk C storage reservoirs" & bulkOceanC.rawDat$Subcategory=="Ocean DIC" & bulkOceanC.rawDat$Biomass.or.export!="Anthro DIC",c("Reference")], collapse = "; "),
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Bulk C storage reservoirs" & bulkOceanC.rawDat$Subcategory=="Ocean DIC" & bulkOceanC.rawDat$Biomass.or.export!="Anthro DIC",c("refNum")], collapse = ", "))

bulkOceanC.plotDat[bulkOceanC.plotDat$Reservoir=="Anthropogenic DIC",c(2:4)] <-
  c(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Bulk C storage reservoirs" & bulkOceanC.rawDat$Biomass.or.export=="Anthro DIC",c("Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.")],
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Bulk C storage reservoirs" & bulkOceanC.rawDat$Biomass.or.export=="Anthro DIC",c("Reference")], collapse = "; "),
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Overall.category=="Bulk C storage reservoirs" & bulkOceanC.rawDat$Biomass.or.export=="Anthro DIC",c("refNum")], collapse = ", "))

bulkOceanC.plotDat[bulkOceanC.plotDat$Reservoir=="Sediment C_org (top 1 m)",c(2:4)] <-
  c(sum(bulkOceanC.rawDat[bulkOceanC.rawDat$Subcategory=="Ocean sediment C_org - top 1 m",
                          c("Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.")]),
    paste(unique(bulkOceanC.rawDat[bulkOceanC.rawDat$Subcategory=="Ocean sediment C_org - top 1 m",
                            c("Reference")], collapse = "; ")),
    paste(unique(bulkOceanC.rawDat[bulkOceanC.rawDat$Subcategory=="Ocean sediment C_org - top 1 m",
                                   c("refNum")], collapse = "; ")))

bulkOceanC.plotDat[bulkOceanC.plotDat$Reservoir=="Sediment C_org (1-5 cm)",c(2:4)] <-
  c(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="1-5 cm",c("Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.")],
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="1-5 cm",c("Reference")], collapse = "; "),
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="1-5 cm",c("refNum")], collapse = ", "))

bulkOceanC.plotDat[bulkOceanC.plotDat$Reservoir=="Sediment C_org (5-30 cm)",c(2:4)] <-
  c(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="5-30 cm",c("Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.")],
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="5-30 cm",c("Reference")], collapse = "; "),
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="5-30 cm",c("refNum")], collapse = ", "))

bulkOceanC.plotDat[bulkOceanC.plotDat$Reservoir=="Sediment C_org (30-100 cm)",c(2:4)] <-
  c(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="30-100 cm",c("Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.")],
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="30-100 cm",c("Reference")], collapse = "; "),
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="30-100 cm",c("refNum")], collapse = ", "))

bulkOceanC.plotDat[bulkOceanC.plotDat$Reservoir=="Sediment C_org (all reactive sediments)",c(2:4)] <-
  c(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="All reactive sediments",c("Current.quantity.or.height.of.whaling.fishing.minimum..Pg.C.or.Pg.C.yr.1.")],
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="All reactive sediments",c("Reference")], collapse = "; "),
    paste(bulkOceanC.rawDat[bulkOceanC.rawDat$Biomass.or.export=="All reactive sediments",c("refNum")], collapse = ", "))

bulkOceanC.plotDat$`Storage (Pg C)` <- as.numeric(bulkOceanC.plotDat$`Storage (Pg C)`)

masterSlices <- as.numeric(c(bulkOceanC.plotDat$`Storage (Pg C)`[c(1:3,9)]))
sedSlices <- as.numeric(c(bulkOceanC.plotDat$`Storage (Pg C)`[c(1:3)],
                          (bulkOceanC.plotDat$`Storage (Pg C)`[9]*
                            (bulkOceanC.plotDat$`Storage (Pg C)`[6]/bulkOceanC.plotDat$`Storage (Pg C)`[5])),
                          (bulkOceanC.plotDat$`Storage (Pg C)`[9]*
                            (bulkOceanC.plotDat$`Storage (Pg C)`[7]/bulkOceanC.plotDat$`Storage (Pg C)`[5])),
                          (bulkOceanC.plotDat$`Storage (Pg C)`[9]*
                            (bulkOceanC.plotDat$`Storage (Pg C)`[8]/bulkOceanC.plotDat$`Storage (Pg C)`[5]))
                          ))
anthroSlices <- c(bulkOceanC.plotDat$`Storage (Pg C)`[c(1:2)],
                  bulkOceanC.plotDat$`Storage (Pg C)`[3]-bulkOceanC.plotDat$`Storage (Pg C)`[4],
                  bulkOceanC.plotDat$`Storage (Pg C)`[4],
                  bulkOceanC.plotDat$`Storage (Pg C)`[9])

masterLabels <- c(paste(c(bulkOceanC.plotDat$Reservoir[c(1:3,9)])))

cols <- c(brewer.pal(n = 3, name = 'Greens'),
          brewer.pal(n = 3, name = 'Reds'),
          brewer.pal(n = 9, name = 'Blues'),
          brewer.pal(n = 9, name = 'Purples'))

cols.base <- cols[c(3,6,13,19)]
density.base <- c(NA,NA,NA,50)
cols.seds <- c(rep("#FF000000",3),cols[18])
cols.sedBD <- c(rep("#FF000000",3),cols[c(20,22,24)])
cols.anthro <- c(rep("#FF000000",3),cols[11],rep("#FF000000",1))

# ratio for sediments

sedRat <- as.numeric(bulkOceanC.plotDat$`Storage (Pg C)`[5])/as.numeric(bulkOceanC.plotDat$`Storage (Pg C)`[9])
slicerat <- masterSlices[4]/sum(masterSlices)

pie(masterSlices, masterLabels, col = cols.base, density = density.base, radius = -1)
par(new=TRUE)
pie(masterSlices, col = cols.seds, radius = -sqrt(sedRat))
par(new=TRUE)
pie(sedSlices, col = cols.sedBD, radius = -sqrt(sedRat))
par(new=TRUE)
pie(anthroSlices, col = cols.anthro, radius = -1)

