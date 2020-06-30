# Copyright (c) 2018, Sam Lambrick
# All rights reserverd
# This file is part of the SHeM Ray Tracing Simulation, subject to the
# GNU/GPL-3.0-or-later.
#
# This script produces line and scatter plots for 1D scans. Produces plots 
# comparing the three different contributions: single scattering, multiple 
# scattering and scattering from the diffuse beam. Also produces plots of just 
# the total detections and ratios between different detections.
#
# The location and name of the data file (.csv), the direction of the line scan
# and whether or not to zero the y axis of the plots can be given as either 
# parameters at the top of the script. To call with command line arguments
# specify zeroed as either true/false:
#  Rscript --vanilla line_plotting.R directory scan_direction zeroed
#
# Plots are produced to be a size useful for display on a computer or for
# printing. They are not designed to be publication quality. Modify the ggsave
# functions at the bottom to change the size of the figure (default units in 
# inches). The font size is set to 10pt, this may also be changed.

# Use the tidyverse
suppressPackageStartupMessages(library(tidyverse))

#---------------------------- Arguments/parameters ----------------------------#

# Take three command line arguments or non, if non the arguments are specified
# below
args <- commandArgs(trailingOnly=TRUE)

# Specify the directory the data is in
directory <- "002_rodSamples/03_noGrid_uniform"

# Specify the direction of the line scan
scan_direction <- "y"

# Make the plots go all the way to zero on the y axis
zeroed <- TRUE


# Read command line arguments if there are any
if (length(args) == 3) {
    directory <- args[1]
    scan_direction <- args[2]
    zeroed <- args[3]
} else if(length(args) != 0) {
    stop("Three (or zero) arguments must be supplied to line_plotting.R")
}

# Check the second and third argument
if (scan_direction != "x" && scan_direction != "y" && scan_direction != "z") {
    stop("The scan direction must be 'x' 'y' or 'z'")
}

if (zeroed == "TRUE" || zeroed == "True" || zeroed == "true") {
    zeroed <- TRUE
} else if (zeroed == "FALSE" || zeroed == "False" || zeroed == "false") {
    zeroed <- FALSE
} else {
    stop("The third argument (zero the plots) must be true or false.")
}

# Check to see if the file exists
fname <- paste(directory, "/data_for_plotting.csv", sep="")
if (!file.exists(fname)) {
    stop("The data file does not exsist.")
}

#-------------------------------- Data analysis -------------------------------#


# Read in the data
theData <- suppressMessages(read_csv(file=fname))

# Rescale the postions
if (scan_direction == "y") {
    # We plot as a function of the distance from the pinhole plate in the case
    # of a y direction scan
    theData$Positions <- theData$Positions + 2.121
} else {
    # Plot in micrometers instead of mm
    theData$Positions <- theData$Positions*1000
}

# Get the correct variables for plotting
theData$Effuse <- theData$EffuseSingle + theData$EffuseMultiple
theData$All <- theData$Single + theData$Multiple + theData$Effuse
theData$snglOverMltp <- theData$Single/theData$Multiple
theData$proportionEffuse <- theData$Effuse/theData$All
theData$proportionSingle <- theData$Single/theData$All

# Trandform the tibble into long format for plotting
if (sum(theData$Effuse) == 0) {
    # If there was no effuse beam then don't plot the effuse components
    longData <- gather(theData, "Scattering", "Counts", Single, Multiple, All,
                       -Position, factor_key=TRUE)
    effuse_present <- FALSE
} else {
    # Include the effuse contribution and create a seperate data with
    # information about the multiple/single scattering of the effuse beam
    longData <- gather(theData, "Scattering", "Counts", Single, Multiple, Effuse,
                       All, -Positions, factor_key=TRUE)
    longEverything <- gather(theData, "Scattering", "Counts", Single, Multiple, 
                             EffuseSingle, EffuseMultiple, Effuse, All,
                             -Positions, factor_key=TRUE)
    effuseData <- gather(theData, "Scattering", "Counts", EffuseSingle, 
                         EffuseMultiple, Effuse, -Positions, factor_key=TRUE)
    effuse_present <- TRUE
}

#-------------------------------- Produce plots -------------------------------#

# Simple plot of the total line scan
pltAll <- ggplot(data=theData, aes(x=Positions, y=All)) +
    ylab("Number of counts")

# A plot comparing the different contributions
if (effuse_present) {
    theLabs <- c("Single", "Multiple", "Effuse", "Total")
} else {
    theLabs <- c("Single", "Multiple", "Total")
}
pltCmp <- ggplot(data=longData, aes(x=Positions, y=Counts, colour=Scattering)) +
    scale_colour_hue(labels=theLabs) + 
    ylab("Number of counts") +
    theme(legend.title=element_blank())

# Plot the ratio of single scattering to multiple scatteing
pltRatio <- ggplot(data=theData, aes(x=Positions, y=snglOverMltp)) +
    ylab("Single/Multiple scattering") + 
    geom_point(size=0.9, colour="firebrick4")

if (effuse_present) {
    # Plot the proportion of detected rays that are from the effuse beam
    pltEffuseRatio <- ggplot(data=theData, aes(x=Positions, y=proportionEffuse)) +
        ylab("Proportion of detections\n from effuse beam") +
        geom_point(size=0.9, colour="firebrick4") +
        ylim(0, max(theData$proportionEffuse))
    
    # Do a plot of all the contributions, with single and multiple effuse
    # scattering seperate, need to create a new tibble for this and create a 
    # new set of data
    moreLabs <- c("Single", "Multiple", "Effuse single", "Effuse multiple",
                  "Total effuse", "Total")
    pltEverything <- ggplot(longEverything, aes(x=Positions, y=Counts, 
                            colour=Scattering)) +
        scale_colour_hue(labels=moreLabs) +
        ylab("Number of counts") + 
        geom_line() +
        theme(legend.title=element_blank())
    
    # Do a plot of the effuse contributions
    effuseLabs <- c("Effuse single", "Effuse multiple", "Total effuse")
    pltEffuseComp <- ggplot(effuseData, aes(x=Positions, y=Counts, colour=Scattering)) +
        scale_colour_hue(labels=effuseLabs) +
        ylab("Number of counts") +
        geom_line() +
        theme(legend.title=element_blank())
}

# Plot the proportion of detected rays that are from single scattering
pltSingle <- ggplot(data=theData, aes(x=Positions, y=proportionSingle)) +
    ylab("Proportion single scattering") +
    geom_line(colour="dodgerblue4") + 
    ylim(0,1)

#------------------------------- Annotate plots -------------------------------#
# A function that adds the x labels to the plots
addLab <- function(plt) {
    if (scan_direction == "y") {
        plt <- plt + xlab("Distance from pinhole plate/mm")
    }
    if (scan_direction == "x") {
        plt <- plt + xlab(expression(paste("x", "/", mu, "m")))
    }
    if (scan_direction == "z") {
        plt <- plt + xlab(expression(paste("z", "/", mu, "m")))
    }
    return(plt)
}

# Add the x labels to each of the plots
pltCmp <- addLab(pltCmp)
pltAll <- addLab(pltAll)
pltRatio <- addLab(pltRatio)
pltSingle <- addLab(pltSingle)

if (effuse_present) {
    pltEffuseRatio <- addLab(pltEffuseRatio)
    pltEverything <- addLab(pltEverything)
    pltEffuseComp <- addLab(pltEffuseComp)
}

# Force the y scale to go to zero
if (zeroed) {
    pltCmp <- pltCmp + ylim(0, max(theData$All))
    pltAll <- pltAll + ylim(0, max(theData$All))
}

# Change the font size in the plots to 10 (from default of 11)
pltCmp <- pltCmp + theme(text=element_text(size=10))
pltAll <- pltAll + theme(text=element_text(size=10))
pltRatio <- pltRatio + theme(text=element_text(size=10))
if (effuse_present) {
    pltEffuseRatio <- pltEffuseRatio + theme(text=element_text(size=10))
    pltEverything <- pltEverything + theme(text=element_text(size=10))
    pltEffuseComp <- pltEffuseComp + theme(text=element_text(size=10))
}
pltSingle <- pltSingle + theme(text=element_text(size=10))

# We make both line and scatter plots for the two main comparisons
pltCmp1 <- pltCmp + geom_line()
pltCmp2 <- pltCmp + geom_point(size=0.9)

pltAll1 <- pltAll + geom_line()
pltAll2 <- pltAll + geom_point(size=0.9)

#--------------------------------- Save plots ---------------------------------#
# Save the plots to the directory provided

ggsave(paste(directory, "/comparison_plot1.eps", sep=""), pltCmp1, 
       width=6, height=4)
ggsave(paste(directory, "/comparison_plot2.eps", sep=""), pltCmp2,
       width=6, height=4)
ggsave(paste(directory, "/all_plot1.eps", sep=""), pltAll1, 
       width=4.5, height=4)
ggsave(paste(directory, "/all_plot2.eps", sep=""), pltAll2,
       width=4.5, height=4)
ggsave(paste(directory, "/single_vs_multiple.eps", sep=""), pltRatio, 
       width=3.5, height=3)
if (effuse_present) {
    ggsave(paste(directory, "/proportion_effuse.eps", sep=""), pltEffuseRatio,
           width=3.5, height=3)
    ggsave(paste(directory, "/everything_comparison.eps", sep=""), 
           pltEverything, width=6, height=4)
    ggsave(paste(directory, "/effuse_contributions.eps", sep=""), pltEffuseComp,
           width=6, height=4)
}
ggsave(paste(directory, "/proportion_single.eps", sep=""), pltSingle, 
       width=3.5, height=3)

