# Libraries ####
library(dplR)
library(ggplot2)
library(reshape2)

# Reproducting the bug ####

# Load bugged series
ngxp_bugged <- read.rwl("ngxp_bugged.rwl")

# Attempt to detrend, showing plots
std_series <- detrend(ngxp_bugged, make.plot=TRUE)

# Reconstruct values of the detrending curve
# Raw data / detrending curve = standardized
# Raw data /  standardized = detrending curve
detrending_curve <- ngxp_bugged / std_series[['101670']]$ModNegExp

# Show raw data and detrending curve in a clearer graph
plot_df <- melt(data.frame(Raw=ngxp_bugged[[1]], Detrending=detrending_curve[[1]], Year=as.numeric(rownames(ngxp_bugged))), id.vars="Year")
ggplot(data=plot_df,aes(x=Year, y=value, colour=variable)) + geom_line(size=2) + theme_bw() + geom_hline(x=0) + ylab("Ring Width (mm)") + scale_colour_discrete(name="Curve")

# Opening up the relevant sections 