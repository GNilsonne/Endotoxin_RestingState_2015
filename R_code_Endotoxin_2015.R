# Script to analyse behavioral and physiological data in the endotoxin resting state project 2015
# Gustav Nilsonne 2015-04-04
# Code is free

# REQUIRE PACKAGES
require(RCurl) # To read data from GitHub

# READ DATA
url <- getURL("https://raw.githubusercontent.com/GNilsonne/Endotoxin_RestingState_2015/master/BehavioralAndPhysiologicalData.csv", ssl.verifypeer = FALSE)
data <- read.table(text = url, header = TRUE)

# MAKE DESCRIPTIVE PLOTS FOR THE DIFFERENT VARIABLES
# 1. Heart rate
HR_data <- data[, c("subject", "group", "HR_0", "HR_1", "HR_2", "HR_3", "HR_4")]
HR_data <- reshape(HR_data, direction = 'long', varying = 3:7, v.names = "HR", timevar = "time", sep = "_")
HR_agg <- aggregate(HR_data, by=list(as.logical(as.integer(HR_data$group)-1), as.numeric(HR_data$time)), FUN=mean, na.rm=TRUE)

hist(HR_data$HR)
hist(sqrt(HR_data$HR))
hist(1/(HR_data$HR))
hist(log(HR_data$HR))
# Note: Data should probably be analysed after inverse-transformation to better approximate a normal distribution. Jennings 1981 provides guidelines for heart rate analysis, but says nothing about transformation.
# I plot data on the original scale nonetheless, since that will be easier to understand and interpret.

pdf("Fig_HR.pdf")
par(mar=c(6, 5, 1, 1))
plot(HR ~ time, data = HR_data, pch = "", xlab = "Time of measurement", ylab = "Beats per minute", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0))
for (i in unique(HR_data$subject[HR_data$group == "Placebo"])){
  lines(HR ~ time, data = HR_data[HR_data$subject == i, ], col = "lightblue")
}
for (i in unique(HR_data$subject[HR_data$group == "Endotoxin"])){
  lines(HR ~ time, data = HR_data[HR_data$subject == i, ], col = "pink")
}
lines(HR ~ time, data = HR_agg[HR_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(HR ~ time, data = HR_agg[HR_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
dev.off()

# 2. Body temperature