# Script to analyse behavioral and physiological data in the endotoxin resting state project 2015
# Gustav Nilsonne 2015-04-04
# Code is free

# REQUIRE PACKAGES
require(RCurl) # To read data from GitHub
require(compute.es) # To compute effect sizes

# READ DATA
url <- getURL("https://raw.githubusercontent.com/GNilsonne/Endotoxin_RestingState_2015/master/BehavioralAndPhysiologicalData.csv", ssl.verifypeer = FALSE)
data <- read.csv(text = url, header = TRUE)

# MAKE DESCRIPTIVE ANALYSES AND PLOTS FOR THE DIFFERENT VARIABLES
# 1. Heart rate
HR_data <- data[, c("subject", "group", "HR_0", "HR_1", "HR_2", "HR_3", "HR_4")]
HR_data <- reshape(HR_data, direction = 'long', varying = 3:7, v.names = "HR", timevar = "time", sep = "_")
HR_data <- HR_data[, -5]
HR_agg <- aggregate(HR_data, by=list(as.logical(as.integer(HR_data$group)-1), as.numeric(HR_data$time)), FUN=mean, na.rm=TRUE)

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

data$HR_responses <- rowMeans(cbind(data$HR_2, data$HR_3), na.rm = T) - data$HR_0 # We use the average of the 3rd and 4th measures to estimate effect of endotoxin
boxplot(HR_responses ~ group, data = data, frame.plot = F, main = "Heart rate difference", ylab = "Beats per minute")
ttest_HR <- t.test(HR_responses ~ group, data = data)
es_HR <- tes(t = ttest_HR$statistic, n.1 = length(data$HR_responses[data$group == "Endotoxin"]), n.2 = length(data$HR_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_HR$d, " [", es_HR$l.d, ", ", es_HR$u.d, "]", sep = ""), bty = "n")


# 2. Body temperature
temp_data <- data[, c("subject", "group", "temp_0", "temp_1", "temp_2", "temp_3", "temp_4")]
temp_data <- reshape(temp_data, direction = 'long', varying = 3:7, v.names = "temp", timevar = "time", sep = "_")
temp_data <- temp_data[, -5]
temp_agg <- aggregate(temp_data, by=list(as.logical(as.integer(temp_data$group)-1), as.numeric(temp_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_temp.pdf")
par(mar=c(6, 5, 1, 1))
plot(temp ~ time, data = temp_data, pch = "", xlab = "Time of measurement", ylab = "Degrees C", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0))
for (i in unique(temp_data$subject[temp_data$group == "Placebo"])){
  lines(temp ~ time, data = temp_data[temp_data$subject == i, ], col = "lightblue")
}
for (i in unique(temp_data$subject[temp_data$group == "Endotoxin"])){
  lines(temp ~ time, data = temp_data[temp_data$subject == i, ], col = "pink")
}
lines(temp ~ time, data = temp_agg[temp_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(temp ~ time, data = temp_agg[temp_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
dev.off()

data$temp_responses <- rowMeans(cbind(data$temp_2, data$temp_3), na.rm = T) - data$temp_0 # We use the average of the 3rd and 4th measures to estimate effect of endotoxin
boxplot(temp_responses ~ group, data = data, frame.plot = F, main = "Temperature")
ttest_temp <- t.test(temp_responses ~ group, data = data)
es_temp <- tes(t = ttest_temp$statistic, n.1 = length(data$temp_responses[data$group == "Endotoxin"]), n.2 = length(data$temp_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_temp$d), bty = "n")




# X. Headache
headache_data <- data[, c("subject", "group", "headache_0", "headache_1", "headache_2", "headache_3", "headache_4")]
headache_data <- reshape(headache_data, direction = 'long', varying = 3:7, v.names = "headache", timevar = "time", sep = "_")
headache_data <- headache_data[, -5]
headache_agg <- aggregate(headache_data, by=list(as.logical(as.integer(headache_data$group)-1), as.numeric(headache_data$time)), FUN=mean, na.rm=TRUE)

hist(headache_data$headache)

pdf("Fig_headache.pdf")
par(mar=c(6, 5, 1, 1))
plot(headache ~ time, data = headache_data, pch = "", xlab = "Time of measurement", ylab = "Rated pain", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), ylim = c(0, 7))
for (i in unique(headache_data$subject[headache_data$group == "Endotoxin"])){
  lines(headache ~ time, data = headache_data[headache_data$subject == i, ], col = "pink")
}
for (i in unique(headache_data$subject[headache_data$group == "Placebo"])){
  lines(headache ~ time, data = headache_data[headache_data$subject == i, ], col = "lightblue")
}
lines(headache ~ time, data = headache_agg[headache_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(headache ~ time, data = headache_agg[headache_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
dev.off()

# X. Back pain
back_pain_data <- data[, c("subject", "group", "back_pain_0", "back_pain_1", "back_pain_2", "back_pain_3", "back_pain_4")]
back_pain_data <- reshape(back_pain_data, direction = 'long', varying = 3:7, v.names = "back_pain", timevar = "time", sep = "_")
back_pain_data <- back_pain_data[, -5]
back_pain_agg <- aggregate(back_pain_data, by=list(as.logical(as.integer(back_pain_data$group)-1), as.numeric(back_pain_data$time)), FUN=mean, na.rm=TRUE)

hist(back_pain_data$back_pain)

pdf("Fig_back_pain.pdf")
par(mar=c(6, 5, 1, 1))
plot(back_pain ~ time, data = back_pain_data, pch = "", xlab = "Time of measurement", ylab = "Rated pain", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), ylim = c(0, 7))
for (i in unique(back_pain_data$subject[back_pain_data$group == "Endotoxin"])){
  lines(back_pain ~ time, data = back_pain_data[back_pain_data$subject == i, ], col = "pink")
}
for (i in unique(back_pain_data$subject[back_pain_data$group == "Placebo"])){
  lines(back_pain ~ time, data = back_pain_data[back_pain_data$subject == i, ], col = "lightblue")
}
lines(back_pain ~ time, data = back_pain_agg[back_pain_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(back_pain ~ time, data = back_pain_agg[back_pain_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
dev.off()

# X. Nausea
nausea_data <- data[, c("subject", "group", "nausea_0", "nausea_1", "nausea_2", "nausea_3", "nausea_4")]
nausea_data <- reshape(nausea_data, direction = 'long', varying = 3:7, v.names = "nausea", timevar = "time", sep = "_")
nausea_data <- nausea_data[, -5]
nausea_agg <- aggregate(nausea_data, by=list(as.logical(as.integer(nausea_data$group)-1), as.numeric(nausea_data$time)), FUN=mean, na.rm=TRUE)

hist(nausea_data$nausea)

pdf("Fig_nausea.pdf")
par(mar=c(6, 5, 1, 1))
plot(nausea ~ time, data = nausea_data, pch = "", xlab = "Time of measurement", ylab = "Rated pain", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), ylim = c(0, 7))
for (i in unique(nausea_data$subject[nausea_data$group == "Endotoxin"])){
  lines(nausea ~ time, data = nausea_data[nausea_data$subject == i, ], col = "pink")
}
for (i in unique(nausea_data$subject[nausea_data$group == "Placebo"])){
  lines(nausea ~ time, data = nausea_data[nausea_data$subject == i, ], col = "lightblue")
}
lines(nausea ~ time, data = nausea_agg[nausea_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(nausea ~ time, data = nausea_agg[nausea_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
dev.off()

# WRITE DATA IN LONG FORMAT
PhysiologicalDataLong <- merge(HR_data, by = c("subject", "time"))

SelfRatedDataLong <- merge(headache_data, back_pain_data, by = c("subject", "time"))
SelfRatedDataLong <- merge(SelfRatedDataLong, nausea_data, by = c("subject", "time"))
SelfRatedDataLong <- SelfRatedDataLong[, c("subject", "group", "time", "headache", "back_pain", "nausea")]
SelfRatedDataLong <- SelfRatedDataLong[order(SelfRatedDataLong$subject), ]
SelfRatedDataLong$time <- SelfRatedDataLong$time -1
write.csv(SelfRatedDataLong, file = "C:/Users/Gustav Nilsonne/Box Sync/Gustavs_arbete/Pek/Endotoxin resting state/Endotoxin_RestingState_2015/SelfRatedDataLong.csv", row.names=FALSE)
