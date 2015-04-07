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
boxplot(temp_responses ~ group, data = data, frame.plot = F, main = "Temperature difference", ylab = "Degrees C")
ttest_temp <- t.test(temp_responses ~ group, data = data)
es_temp <- tes(t = ttest_temp$statistic, n.1 = length(data$temp_responses[data$group == "Endotoxin"]), n.2 = length(data$temp_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_temp$d, " [", es_temp$l.d, ", ", es_temp$u.d, "]", sep = ""), bty = "n")

# 3. Blood pressure (systolic)
# Diastolic blood pressure is not analysed
# One data point recorded a systolic blood pressure of 55 mmHg (participant 27, measurement 3). This recording was judged to be erroneous and was excluded.
data$BP_systolic_2[data$BP_systolic_2 == 55] <- NA
BP_systolic_data <- data[, c("subject", "group", "BP_systolic_0", "BP_systolic_1", "BP_systolic_2", "BP_systolic_3", "BP_systolic_4")]
BP_systolic_data <- reshape(BP_systolic_data, direction = 'long', varying = 3:7, v.names = "BP_systolic", timevar = "time", sep = "_")
BP_systolic_data <- BP_systolic_data[, -5]
BP_systolic_agg <- aggregate(BP_systolic_data, by=list(as.logical(as.integer(BP_systolic_data$group)-1), as.numeric(BP_systolic_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_BP_systolic.pdf")
par(mar=c(6, 5, 1, 1))
plot(BP_systolic ~ time, data = BP_systolic_data, pch = "", xlab = "Time of measurement", ylab = "mm Hg", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0))
for (i in unique(BP_systolic_data$subject[BP_systolic_data$group == "Placebo"])){
  lines(BP_systolic ~ time, data = BP_systolic_data[BP_systolic_data$subject == i, ], col = "lightblue")
}
for (i in unique(BP_systolic_data$subject[BP_systolic_data$group == "Endotoxin"])){
  lines(BP_systolic ~ time, data = BP_systolic_data[BP_systolic_data$subject == i, ], col = "pink")
}
lines(BP_systolic ~ time, data = BP_systolic_agg[BP_systolic_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(BP_systolic ~ time, data = BP_systolic_agg[BP_systolic_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
dev.off()

data$BP_systolic_responses <- data$BP_systolic_2 - data$BP_systolic_0 # We use only the 3rd measure to estimate effect of endotoxin
boxplot(BP_systolic_responses ~ group, data = data, frame.plot = F, main = "Blood pressure difference", ylab = "Degrees C")
ttest_BP_systolic <- t.test(BP_systolic_responses ~ group, data = data)
es_BP_systolic <- tes(t = ttest_BP_systolic$statistic, n.1 = length(data$BP_systolic_responses[data$group == "Endotoxin"]), n.2 = length(data$BP_systolic_responses[data$group == "Placebo"]))
legend("topleft", legend = paste("d = ", es_BP_systolic$d, " [", es_BP_systolic$l.d, ", ", es_BP_systolic$u.d, "]", sep = ""), bty = "n")


# 4. IL-6
IL6_data <- data[, c("subject", "group", "IL6_1", "IL6_2", "IL6_3", "IL6_4")]
IL6_data <- reshape(IL6_data, direction = 'long', varying = 3:6, v.names = "IL6", timevar = "time", sep = "_")
low_values_IL6 <- which(IL6_data$IL6 < 0.9) # Values under the assay detection limit were conservatively set to that limit
IL6_data$IL6[low_values_IL6] <- 0.9
IL6_data$logIL6 <- log(IL6_data$IL6)
IL6_agg <- aggregate(IL6_data, by=list(as.logical(as.integer(IL6_data$group)-1), as.numeric(IL6_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_IL6.pdf")
par(mar=c(6, 5, 1, 1))
plot(logIL6 ~ time, data = IL6_data, pch = "", xlab = "Time of measurement", ylab = "log IL6, pg/ml", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(IL6_data$subject[IL6_data$group == "Placebo"])){
  lines(logIL6 ~ time, data = IL6_data[IL6_data$subject == i, ], col = "lightblue")
}
for (i in unique(IL6_data$subject[IL6_data$group == "Endotoxin"])){
  lines(logIL6 ~ time, data = IL6_data[IL6_data$subject == i, ], col = "pink")
}
lines(logIL6 ~ time, data = IL6_agg[IL6_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(logIL6 ~ time, data = IL6_agg[IL6_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$logIL6_1 <- data$IL6_1
data$logIL6_1[data$logIL6_1 < 0.9] <- 0.9
data$logIL6_1 <- log(data$logIL6_1)
data$logIL6_2 <- data$IL6_2
data$logIL6_2[data$logIL6_2 < 0.9] <- 0.9
data$logIL6_2 <- log(data$logIL6_2)
data$logIL6_3 <- data$IL6_3
data$logIL6_3[data$logIL6_3 < 0.9] <- 0.9
data$logIL6_3 <- log(data$logIL6_3)
data$IL6_responses <- rowMeans(cbind(data$logIL6_2, data$logIL6_3), na.rm = T) - data$logIL6_1 # We use the mean of the 2nd and 3rd measures to estimate effect of endotoxin
boxplot(IL6_responses ~ group, data = data, frame.plot = F, main = "IL6 difference", ylab = "log pg/ml")
ttest_IL6 <- t.test(IL6_responses ~ group, data = data)
es_IL6 <- tes(t = ttest_IL6$statistic, n.1 = length(data$IL6_responses[data$group == "Endotoxin"]), n.2 = length(data$IL6_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_IL6$d, " [", es_IL6$l.d, ", ", es_IL6$u.d, "]", sep = ""), bty = "n")


# 5. IL-8
IL8_data <- data[, c("subject", "group", "IL8_1", "IL8_2", "IL8_3", "IL8_4")]
IL8_data <- reshape(IL8_data, direction = 'long', varying = 3:6, v.names = "IL8", timevar = "time", sep = "_")
low_values_IL8 <- which(IL8_data$IL8 < 0.4) # Values under the assay detection limit were conservatively set to that limit
IL8_data$IL8[low_values_IL8] <- 0.4
IL8_data$logIL8 <- log(IL8_data$IL8)
IL8_agg <- aggregate(IL8_data, by=list(as.logical(as.integer(IL8_data$group)-1), as.numeric(IL8_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_IL8.pdf")
par(mar=c(6, 5, 1, 1))
plot(logIL8 ~ time, data = IL8_data, pch = "", xlab = "Time of measurement", ylab = "log IL8, pg/ml", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(IL8_data$subject[IL8_data$group == "Placebo"])){
  lines(logIL8 ~ time, data = IL8_data[IL8_data$subject == i, ], col = "lightblue")
}
for (i in unique(IL8_data$subject[IL8_data$group == "Endotoxin"])){
  lines(logIL8 ~ time, data = IL8_data[IL8_data$subject == i, ], col = "pink")
}
lines(logIL8 ~ time, data = IL8_agg[IL8_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(logIL8 ~ time, data = IL8_agg[IL8_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$logIL8_1 <- data$IL8_1
data$logIL8_1[data$logIL8_1 < 0.4] <- 0.4
data$logIL8_1 <- log(data$logIL8_1)
data$logIL8_2 <- data$IL8_2
data$logIL8_2[data$logIL8_2 < 0.4] <- 0.4
data$logIL8_2 <- log(data$logIL8_2)
data$logIL8_3 <- data$IL8_3
data$logIL8_3[data$logIL8_3 < 0.4] <- 0.4
data$logIL8_3 <- log(data$logIL8_3)
data$IL8_responses <- rowMeans(cbind(data$logIL8_2, data$logIL8_3), na.rm = T) - data$logIL8_1 # We use the mean of the 2nd and 3rd measures to estimate effect of endotoxin
boxplot(IL8_responses ~ group, data = data, frame.plot = F, main = "IL8 difference", ylab = "log pg/ml")
ttest_IL8 <- t.test(IL8_responses ~ group, data = data)
es_IL8 <- tes(t = ttest_IL8$statistic, n.1 = length(data$IL8_responses[data$group == "Endotoxin"]), n.2 = length(data$IL8_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_IL8$d, " [", es_IL8$l.d, ", ", es_IL8$u.d, "]", sep = ""), bty = "n")


# 6. IL-10
IL10_data <- data[, c("subject", "group", "IL10_1", "IL10_2", "IL10_3", "IL10_4")]
IL10_data <- reshape(IL10_data, direction = 'long', varying = 3:6, v.names = "IL10", timevar = "time", sep = "_")
low_values_IL10 <- which(IL10_data$IL10 < 1.1) # Values under the assay detection limit were conservatively set to that limit
IL10_data$IL10[low_values_IL10] <- 1.1
IL10_data$logIL10 <- log(IL10_data$IL10)
IL10_agg <- aggregate(IL10_data, by=list(as.logical(as.integer(IL10_data$group)-1), as.numeric(IL10_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_IL10.pdf")
par(mar=c(6, 5, 1, 1))
plot(logIL10 ~ time, data = IL10_data, pch = "", xlab = "Time of measurement", ylab = "log IL10, pg/ml", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(IL10_data$subject[IL10_data$group == "Placebo"])){
  lines(logIL10 ~ time, data = IL10_data[IL10_data$subject == i, ], col = "lightblue")
}
for (i in unique(IL10_data$subject[IL10_data$group == "Endotoxin"])){
  lines(logIL10 ~ time, data = IL10_data[IL10_data$subject == i, ], col = "pink")
}
lines(logIL10 ~ time, data = IL10_agg[IL10_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(logIL10 ~ time, data = IL10_agg[IL10_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$logIL10_1 <- data$IL10_1
data$logIL10_1[data$logIL10_1 < 1.1] <- 1.1
data$logIL10_1 <- log(data$logIL10_1)
data$logIL10_2 <- data$IL10_2
data$logIL10_2[data$logIL10_2 < 1.1] <- 1.1
data$logIL10_2 <- log(data$logIL10_2)
data$logIL10_3 <- data$IL10_3
data$logIL10_3[data$logIL10_3 < 1.1] <- 1.1
data$logIL10_3 <- log(data$logIL10_3)
data$IL10_responses <- rowMeans(cbind(data$logIL10_2, data$logIL10_3), na.rm = T) - data$logIL10_1 # We use the mean of the 2nd and 3rd measures to estimate effect of endotoxin
boxplot(IL10_responses ~ group, data = data, frame.plot = F, main = "IL10 difference", ylab = "log pg/ml")
ttest_IL10 <- t.test(IL10_responses ~ group, data = data)
es_IL10 <- tes(t = ttest_IL10$statistic, n.1 = length(data$IL10_responses[data$group == "Endotoxin"]), n.2 = length(data$IL10_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_IL10$d, " [", es_IL10$l.d, ", ", es_IL10$u.d, "]", sep = ""), bty = "n")


# 7. IL-13
IL13_data <- data[, c("subject", "group", "IL13_1", "IL13_2", "IL13_3", "IL13_4")]
IL13_data <- reshape(IL13_data, direction = 'long', varying = 3:6, v.names = "IL13", timevar = "time", sep = "_")
low_values_IL13 <- which(IL13_data$IL13 < 1.3) # Values under the assay detection limit were conservatively set to that limit
IL13_data$IL13[low_values_IL13] <- 1.3
IL13_data$logIL13 <- log(IL13_data$IL13)
IL13_agg <- aggregate(IL13_data, by=list(as.logical(as.integer(IL13_data$group)-1), as.numeric(IL13_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_IL13.pdf")
par(mar=c(6, 5, 1, 1))
plot(logIL13 ~ time, data = IL13_data, pch = "", xlab = "Time of measurement", ylab = "log IL13, pg/ml", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(IL13_data$subject[IL13_data$group == "Placebo"])){
  lines(logIL13 ~ time, data = IL13_data[IL13_data$subject == i, ], col = "lightblue")
}
for (i in unique(IL13_data$subject[IL13_data$group == "Endotoxin"])){
  lines(logIL13 ~ time, data = IL13_data[IL13_data$subject == i, ], col = "pink")
}
lines(logIL13 ~ time, data = IL13_agg[IL13_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(logIL13 ~ time, data = IL13_agg[IL13_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$logIL13_1 <- data$IL13_1
data$logIL13_1[data$logIL13_1 < 1.3] <- 1.3
data$logIL13_1 <- log(data$logIL13_1)
data$logIL13_2 <- data$IL13_2
data$logIL13_2[data$logIL13_2 < 1.3] <- 1.3
data$logIL13_2 <- log(data$logIL13_2)
data$logIL13_3 <- data$IL13_3
data$logIL13_3[data$logIL13_3 < 1.3] <- 1.3
data$logIL13_3 <- log(data$logIL13_3)
data$IL13_responses <- rowMeans(cbind(data$logIL13_2, data$logIL13_3), na.rm = T) - data$logIL13_1 # We use the mean of the 2nd and 3rd measures to estimate effect of endotoxin
boxplot(IL13_responses ~ group, data = data, frame.plot = F, main = "IL13 difference", ylab = "log pg/ml")
ttest_IL13 <- t.test(IL13_responses ~ group, data = data)
es_IL13 <- tes(t = ttest_IL13$statistic, n.1 = length(data$IL13_responses[data$group == "Endotoxin"]), n.2 = length(data$IL13_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_IL13$d, " [", es_IL13$l.d, ", ", es_IL13$u.d, "]", sep = ""), bty = "n")


# 8. TNF
TNF_data <- data[, c("subject", "group", "TNF_1", "TNF_2", "TNF_3", "TNF_4")]
TNF_data <- reshape(TNF_data, direction = 'long', varying = 3:6, v.names = "TNF", timevar = "time", sep = "_")
low_values_TNF <- which(TNF_data$TNF < 0.7) # Values under the assay detection limit were conservatively set to that limit
TNF_data$TNF[low_values_TNF] <- 0.7
TNF_data$logTNF <- log(TNF_data$TNF)
TNF_agg <- aggregate(TNF_data, by=list(as.logical(as.integer(TNF_data$group)-1), as.numeric(TNF_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_TNF.pdf")
par(mar=c(6, 5, 1, 1))
plot(logTNF ~ time, data = TNF_data, pch = "", xlab = "Time of measurement", ylab = "log TNF, pg/ml", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(TNF_data$subject[TNF_data$group == "Placebo"])){
  lines(logTNF ~ time, data = TNF_data[TNF_data$subject == i, ], col = "lightblue")
}
for (i in unique(TNF_data$subject[TNF_data$group == "Endotoxin"])){
  lines(logTNF ~ time, data = TNF_data[TNF_data$subject == i, ], col = "pink")
}
lines(logTNF ~ time, data = TNF_agg[TNF_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(logTNF ~ time, data = TNF_agg[TNF_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$logTNF_1 <- data$TNF_1
data$logTNF_1[data$logTNF_1 < 0.7] <- 0.7
data$logTNF_1 <- log(data$logTNF_1)
data$logTNF_2 <- data$TNF_2
data$logTNF_2[data$logTNF_2 < 0.7] <- 0.7
data$logTNF_2 <- log(data$logTNF_2)
data$logTNF_3 <- data$TNF_3
data$logTNF_3[data$logTNF_3 < 0.7] <- 0.7
data$logTNF_3 <- log(data$logTNF_3)
data$TNF_responses <- rowMeans(cbind(data$logTNF_2, data$logTNF_3), na.rm = T) - data$logTNF_1 # We use the mean of the 2nd and 3rd measures to estimate effect of endotoxin
boxplot(TNF_responses ~ group, data = data, frame.plot = F, main = "TNF difference", ylab = "log pg/ml")
ttest_TNF <- t.test(TNF_responses ~ group, data = data)
es_TNF <- tes(t = ttest_TNF$statistic, n.1 = length(data$TNF_responses[data$group == "Endotoxin"]), n.2 = length(data$TNF_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_TNF$d, " [", es_TNF$l.d, ", ", es_TNF$u.d, "]", sep = ""), bty = "n")


# 9. SQ total
sq_total_data <- data[, c("subject", "group", "sq_total_1", "sq_total_2", "sq_total_4")]
sq_total_data <- reshape(sq_total_data, direction = 'long', varying = 3:5, v.names = "sq_total", timevar = "time", sep = "_")
sq_total_agg <- aggregate(sq_total_data, by=list(as.logical(as.integer(sq_total_data$group)-1), as.numeric(sq_total_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_sq_total.pdf")
par(mar=c(6, 5, 1, 1))
plot(sq_total ~ time, data = sq_total_data, pch = "", xlab = "Time of measurement", ylab = "Score", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(sq_total_data$subject[sq_total_data$group == "Placebo"])){
  lines(sq_total ~ time, data = sq_total_data[sq_total_data$subject == i, ], col = "lightblue")
}
for (i in unique(sq_total_data$subject[sq_total_data$group == "Endotoxin"])){
  lines(sq_total ~ time, data = sq_total_data[sq_total_data$subject == i, ], col = "pink")
}
lines(sq_total ~ time, data = sq_total_agg[sq_total_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(sq_total ~ time, data = sq_total_agg[sq_total_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$sq_total_responses <- data$sq_total_2 - data$sq_total_1 # We use the 2nd vs the 1st measure to estimate effect of endotoxin
boxplot(sq_total_responses ~ group, data = data, frame.plot = F, main = "sq_total difference", ylab = " pg/ml")
ttest_sq_total <- t.test(sq_total_responses ~ group, data = data)
es_sq_total <- tes(t = ttest_sq_total$statistic, n.1 = length(data$sq_total_responses[data$group == "Endotoxin"]), n.2 = length(data$sq_total_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_sq_total$d, " [", es_sq_total$l.d, ", ", es_sq_total$u.d, "]", sep = ""), bty = "n")


# 10. SQ pain
sq_pain_data <- data[, c("subject", "group", "sq_pain_1", "sq_pain_2", "sq_pain_4")]
sq_pain_data <- reshape(sq_pain_data, direction = 'long', varying = 3:5, v.names = "sq_pain", timevar = "time", sep = "_")
sq_pain_agg <- aggregate(sq_pain_data, by=list(as.logical(as.integer(sq_pain_data$group)-1), as.numeric(sq_pain_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_sq_pain.pdf")
par(mar=c(6, 5, 1, 1))
plot(sq_pain ~ time, data = sq_pain_data, pch = "", xlab = "Time of measurement", ylab = "Score", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(sq_pain_data$subject[sq_pain_data$group == "Placebo"])){
  lines(sq_pain ~ time, data = sq_pain_data[sq_pain_data$subject == i, ], col = "lightblue")
}
for (i in unique(sq_pain_data$subject[sq_pain_data$group == "Endotoxin"])){
  lines(sq_pain ~ time, data = sq_pain_data[sq_pain_data$subject == i, ], col = "pink")
}
lines(sq_pain ~ time, data = sq_pain_agg[sq_pain_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(sq_pain ~ time, data = sq_pain_agg[sq_pain_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$sq_pain_responses <- data$sq_pain_2 - data$sq_pain_1 # We use the 2nd vs the 1st measure to estimate effect of endotoxin
boxplot(sq_pain_responses ~ group, data = data, frame.plot = F, main = "sq_pain difference", ylab = " pg/ml")
ttest_sq_pain <- t.test(sq_pain_responses ~ group, data = data)
es_sq_pain <- tes(t = ttest_sq_pain$statistic, n.1 = length(data$sq_pain_responses[data$group == "Endotoxin"]), n.2 = length(data$sq_pain_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_sq_pain$d, " [", es_sq_pain$l.d, ", ", es_sq_pain$u.d, "]", sep = ""), bty = "n")


# 11. SQ fatigue
sq_fatigue_data <- data[, c("subject", "group", "sq_energy_1", "sq_energy_2", "sq_energy_4")]
sq_fatigue_data <- reshape(sq_fatigue_data, direction = 'long', varying = 3:5, v.names = "sq_fatigue", timevar = "time", sep = "_")
sq_fatigue_agg <- aggregate(sq_fatigue_data, by=list(as.logical(as.integer(sq_fatigue_data$group)-1), as.numeric(sq_fatigue_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_sq_fatigue.pdf")
par(mar=c(6, 5, 1, 1))
plot(sq_fatigue ~ time, data = sq_fatigue_data, pch = "", xlab = "Time of measurement", ylab = "Score", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(sq_fatigue_data$subject[sq_fatigue_data$group == "Placebo"])){
  lines(sq_fatigue ~ time, data = sq_fatigue_data[sq_fatigue_data$subject == i, ], col = "lightblue")
}
for (i in unique(sq_fatigue_data$subject[sq_fatigue_data$group == "Endotoxin"])){
  lines(sq_fatigue ~ time, data = sq_fatigue_data[sq_fatigue_data$subject == i, ], col = "pink")
}
lines(sq_fatigue ~ time, data = sq_fatigue_agg[sq_fatigue_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(sq_fatigue ~ time, data = sq_fatigue_agg[sq_fatigue_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$sq_fatigue_responses <- data$sq_energy_2 - data$sq_energy_1 # We use the 2nd vs the 1st measure to estimate effect of endotoxin
boxplot(sq_fatigue_responses ~ group, data = data, frame.plot = F, main = "sq_fatigue difference", ylab = " pg/ml")
ttest_sq_fatigue <- t.test(sq_fatigue_responses ~ group, data = data)
es_sq_fatigue <- tes(t = ttest_sq_fatigue$statistic, n.1 = length(data$sq_fatigue_responses[data$group == "Endotoxin"]), n.2 = length(data$sq_fatigue_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_sq_fatigue$d, " [", es_sq_fatigue$l.d, ", ", es_sq_fatigue$u.d, "]", sep = ""), bty = "n")


# 12. SQ affect
sq_affect_data <- data[, c("subject", "group", "sq_affect_1", "sq_affect_2", "sq_affect_4")]
sq_affect_data <- reshape(sq_affect_data, direction = 'long', varying = 3:5, v.names = "sq_affect", timevar = "time", sep = "_")
sq_affect_agg <- aggregate(sq_affect_data, by=list(as.logical(as.integer(sq_affect_data$group)-1), as.numeric(sq_affect_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_sq_affect.pdf")
par(mar=c(6, 5, 1, 1))
plot(sq_affect ~ time, data = sq_affect_data, pch = "", xlab = "Time of measurement", ylab = "Score", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(sq_affect_data$subject[sq_affect_data$group == "Placebo"])){
  lines(sq_affect ~ time, data = sq_affect_data[sq_affect_data$subject == i, ], col = "lightblue")
}
for (i in unique(sq_affect_data$subject[sq_affect_data$group == "Endotoxin"])){
  lines(sq_affect ~ time, data = sq_affect_data[sq_affect_data$subject == i, ], col = "pink")
}
lines(sq_affect ~ time, data = sq_affect_agg[sq_affect_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(sq_affect ~ time, data = sq_affect_agg[sq_affect_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$sq_affect_responses <- data$sq_affect_2 - data$sq_affect_1 # We use the 2nd vs the 1st measure to estimate effect of endotoxin
boxplot(sq_affect_responses ~ group, data = data, frame.plot = F, main = "sq_affect difference", ylab = " pg/ml")
ttest_sq_affect <- t.test(sq_affect_responses ~ group, data = data)
es_sq_affect <- tes(t = ttest_sq_affect$statistic, n.1 = length(data$sq_affect_responses[data$group == "Endotoxin"]), n.2 = length(data$sq_affect_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_sq_affect$d, " [", es_sq_affect$l.d, ", ", es_sq_affect$u.d, "]", sep = ""), bty = "n")


# 13. Self-rated health now
health_now_data <- data[, c("subject", "group", "health_now_1", "health_now_2", "health_now_4")]
health_now_data <- reshape(health_now_data, direction = 'long', varying = 3:5, v.names = "health_now", timevar = "time", sep = "_")
health_now_agg <- aggregate(health_now_data, by=list(as.logical(as.integer(health_now_data$group)-1), as.numeric(health_now_data$time)), FUN=mean, na.rm=TRUE)

pdf("Fig_health_now.pdf")
par(mar=c(6, 5, 1, 1))
plot(health_now ~ time, data = health_now_data, pch = "", xlab = "Time of measurement", ylab = "Score", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0), xaxt = "n")
for (i in unique(health_now_data$subject[health_now_data$group == "Placebo"])){
  lines(health_now ~ time, data = health_now_data[health_now_data$subject == i, ], col = "lightblue")
}
for (i in unique(health_now_data$subject[health_now_data$group == "Endotoxin"])){
  lines(health_now ~ time, data = health_now_data[health_now_data$subject == i, ], col = "pink")
}
lines(health_now ~ time, data = health_now_agg[health_now_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(health_now ~ time, data = health_now_agg[health_now_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
axis(1, at = c(1, 2, 3, 4), cex.axis = 2, cex.lab = 2)
dev.off()

data$health_now_responses <- data$health_now_2 - data$health_now_1 # We use the 2nd vs the 1st measure to estimate effect of endotoxin
boxplot(health_now_responses ~ group, data = data, frame.plot = F, main = "health_now difference", ylab = " pg/ml")
ttest_health_now <- t.test(health_now_responses ~ group, data = data)
es_health_now <- tes(t = ttest_health_now$statistic, n.1 = length(data$health_now_responses[data$group == "Endotoxin"]), n.2 = length(data$health_now_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_health_now$d, " [", es_health_now$l.d, ", ", es_health_now$u.d, "]", sep = ""), bty = "n")


# 14. Headache
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

data$headache_responses <- rowMeans(cbind(data$headache_2, data$headache_3), na.rm = T) - data$headache_0 # We use the mean of the 3rd and 4th measures to estimate effect of endotoxin
boxplot(headache_responses ~ group, data = data, frame.plot = F, main = "headache difference", ylab = " pg/ml")
ttest_headache <- t.test(headache_responses ~ group, data = data)
es_headache <- tes(t = ttest_headache$statistic, n.1 = length(data$headache_responses[data$group == "Endotoxin"]), n.2 = length(data$headache_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_headache$d, " [", es_headache$l.d, ", ", es_headache$u.d, "]", sep = ""), bty = "n")


# 15. Back pain
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

data$back_pain_responses <- rowMeans(cbind(data$back_pain_2, data$back_pain_3), na.rm = T) - data$back_pain_0 # We use the mean of the 3rd and 4th measures to estimate effect of endotoxin
boxplot(back_pain_responses ~ group, data = data, frame.plot = F, main = "back_pain difference", ylab = " pg/ml")
ttest_back_pain <- t.test(back_pain_responses ~ group, data = data)
es_back_pain <- tes(t = ttest_back_pain$statistic, n.1 = length(data$back_pain_responses[data$group == "Endotoxin"]), n.2 = length(data$back_pain_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_back_pain$d, " [", es_back_pain$l.d, ", ", es_back_pain$u.d, "]", sep = ""), bty = "n")


# 16. Nausea
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

data$nausea_responses <- rowMeans(cbind(data$nausea_2, data$nausea_3), na.rm = T) - data$nausea_0 # We use the mean of the 3rd and 4th measures to estimate effect of endotoxin
boxplot(nausea_responses ~ group, data = data, frame.plot = F, main = "nausea difference", ylab = " pg/ml")
ttest_nausea <- t.test(nausea_responses ~ group, data = data)
es_nausea <- tes(t = ttest_nausea$statistic, n.1 = length(data$nausea_responses[data$group == "Endotoxin"]), n.2 = length(data$nausea_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_nausea$d, " [", es_nausea$l.d, ", ", es_nausea$u.d, "]", sep = ""), bty = "n")


# 17. Thumb pressure
pain_threshold_data <- data[, c("subject", "group", "pain_threshold_baseline", "pain_threshold_diff")]
pain_threshold_data2 <- pain_threshold_data
pain_threshold_data2$pain_threshold_1 <- pain_threshold_data2$pain_threshold_baseline
pain_threshold_data2$pain_threshold_2 <- pain_threshold_data2$pain_threshold_1 + pain_threshold_data2$pain_threshold_diff
pain_threshold_data2 <- pain_threshold_data2[, c("subject", "group", "pain_threshold_1", "pain_threshold_2")]
pain_threshold_data2 <- reshape(pain_threshold_data2, direction = 'long', varying = 3:4, v.names = "pain_threshold", timevar = "time", sep = "_")
pain_threshold_data2 <- pain_threshold_data2[, -5]
pain_threshold_agg <- aggregate(pain_threshold_data2, by=list(as.logical(as.integer(pain_threshold_data2$group)-1), as.numeric(pain_threshold_data2$time)), FUN=mean, na.rm=TRUE)

hist(pain_threshold_data2$pain_threshold)

pdf("Fig_pain_threshold.pdf")
par(mar=c(6, 5, 1, 1))
plot(pain_threshold ~ time, data = pain_threshold_data2, pch = "", xlab = "Time of measurement", ylab = "Pain threshold", main = "", bty = "n", cex.axis = 2, cex.lab = 2, mgp = c(3.5, 1.5, 0))
for (i in unique(pain_threshold_data2$subject[pain_threshold_data2$group == "Endotoxin"])){
  lines(pain_threshold ~ time, data = pain_threshold_data2[pain_threshold_data2$subject == i, ], col = "pink")
}
for (i in unique(pain_threshold_data2$subject[pain_threshold_data2$group == "Placebo"])){
  lines(pain_threshold ~ time, data = pain_threshold_data2[pain_threshold_data2$subject == i, ], col = "lightblue")
}
lines(pain_threshold ~ time, data = pain_threshold_agg[pain_threshold_agg$Group.1 == TRUE, ], col = "blue", lwd = 3)
lines(pain_threshold ~ time, data = pain_threshold_agg[pain_threshold_agg$Group.1 == FALSE, ], col = "red", lwd = 3)
legend("topleft", legend = c("Endotoxin", "Placebo"), col = c("red", "blue"), lty = 1, lwd = 3, bty = "n", cex = 2)
dev.off()

data$pain_threshold_responses <- data$pain_threshold_diff
boxplot(pain_threshold_responses ~ group, data = data, frame.plot = F, main = "pain_threshold difference", ylab = " pg/ml")
ttest_pain_threshold <- t.test(pain_threshold_responses ~ group, data = data)
es_pain_threshold <- tes(t = ttest_pain_threshold$statistic, n.1 = length(data$pain_threshold_responses[data$group == "Endotoxin"]), n.2 = length(data$pain_threshold_responses[data$group == "Placebo"]))
legend("topright", legend = paste("d = ", es_pain_threshold$d, " [", es_pain_threshold$l.d, ", ", es_pain_threshold$u.d, "]", sep = ""), bty = "n")


# WRITE DATA IN LONG FORMAT
PhysiologicalDataLong <- merge(HR_data, temp_data, by = c("subject", "time"))
PhysiologicalDataLong <- merge(PhysiologicalDataLong, BP_systolic_data, by = c("subject", "time"))
PhysiologicalDataLong <- PhysiologicalDataLong[, c("subject", "group", "time", "HR", "temp", "BP_systolic")]
PhysiologicalDataLong <- PhysiologicalDataLong[order(PhysiologicalDataLong$subject), ]
PhysiologicalDataLong$time <- PhysiologicalDataLong$time -1
write.csv(PhysiologicalDataLong, file = "C:/Users/Gustav Nilsonne/Box Sync/Gustavs_arbete/Pek/Endotoxin resting state/Endotoxin_RestingState_2015/PhysiologicalDataLong.csv", row.names=FALSE)

SelfRatedDataLong <- merge(headache_data, back_pain_data, by = c("subject", "time"))
SelfRatedDataLong <- merge(SelfRatedDataLong, nausea_data, by = c("subject", "time"))
SelfRatedDataLong <- SelfRatedDataLong[, c("subject", "group", "time", "headache", "back_pain", "nausea")]
SelfRatedDataLong <- SelfRatedDataLong[order(SelfRatedDataLong$subject), ]
SelfRatedDataLong$time <- SelfRatedDataLong$time -1
write.csv(SelfRatedDataLong, file = "C:/Users/Gustav Nilsonne/Box Sync/Gustavs_arbete/Pek/Endotoxin resting state/Endotoxin_RestingState_2015/SelfRatedDataLong.csv", row.names=FALSE)


# TABULATE EFFECTS OF ENDOTOXIN
effects_table <- data.frame(variable = c("Heart rate", "Temperature", "Systolic blood pressure", "log IL-6", "log IL-8", "log IL-10", "log IL-13", "log TNF",
                                         "SQ total", "SQ pain", "SQ fatigue", "SQ affect", "Self-rated health now", "Headache", "Back pain", "Nausea", "Pain threshold"))
effects_table$timepoints_for_comparison <- c("3 and 4 vs 1", "3 and 4 vs 1", "3 vs 1", "2 and 3 vs 1", "2 and 3 vs 1", "2 and 3 vs 1", "2 and 3 vs 1", "2 and 3 vs 1",
                              "2 vs 1", "2 vs 1", "2 vs 1", "2 vs 1", "2 vs 1", "3 and 4 vs 1", "3 and 4 vs 1", "3 and 4 vs 1", "2 vs 1")
effects_table$effect_orig_scale <- c(ttest_HR$estimate[1] - ttest_HR$estimate[2], 
                                     ttest_temp$estimate[1] - ttest_temp$estimate[2],
                                     ttest_BP_systolic$estimate[1] - ttest_BP_systolic$estimate[2],
                                     ttest_IL6$estimate[1] - ttest_IL6$estimate[2],
                                     ttest_IL8$estimate[1] - ttest_IL8$estimate[2],
                                     ttest_IL10$estimate[1] - ttest_IL10$estimate[2],
                                     ttest_IL13$estimate[1] - ttest_IL13$estimate[2],
                                     ttest_TNF$estimate[1] - ttest_TNF$estimate[2],
                                     ttest_sq_total$estimate[1] - ttest_sq_total$estimate[2],
                                     ttest_sq_pain$estimate[1] - ttest_sq_pain$estimate[2],
                                     ttest_sq_fatigue$estimate[1] - ttest_sq_fatigue$estimate[2],
                                     ttest_sq_affect$estimate[1] - ttest_sq_affect$estimate[2],
                                     ttest_health_now$estimate[1] - ttest_health_now$estimate[2],
                                     ttest_headache$estimate[1] - ttest_headache$estimate[2],
                                     ttest_back_pain$estimate[1] - ttest_back_pain$estimate[2],
                                     ttest_nausea$estimate[1] - ttest_nausea$estimate[2],
                                     ttest_pain_threshold$estimate[1] - ttest_pain_threshold$estimate[2])
effects_table$CI_orig_scale_lower <- c(ttest_HR$conf.int[1], 
                                       ttest_temp$conf.int[1],
                                 ttest_BP_systolic$conf.int[1],
                                 ttest_IL6$conf.int[1],
                                 ttest_IL8$conf.int[1],
                                 ttest_IL10$conf.int[1],
                                 ttest_IL13$conf.int[1],
                                 ttest_TNF$conf.int[1],
                                 ttest_sq_total$conf.int[1],
                                 ttest_sq_pain$conf.int[1],
                                 ttest_sq_fatigue$conf.int[1],
                                 ttest_sq_affect$conf.int[1],
                                 ttest_health_now$conf.int[1],
                                 ttest_headache$conf.int[1],
                                 ttest_back_pain$conf.int[1],
                                 ttest_nausea$conf.int[1],
                                 ttest_pain_threshold$conf.int[1])
effects_table$CI_orig_scale_upper <- c(ttest_HR$conf.int[2], 
                                       ttest_temp$conf.int[2],
                                       ttest_BP_systolic$conf.int[2],
                                       ttest_IL6$conf.int[2],
                                       ttest_IL8$conf.int[2],
                                       ttest_IL10$conf.int[2],
                                       ttest_IL13$conf.int[2],
                                       ttest_TNF$conf.int[2],
                                       ttest_sq_total$conf.int[2],
                                       ttest_sq_pain$conf.int[2],
                                       ttest_sq_fatigue$conf.int[2],
                                       ttest_sq_affect$conf.int[2],
                                       ttest_health_now$conf.int[2],
                                       ttest_headache$conf.int[2],
                                       ttest_back_pain$conf.int[2],
                                       ttest_nausea$conf.int[2],
                                       ttest_pain_threshold$conf.int[2])
effects_table$d <- c(es_HR$d, es_temp$d, es_BP_systolic$d, es_IL6$d, es_IL8$d, es_IL10$d, es_IL13$d, es_TNF$d,
                     es_sq_total$d, es_sq_pain$d, es_sq_fatigue$d, es_sq_affect$d, es_health_now$d, es_headache$d, es_back_pain$d, es_nausea$d, es_pain_threshold$d)
effects_table$CI_d_lower <- c(es_HR$l.d, es_temp$l.d, es_BP_systolic$l.d, es_IL6$l.d, es_IL8$l.d, es_IL10$l.d, es_IL13$l.d, es_TNF$l.d,
                              es_sq_total$l.d, es_sq_pain$l.d, es_sq_fatigue$l.d, es_sq_affect$l.d, es_health_now$l.d, es_headache$l.d, es_back_pain$l.d, es_nausea$l.d, es_pain_threshold$l.d)
effects_table$CI_d_upper <- c(es_HR$u.d, es_temp$u.d, es_BP_systolic$u.d, es_IL6$u.d, es_IL8$u.d, es_IL10$u.d, es_IL13$u.d, es_TNF$u.d,
                              es_sq_total$u.d, es_sq_pain$u.d, es_sq_fatigue$u.d, es_sq_affect$u.d, es_health_now$u.d, es_headache$u.d, es_back_pain$u.d, es_nausea$u.d, es_pain_threshold$u.d)
effects_table$t <- c(ttest_HR$statistic, ttest_temp$statistic,
                     ttest_BP_systolic$statistic,
                     ttest_IL6$statistic,
                     ttest_IL8$statistic,
                     ttest_IL10$statistic,
                     ttest_IL13$statistic,
                     ttest_TNF$statistic,
                     ttest_sq_total$statistic,
                     ttest_sq_pain$statistic,
                     ttest_sq_fatigue$statistic,
                     ttest_sq_affect$statistic,
                     ttest_health_now$statistic,
                     ttest_headache$statistic,
                     ttest_back_pain$statistic,
                     ttest_nausea$statistic,
                     ttest_pain_threshold$statistic)
effects_table$df <- c(ttest_HR$parameter, ttest_temp$parameter,
                      ttest_BP_systolic$parameter,
                      ttest_IL6$parameter,
                      ttest_IL8$parameter,
                      ttest_IL10$parameter,
                      ttest_IL13$parameter,
                      ttest_TNF$parameter,
                      ttest_sq_total$parameter,
                      ttest_sq_pain$parameter,
                      ttest_sq_fatigue$parameter,
                      ttest_sq_affect$parameter,
                      ttest_health_now$parameter,
                      ttest_headache$parameter,
                      ttest_back_pain$parameter,
                      ttest_nausea$parameter,
                      ttest_pain_threshold$parameter)
effects_table$p <- c(ttest_HR$p.value, ttest_temp$p.value,
                     ttest_BP_systolic$p.value,
                     ttest_IL6$p.value,
                     ttest_IL8$p.value,
                     ttest_IL10$p.value,
                     ttest_IL13$p.value,
                     ttest_TNF$p.value,
                     ttest_sq_total$p.value,
                     ttest_sq_pain$p.value,
                     ttest_sq_fatigue$p.value,
                     ttest_sq_affect$p.value,
                     ttest_health_now$p.value,
                     ttest_headache$p.value,
                     ttest_back_pain$p.value,
                     ttest_nausea$p.value,
                     ttest_pain_threshold$p.value)
write.csv(effects_table, file = "C:/Users/Gustav Nilsonne/Box Sync/Gustavs_arbete/Pek/Endotoxin resting state/Endotoxin_RestingState_2015/EffectsTable.csv", row.names=FALSE)


# ANALYSE AND PLOT BIVARIATE RELATIONSHIPS TO BRAIN CONNECTIVITY
pdf("Fig_IL6corr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ IL6_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "log IL-6, pg/ml", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_IL6 <- lm(midcingulate_insula_conn ~ IL6_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$IL6_responses, na.rm = T), max(data$IL6_responses, na.rm = T))
Y <- predict(lm_IL6, newdata=data.frame(IL6_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_IL6$coefficients[2], 2)
#coefLower <- round(confint(lm_IL6)[2, 1], 2)
#coefUpper <- round(confint(lm_IL6)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_IL6)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_IL6)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_IL8corr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ IL8_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "log IL-8, pg/ml", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_IL8 <- lm(midcingulate_insula_conn ~ IL8_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$IL8_responses, na.rm = T), max(data$IL8_responses, na.rm = T))
Y <- predict(lm_IL8, newdata=data.frame(IL8_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_IL8$coefficients[2], 2)
#coefLower <- round(confint(lm_IL8)[2, 1], 2)
#coefUpper <- round(confint(lm_IL8)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_IL8)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_IL8)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_IL10corr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ IL10_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "log IL-10, pg/ml", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_IL10 <- lm(midcingulate_insula_conn ~ IL10_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$IL10_responses, na.rm = T), max(data$IL10_responses, na.rm = T))
Y <- predict(lm_IL10, newdata=data.frame(IL10_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_IL10$coefficients[2], 2)
#coefLower <- round(confint(lm_IL10)[2, 1], 2)
#coefUpper <- round(confint(lm_IL10)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_IL10)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_IL10)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_IL13corr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ IL13_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "log IL-13, pg/ml", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_IL13 <- lm(midcingulate_insula_conn ~ IL13_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$IL13_responses, na.rm = T), max(data$IL13_responses, na.rm = T))
Y <- predict(lm_IL13, newdata=data.frame(IL13_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_IL13$coefficients[2], 2)
#coefLower <- round(confint(lm_IL13)[2, 1], 2)
#coefUpper <- round(confint(lm_IL13)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_IL13)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_IL13)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_TNFcorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ TNF_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "log TNF, pg/ml", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_TNF <- lm(midcingulate_insula_conn ~ TNF_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$TNF_responses, na.rm = T), max(data$TNF_responses, na.rm = T))
Y <- predict(lm_TNF, newdata=data.frame(TNF_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_TNF$coefficients[2], 2)
#coefLower <- round(confint(lm_TNF)[2, 1], 2)
#coefUpper <- round(confint(lm_TNF)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_TNF)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_TNF)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_sq_totalcorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ sq_total_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "Score", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_sq_total <- lm(midcingulate_insula_conn ~ sq_total_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$sq_total_responses, na.rm = T), max(data$sq_total_responses, na.rm = T))
Y <- predict(lm_sq_total, newdata=data.frame(sq_total_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_sq_total$coefficients[2], 2)
#coefLower <- round(confint(lm_sq_total)[2, 1], 2)
#coefUpper <- round(confint(lm_sq_total)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_sq_total)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_sq_total)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_sq_paincorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ sq_pain_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "Score", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_sq_pain <- lm(midcingulate_insula_conn ~ sq_pain_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$sq_pain_responses, na.rm = T), max(data$sq_pain_responses, na.rm = T))
Y <- predict(lm_sq_pain, newdata=data.frame(sq_pain_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_sq_pain$coefficients[2], 2)
#coefLower <- round(confint(lm_sq_pain)[2, 1], 2)
#coefUpper <- round(confint(lm_sq_pain)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_sq_pain)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_sq_pain)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_sq_fatiguecorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ sq_fatigue_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "Score", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_sq_fatigue <- lm(midcingulate_insula_conn ~ sq_fatigue_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$sq_fatigue_responses, na.rm = T), max(data$sq_fatigue_responses, na.rm = T))
Y <- predict(lm_sq_fatigue, newdata=data.frame(sq_fatigue_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_sq_fatigue$coefficients[2], 2)
#coefLower <- round(confint(lm_sq_fatigue)[2, 1], 2)
#coefUpper <- round(confint(lm_sq_fatigue)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_sq_fatigue)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_sq_fatigue)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_sq_affectcorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ sq_affect_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "Score", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_sq_affect <- lm(midcingulate_insula_conn ~ sq_affect_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$sq_affect_responses, na.rm = T), max(data$sq_affect_responses, na.rm = T))
Y <- predict(lm_sq_affect, newdata=data.frame(sq_affect_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_sq_affect$coefficients[2], 2)
#coefLower <- round(confint(lm_sq_affect)[2, 1], 2)
#coefUpper <- round(confint(lm_sq_affect)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_sq_affect)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_sq_affect)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_health_nowcorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ health_now_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "Score", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_health_now <- lm(midcingulate_insula_conn ~ health_now_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$health_now_responses, na.rm = T), max(data$health_now_responses, na.rm = T))
Y <- predict(lm_health_now, newdata=data.frame(health_now_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_health_now$coefficients[2], 2)
#coefLower <- round(confint(lm_health_now)[2, 1], 2)
#coefUpper <- round(confint(lm_health_now)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_health_now)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_health_now)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_headachecorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ headache_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "Score", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_headache <- lm(midcingulate_insula_conn ~ headache_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$headache_responses, na.rm = T), max(data$headache_responses, na.rm = T))
Y <- predict(lm_headache, newdata=data.frame(headache_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_headache$coefficients[2], 2)
#coefLower <- round(confint(lm_headache)[2, 1], 2)
#coefUpper <- round(confint(lm_headache)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_headache)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_headache)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_back_paincorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ back_pain_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "Score", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_back_pain <- lm(midcingulate_insula_conn ~ back_pain_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$back_pain_responses, na.rm = T), max(data$back_pain_responses, na.rm = T))
Y <- predict(lm_back_pain, newdata=data.frame(back_pain_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_back_pain$coefficients[2], 2)
#coefLower <- round(confint(lm_back_pain)[2, 1], 2)
#coefUpper <- round(confint(lm_back_pain)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_back_pain)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_back_pain)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_nauseacorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ nausea_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "Score", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_nausea <- lm(midcingulate_insula_conn ~ nausea_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$nausea_responses, na.rm = T), max(data$nausea_responses, na.rm = T))
Y <- predict(lm_nausea, newdata=data.frame(nausea_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_nausea$coefficients[2], 2)
#coefLower <- round(confint(lm_nausea)[2, 1], 2)
#coefUpper <- round(confint(lm_nausea)[2, 2], 2)
legend("bottomright", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_nausea)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_nausea)$coefficients[2, 4], 2))))))
dev.off()

pdf("Fig_pain_thresholdcorr.pdf")
par(mar=c(6, 5, 0, 1))
plot(midcingulate_insula_conn ~ pain_threshold_responses, data = data[data$group == "Endotoxin", ], bty = "n", xlab = "Score", ylab = "Connectivity (unit)", cex.axis = 2.5, cex.lab = 2.5, pch = 16, cex = 2, col = "red", yaxt = "n")
axis(2, at = c(-0.3, -0.1, 0.1, 0.3), cex.axis = 2.5)
lm_pain_threshold <- lm(midcingulate_insula_conn ~ pain_threshold_responses, data = data[data$group == "Endotoxin", ])
X <- c(min(data$pain_threshold_responses, na.rm = T), max(data$pain_threshold_responses, na.rm = T))
Y <- predict(lm_pain_threshold, newdata=data.frame(pain_threshold_responses=X))
lines(x=X, y=Y, col = "red", lwd = 2)
coef <- round(lm_pain_threshold$coefficients[2], 2)
#coefLower <- round(confint(lm_pain_threshold)[2, 1], 2)
#coefUpper <- round(confint(lm_pain_threshold)[2, 2], 2)
legend("bottomleft", bty = "n", cex = 2.5, inset = c(-0.07, -0.02), legend = c(  
  as.expression(bquote(beta == .(coef))),
  as.expression(bquote(r[2] == .(round(summary(lm_pain_threshold)$r.squared, 2)))),
  as.expression(bquote(p == .(round(summary(lm_pain_threshold)$coefficients[2, 4], 2))))))
dev.off()