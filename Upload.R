#Read data
rad <- read.csv("D:/Articles/Article12 pathomics/radiomics/radiomics.csv", sep = ",", header = TRUE)
#Standardized and reshape
rad1 <- scale(rad[,6:2291], center = TRUE, scale = TRUE)
rad2 <- cbind(rad[,1:5], rad1)
#LASSO/ridge/elastic net feature selection
library(survival)
library(glmnet)
surv <- Surv(rad$RFS, rad$RFS.status)
cv.fit <- cv.glmnet(rad1, surv, family = "cox", nfolds = 10, type.measure = "C", alpha = 1) #alpha = 1, LASSO； alpha = 0~1, elastic net; alpha = 0, ridge
cv.fit$lambda.min
plot(cv.fit)
coefficient <- coef(cv.fit, s= cv.fit$lambda.min)
selected_index <- which(as.numeric(coefficient) != 0)
selected_features <- names(coefficient[selected_index,])
fit <- glmnet(rad1, surv, family = "cox", alpha = 1)
plot(fit,xvar = "lambda", lwd=1.5, xlim=c(-5,-2), ylim=c(-2,2))
+abline(v=log(0.07399609), col="orange", lwd=2)
+text(log(cv.fit$lambda.min), 1.4, paste("lambda.min=",round(cv.fit$lambda.min,4),"\n", sep="","feature number = 3"),
      col="black", cex=1, pos=4)
#random forest feature selection
library(randomForestSRC)
rfsrcfit <- rfsrc(Surv(RFS, RFS.status)~., data = rad3, ntree = 100, nsplit = 5, importance = TRUE, tree.err = TRUE)
plot(get.tree(rfsrcfit,1))
tune.nodesize(Surv(RFS, RFS.status)~., data = rad3)
rfsrcfit <- rfsrc(Surv(RFS, RFS.status)~., data = rad3, ntree = 100, nsplit = 45, importance = TRUE, tree.err = TRUE)
rfvs <- var.select.rfsrc(rfsrcfit, method = "vh", conservative = "high", nodesize = 45, ntree = 100, K = 10)
#Cox regression
library(rms)
library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2)
library(survival)
library(coin)
library(survminer)
library(Rcpp)
ddist <- datadist(training.set)
options(datadist = ddist)
nomocoxos <- cph(Surv(Survival.months,eventos)~Age+Reg.LN.positive+Tumor.sizes+Mitosis,x=T,y=T,data=training.set,surv=T)
surv.cox <- Survival(nomocoxos)
nom.cox <- nomogram(nomocoxos,fun=list(function(x) surv.cox(36, x),function(x) surv.cox(60, x)),funlabel=c("3-Year Sur. Prob.","5-Year Sur. Prob."),lp=F,fun.at=c('0.95','0.9','0.8','0.7','0.5','0.4','0.3','0.2','0.1'))
plot(nom.cox)
pre1 <- predict(nomocoxdss1,newdata = validation.set)
f1 <- cph(Surv(Survival.months, eventdss)~pre1, x=T, y=T, surv=T, data=validation.set, time.inc=36)
validate(f1)
rcorrcens(Surv(Survival.months, eventdss)~pre1,data = validation.set)
#c-index 1
cox <- coxph(Surv(RFS, RFS.status)~, data = )
summary(cox)
#test c-index
c_index <- function(data,indices){
      dat <- rad[indices,]
      vames<-c("boxsigmaimage_glrlm_LongRunHighGrayLevelEmphasis", "boxsigmaimage_glrlm_RunLengthNonUniformity", "wavelet_glszm_wavelet.LHH.GrayLevelNonUniformity")
      FML <- as.formula(paste('Surv(RFS, RFS.status)~',paste(vames, collapse = "+")))
      fit<- coxph(FML,data =dat )
      pr1<-predict(fit,newdata=test.set)
      Cindex=rcorrcens(Surv(RFS, RFS.status) ~ pr1, data = test.set)[1]
      Cindex=1-Cindex
      Cindex 
  }
library(boot)
bootci <- boot(data=test.set, statistic=c_index, R=500)
print(bootci)
plot(bootci)

#stepwise cox
for (direction in c("both", "backward")){
  fit <- step(coxph(Surv(RFS, RFS.status)~pfeature_83+pfeature_431+pfeature_776+pfeature_1183+pfeature_1576+pfeature_69+pfeature_382+pfeature_464+pfeature_741+pfeature_934+pfeature_1298+pfeature_1726, pat), direction = direction)
}
#cutoff
library(survminer)
library(nomogramFormula)
results <- formula_rd(nomogram = nom.cox)
rad$points <- points_cal(formula = results$formula,rd=rad)
test.set$points <- points_cal(formula = results$formula,rd=test.set)
cutoff <- surv_cutpoint(rad, time = "RFS", event = "RFS.status", variables = "points")
summary(cutoff)
plot(cutoff, "points", palette = "npg")
test.set$radscore <- ifelse(test.set$points >= 119.7455, "High", "Low")
table(test.set$radscore)
KMtrain <- survfit(Surv(RFS, RFS.status)~radscore,data=rad)
B <- ggsurvplot(KMtrain, data=rad, ggtheme = theme_bw(), pval = T,break.time.by=12, censor=TRUE, xlim=c(0,60),risk.table = "abs_pct",risk.table.fontsize=2.5,risk.tables.height=0.3,legend.labs=c("High","Low"))+xlab("time")+ylab("rate")
B
#CV validation
folds <- createMultiFolds(y=rad$RFS.status, k=10, times = 200)
cindex_value <- as.numeric()
for (i in 1:2000){
      train <- pat[folds[[i]],]
      test <- pat[-folds[[i]],]
      model <- coxph(Surv(RFS, RFS.status)~pfeature_69+pfeature_83+pfeature_382+pfeature_464+pfeature_741+pfeature_934+pfeature_1298+pfeature_1726, data=train)
      pr1 <- predict(model, newdata = test)
      Cindex <- rcorrcens(Surv(RFS, RFS.status) ~ pr1, data = test)[1]
      Cindex = 1-Cindex
      cindex_value<- append(cindex_value, as.numeric(Cindex))
  }
summary(cindex_value)

#correlation
cor <- round(cor(correlation), 3)
pcor <- cor_pmat(correlation)
ggcorrplot(cor,hc.order = T, ggtheme = ggplot2::theme_void(base_size = 15), colors = c("CornflowerBlue","white","Salmon"), lab = T,lab_size = 5, tl.cex = 15,p.mat = pcor, sig.level = 0.01,pch = 4,pch.cex = 10)
library(corrmorant)
library(ggplot2)
p1 <- ggcorrm(data = , corr_method = c('pearson'))+lotri(geom_point(alpha = 0.3))+lotri(geom_smooth(method = 'lm'))+utri_corrtext(corr_size = FALSE)+dia_names(y_pos = 0.15, size = 3.5)+dia_density(lower = 0.3, upper = 0.98, alpha = 0.3)
#feature map
library(pheatmap)
mapt <- read.csv("D:/Articles/Article12 pathomics/XAI/test.csv", header = T, sep = ",")
mapt$x <- mapt$X/224
mapt$y <- mapt$Y/224
matrix1 <- matrix(0, 500, 500)
for (i in 1:98648){
  x = mapt[i, 12]
  y = mapt[i, 13]
  matrix2[y, x] <- mapt[i, 4]
}
pheatmap(matrix1, color = colorRampPalette(c("white", "orange"))(10), cluster_rows = FALSE, cluster_cols = FALSE)
#bar plot
g <- ggplot(total, aes(Model, mu))+geom_col(aes(fill = Model), width = 0.5)+geom_errorbar(aes(ymax = mu + se, ymin = mu - se), width = 0.3, cex = 1)



library(rms)
library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2)
library(survival)
library(coin)
library(survminer)
library(Rcpp)


#KM
gistKM <- survfit(Surv(Survival.months,eventos)~,data=)
B <- ggsurvplot(gistKM, data=gist, ggtheme = theme_bw(), pval = T,break.time.by=12, censor=TRUE, xlim=c(0,120),risk.table = "abs_pct",risk.table.fontsize=2.5,risk.tables.height=0.3,legend.labs=c(""))+xlab("time")+ylab("rate")
b <- coxph(Surv(gist$Survival.months,gist$eventa)~gist$Age,data=gist)


#nomogram
ddist <- datadist(training.set)
options(datadist = ddist)
nomocoxos <- cph(Surv(Survival.months,eventos)~Age+Reg.LN.positive+Tumor.sizes+Mitosis,x=T,y=T,data=training.set,surv=T)
surv.cox <- Survival(nomocoxos)
nom.cox <- nomogram(nomocoxos,fun=list(function(x) surv.cox(36, x),function(x) surv.cox(60, x)),funlabel=c("3-Year Sur. Prob.","5-Year Sur. Prob."),lp=F,fun.at=c('0.95','0.9','0.8','0.7','0.5','0.4','0.3','0.2','0.1'))
plot(nom.cox)


#DCA
library(ggDCA)
dca(nomocoxos,times=)

#points calculation
library(nomogramFormula)
results <- formula_rd(nomogram = nom.cox)
training.set$pointsos <- points_cal(formula = results$formula,rd=training.set)

#RF Feature selection
library(randomForestSRC)
rfmfsvh <- var.select.rfsrc(Surv(RFS, RFS.status) ~ ., data = RF1t, tree.err = TRUE, importance = TRUE, method = "vh", conservative = "high", nrep = 100, K = 5, nstep = 1)


#XGBoost & SHAP
library(xgboost)
library(SHAPforxgboost)
xgbmodel <- xgboost(data = tumormatrix1, label = train_label, objective='binary:logistic', nrounds = 100)
shap_values <- shap.values(xgb_model = mod, X_train = dataX)
shap_values$mean_shap_score
shap_long <- shap.prep(xgb_model = xgbmodel, X_train = tumormatrix1, top_n = 20)
shap.plot.summary(shap_long)

#Sangi plot
library(ggplot2)
library(ggalluvial)
sangi <- read.csv("D:/Articles/Article13 muscle/ml/sangi.csv", header = T, sep = ",")
sangidata <- to_lodes_form(sangi[,1:ncol(sangi)],axes = 1:ncol(sangi),id = "value")
ggplot(sangidata, aes(x = x, fill=stratum, label=stratum, stratum = stratum, alluvium  = value))+geom_flow(width = 0.3, curve_type = "sine",alpha = 0.5,color = 'white',size = 0.1)+geom_stratum(width = 0.28)+geom_text(stat = 'stratum', size = 2, color = 'black')+theme_void()+theme(legend.position = 'none')

#Cluster
library(factoextra)

#bar
barplot <- ggplot(datalabel, aes(x = musclescore, y = ALB, fill = musclescore)) + geom_boxplot(outlier.shape = NA) + theme_classic()
barplot + stat_compare_means(label = "p.format", label.x = 1.5)

#Bioenvironment pipline
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(factoextra)
library(pheatmap)
b <- read.csv("D:/Articles/Article17 生境/Data/Z0952025/tumor feature.csv", sep = ",", header = F)
cname <- c("patch")
for (i in 1:2048){
   cname <- append(cname, paste("pfeature_", i, sep = ""))
}
colnames(b) <- cname
b1 <- b[, -1]

tsne <- Rtsne(b1)
b2 <- tsne$Y
fviz_nbclust(b2, kmeans, method = "wss")
fviz_nbclust(b2, kmeans, method = "silhouette")
b2 <- as.data.frame(b2)
colnames(b2) <- c("TSNE1", "TSNE2")
cluster <- kmeans(b2, 4)
fviz_cluster(cluster, data = b2, labelsize = 0)
b$cluster <- cluster$cluster

b3 <- cbind(b$patch, b$cluster)
b3 <- as.data.frame(b3)
colnames(b3) <- c("patch", "cluster")
write.csv(b3, "D:/Articles/Article17 生境/Data/Z0952025/map.csv")

mapt <- read.csv("D:/Articles/Article17 生境/Data/Z0952025/map.csv", sep = ",", header = T)
mapt$x <- mapt$X/512
mapt$y <- mapt$Y/512
range(mapt$x)
range(mapt$y)
matrix1 <- matrix(0, 250, 250)
for (i in 1:15978){
  x = mapt[i, 5]
  y = mapt[i, 6]
  matrix1[y, x] <- mapt[i, 2]
}
pheatmap(matrix1, color = c("white", "red", "green", "blue", "purple"), cluster_rows = FALSE, cluster_cols = FALSE)

bcluster1 <- subset(b, b$cluster == "1")
bcluster2 <- subset(b, b$cluster == "2")
bcluster3 <- subset(b, b$cluster == "3")
bcluster4 <- subset(b, b$cluster == "4")
cluster1 <- c()
cluster2 <- c()
cluster3 <- c()
cluster4 <- c()
for (i in 2:2049){
  cluster1 <- append(cluster1, mean(bcluster1[, i]))
}
clusterfeature <- rbind(cluster1, cluster2)
clusterfeature <- rbind(clusterfeature, cluster3)
clusterfeature <- rbind(clusterfeature, cluster4)

#patch predict
#cluster1
results <- formula_rd(nomogram = nom.cox1)
patch1$c1m1 <- points_cal(formula = results$formula,rd=patch1)
results <- formula_rd(nomogram = nom.cox2)
patch1$c1m2 <- points_cal(formula = results$formula,rd=patch1)
results <- formula_rd(nomogram = nom.cox3)
patch1$c1m3 <- points_cal(formula = results$formula,rd=patch1)
results <- formula_rd(nomogram = nom.cox4)
patch1$c1m4 <- points_cal(formula = results$formula,rd=patch1)
results <- formula_rd(nomogram = nom.cox5)
patch1$c1m5 <- points_cal(formula = results$formula,rd=patch1)
results <- formula_rd(nomogram = nom.cox6)
patch1$c1m6 <- points_cal(formula = results$formula,rd=patch1)
results <- formula_rd(nomogram = nom.cox7)
patch1$c1m7 <- points_cal(formula = results$formula,rd=patch1)
results <- formula_rd(nomogram = nom.cox8)
patch1$c1m8 <- points_cal(formula = results$formula,rd=patch1)
results <- formula_rd(nomogram = nom.cox9)
patch1$c1m9 <- points_cal(formula = results$formula,rd=patch1)
results <- formula_rd(nomogram = nom.cox10)
patch1$c1m10 <- points_cal(formula = results$formula,rd=patch1)

patch1$c1m1l <- ifelse(patch1$c1m1>149.83, 1, 0)
patch1$c1m2l <- ifelse(patch1$c1m2>172.02, 1, 0)
patch1$c1m3l <- ifelse(patch1$c1m3>172.58, 1, 0)
patch1$c1m4l <- ifelse(patch1$c1m4>144, 1, 0)
patch1$c1m5l <- ifelse(patch1$c1m5>154.8, 1, 0)
patch1$c1m6l <- ifelse(patch1$c1m6>179.96, 1, 0)
patch1$c1m7l <- ifelse(patch1$c1m7>170.74, 1, 0)
patch1$c1m8l <- ifelse(patch1$c1m8>190.27, 1, 0)
patch1$c1m9l <- ifelse(patch1$c1m9>174.1, 1, 0)
patch1$c1m10l <- ifelse(patch1$c1m10>169.36, 1, 0)

table(patch1$c1m1l)+table(patch1$c1m2l)+table(patch1$c1m3l)+table(patch1$c1m4l)+table(patch1$c1m5l)+table(patch1$c1m6l)+table(patch1$c1m7l)+table(patch1$c1m8l)+table(patch1$c1m9l)+table(patch1$c1m10l)

#cluster2
results <- formula_rd(nomogram = nom.cox11)
patch2$c2m1 <- points_cal(formula = results$formula,rd=patch2)
results <- formula_rd(nomogram = nom.cox12)
patch2$c2m2 <- points_cal(formula = results$formula,rd=patch2)
results <- formula_rd(nomogram = nom.cox13)
patch2$c2m3 <- points_cal(formula = results$formula,rd=patch2)
results <- formula_rd(nomogram = nom.cox14)
patch2$c2m4 <- points_cal(formula = results$formula,rd=patch2)
results <- formula_rd(nomogram = nom.cox15)
patch2$c2m5 <- points_cal(formula = results$formula,rd=patch2)
results <- formula_rd(nomogram = nom.cox16)
patch2$c2m6 <- points_cal(formula = results$formula,rd=patch2)
results <- formula_rd(nomogram = nom.cox17)
patch2$c2m7 <- points_cal(formula = results$formula,rd=patch2)
results <- formula_rd(nomogram = nom.cox18)
patch2$c2m8 <- points_cal(formula = results$formula,rd=patch2)
results <- formula_rd(nomogram = nom.cox19)
patch2$c2m9 <- points_cal(formula = results$formula,rd=patch2)
results <- formula_rd(nomogram = nom.cox20)
patch2$c2m10 <- points_cal(formula = results$formula,rd=patch2)

patch2$c2m1l <- ifelse(patch2$c2m1>151.53, 1, 0)
patch2$c2m2l <- ifelse(patch2$c2m2>177.6, 1, 0)
patch2$c2m3l <- ifelse(patch2$c2m3>155.64, 1, 0)
patch2$c2m4l <- ifelse(patch2$c2m4>153.34, 1, 0)
patch2$c2m5l <- ifelse(patch2$c2m5>150.62, 1, 0)
patch2$c2m6l <- ifelse(patch2$c2m6>153.34, 1, 0)
patch2$c2m7l <- ifelse(patch2$c2m7>149.56, 1, 0)
patch2$c2m8l <- ifelse(patch2$c2m8>153.41, 1, 0)
patch2$c2m9l <- ifelse(patch2$c2m9>155.01, 1, 0)
patch2$c2m10l <- ifelse(patch2$c2m10>152.21, 1, 0)

table(patch2$c2m1l)+table(patch2$c2m2l)+table(patch2$c2m3l)+table(patch2$c2m4l)+table(patch2$c2m5l)+table(patch2$c2m6l)+table(patch2$c2m7l)+table(patch2$c2m8l)+table(patch2$c2m9l)+table(patch2$c2m10l)

#cluster3
results <- formula_rd(nomogram = nom.cox21)
patch3$c3m1 <- points_cal(formula = results$formula,rd=patch3)
results <- formula_rd(nomogram = nom.cox22)
patch3$c3m2 <- points_cal(formula = results$formula,rd=patch3)
results <- formula_rd(nomogram = nom.cox23)
patch3$c3m3 <- points_cal(formula = results$formula,rd=patch3)
results <- formula_rd(nomogram = nom.cox24)
patch3$c3m4 <- points_cal(formula = results$formula,rd=patch3)
results <- formula_rd(nomogram = nom.cox25)
patch3$c3m5 <- points_cal(formula = results$formula,rd=patch3)
results <- formula_rd(nomogram = nom.cox26)
patch3$c3m6 <- points_cal(formula = results$formula,rd=patch3)
results <- formula_rd(nomogram = nom.cox27)
patch3$c3m7 <- points_cal(formula = results$formula,rd=patch3)
results <- formula_rd(nomogram = nom.cox28)
patch3$c3m8 <- points_cal(formula = results$formula,rd=patch3)
results <- formula_rd(nomogram = nom.cox29)
patch3$c3m9 <- points_cal(formula = results$formula,rd=patch3)
results <- formula_rd(nomogram = nom.cox30)
patch3$c3m10 <- points_cal(formula = results$formula,rd=patch3)

patch3$c3m1l <- ifelse(patch3$c3m1>175.14, 1, 0)
patch3$c3m2l <- ifelse(patch3$c3m2>156.95, 1, 0)
patch3$c3m3l <- ifelse(patch3$c3m3>185.79, 1, 0)
patch3$c3m4l <- ifelse(patch3$c3m4>169.02, 1, 0)
patch3$c3m5l <- ifelse(patch3$c3m5>163.39, 1, 0)
patch3$c3m6l <- ifelse(patch3$c3m6>157.67, 1, 0)
patch3$c3m7l <- ifelse(patch3$c3m7>160.59, 1, 0)
patch3$c3m8l <- ifelse(patch3$c3m8>200.74, 1, 0)
patch3$c3m9l <- ifelse(patch3$c3m9>145.23, 1, 0)
patch3$c3m10l <- ifelse(patch3$c3m10>180.47, 1, 0)

table(patch3$c3m1l)+table(patch3$c3m2l)+table(patch3$c3m3l)+table(patch3$c3m4l)+table(patch3$c3m5l)+table(patch3$c3m6l)+table(patch3$c3m7l)+table(patch3$c3m8l)+table(patch3$c3m9l)+table(patch3$c3m10l)

#cluster4
results <- formula_rd(nomogram = nom.cox31)
patch4$c4m1 <- points_cal(formula = results$formula,rd=patch4)
results <- formula_rd(nomogram = nom.cox32)
patch4$c4m2 <- points_cal(formula = results$formula,rd=patch4)
results <- formula_rd(nomogram = nom.cox33)
patch4$c4m3 <- points_cal(formula = results$formula,rd=patch4)
results <- formula_rd(nomogram = nom.cox34)
patch4$c4m4 <- points_cal(formula = results$formula,rd=patch4)
results <- formula_rd(nomogram = nom.cox35)
patch4$c4m5 <- points_cal(formula = results$formula,rd=patch4)
results <- formula_rd(nomogram = nom.cox36)
patch4$c4m6 <- points_cal(formula = results$formula,rd=patch4)
results <- formula_rd(nomogram = nom.cox37)
patch4$c4m7 <- points_cal(formula = results$formula,rd=patch4)
results <- formula_rd(nomogram = nom.cox38)
patch4$c4m8 <- points_cal(formula = results$formula,rd=patch4)
results <- formula_rd(nomogram = nom.cox39)
patch4$c4m9 <- points_cal(formula = results$formula,rd=patch4)
results <- formula_rd(nomogram = nom.cox40)
patch4$c4m10 <- points_cal(formula = results$formula,rd=patch4)

patch4$c4m1l <- ifelse(patch4$c4m1>151.66, 1, 0)
patch4$c4m2l <- ifelse(patch4$c4m2>155.97, 1, 0)
patch4$c4m3l <- ifelse(patch4$c4m3>133.12, 1, 0)
patch4$c4m4l <- ifelse(patch4$c4m4>133.78, 1, 0)
patch4$c4m5l <- ifelse(patch4$c4m5>134.73, 1, 0)
patch4$c4m6l <- ifelse(patch4$c4m6>136.96, 1, 0)
patch4$c4m7l <- ifelse(patch4$c4m7>142.5, 1, 0)
patch4$c4m8l <- ifelse(patch4$c4m8>128.38, 1, 0)
patch4$c4m9l <- ifelse(patch4$c4m9>111.3, 1, 0)
patch4$c4m10l <- ifelse(patch4$c4m10>139.56, 1, 0)

table(patch4$c4m1l)+table(patch4$c4m2l)+table(patch4$c4m3l)+table(patch4$c4m4l)+table(patch4$c4m5l)+table(patch4$c4m6l)+table(patch4$c4m7l)+table(patch4$c4m8l)+table(patch4$c4m9l)+table(patch4$c4m10l)

#cluster5
results <- formula_rd(nomogram = nom.cox41)
patch5$c5m1 <- points_cal(formula = results$formula,rd=patch5)
results <- formula_rd(nomogram = nom.cox42)
patch5$c5m2 <- points_cal(formula = results$formula,rd=patch5)
results <- formula_rd(nomogram = nom.cox43)
patch5$c5m3 <- points_cal(formula = results$formula,rd=patch5)
results <- formula_rd(nomogram = nom.cox44)
patch5$c5m4 <- points_cal(formula = results$formula,rd=patch5)
results <- formula_rd(nomogram = nom.cox45)
patch5$c5m5 <- points_cal(formula = results$formula,rd=patch5)
results <- formula_rd(nomogram = nom.cox46)
patch5$c5m6 <- points_cal(formula = results$formula,rd=patch5)
results <- formula_rd(nomogram = nom.cox47)
patch5$c5m7 <- points_cal(formula = results$formula,rd=patch5)
results <- formula_rd(nomogram = nom.cox48)
patch5$c5m8 <- points_cal(formula = results$formula,rd=patch5)
results <- formula_rd(nomogram = nom.cox49)
patch5$c5m9 <- points_cal(formula = results$formula,rd=patch5)
results <- formula_rd(nomogram = nom.cox50)
patch5$c5m10 <- points_cal(formula = results$formula,rd=patch5)

patch5$c5m1l <- ifelse(patch5$c5m1>183.55, 1, 0)
patch5$c5m2l <- ifelse(patch5$c5m2>232.44, 1, 0)
patch5$c5m3l <- ifelse(patch5$c5m3>228.67, 1, 0)
patch5$c5m4l <- ifelse(patch5$c5m4>173.7, 1, 0)
patch5$c5m5l <- ifelse(patch5$c5m5>191.75, 1, 0)
patch5$c5m6l <- ifelse(patch5$c5m6>221.7, 1, 0)
patch5$c5m7l <- ifelse(patch5$c5m7>214.3, 1, 0)
patch5$c5m8l <- ifelse(patch5$c5m8>205.68, 1, 0)
patch5$c5m9l <- ifelse(patch5$c5m9>230.96, 1, 0)
patch5$c5m10l <- ifelse(patch5$c5m10>226.11, 1, 0)

table(patch5$c5m1l)+table(patch5$c5m2l)+table(patch5$c5m3l)+table(patch5$c5m4l)+table(patch5$c5m5l)+table(patch5$c5m6l)+table(patch5$c5m7l)+table(patch5$c5m8l)+table(patch5$c5m9l)+table(patch5$c5m10l)

#dist
map <- read.csv("F:/PRIVATE DATA/病理生境数据/Z0948231/map.csv", header = T, sep = ",")
mapr <- read.csv("G:/leading tile/Z0948231/patch_info.csv", header = T, sep = ",")
mapr <- subset(mapr, mapr$edge==1)
map1 <- subset(map, cluster==1)
map2 <- subset(map, cluster==2)
map3 <- subset(map, cluster==3)
map4 <- subset(map, cluster==4)
map5 <- subset(map, cluster==5)
vector2 <- mapr[, 6:7]
edv <- c()
for (i in 1:2809){
    ed <- c()
    vector1 <- map1[i, 3:4]
    for (j in 1:556){
        ed <- append(ed, sqrt(sum((vector1 - vector2[j,])^2)))
        edm <- min(ed)
    }
    edv <- append(edv, edm)
}
map1$dist <- edv
mean(map1$dist)
write.csv(map1, "F:/PRIVATE DATA/病理生境数据/Z0948231/map1.csv")
write.csv(map2, "F:/PRIVATE DATA/病理生境数据/Z0948231/map2.csv")
write.csv(map3, "F:/PRIVATE DATA/病理生境数据/Z0948231/map3.csv")
write.csv(map4, "F:/PRIVATE DATA/病理生境数据/Z0948231/map4.csv")
write.csv(map5, "F:/PRIVATE DATA/病理生境数据/Z0948231/map5.csv")

patch <- read.csv("F:/PRIVATE DATA/外部验证特征 省肿瘤/DC1803635/tumor feature.csv", sep = ",", header = F)
colnames(patch) <- cname
map <- read.csv("F:/PRIVATE DATA/外部验证特征 省肿瘤/DC1803635/map.csv", sep = ",", header = T)
patch <- merge(patch, map, by = "patch")
patch1 <- subset(patch, cluster == 1)
patch2 <- subset(patch, cluster == 2)
patch3 <- subset(patch, cluster == 3)
for (i in 2:2049){
    patch1[, i] <- (patch1[, i])/(mean(patch1[, i]))*(clinexters[1, (i-1)])
}
for (i in 2:2049){
    patch2[, i] <- (patch2[, i])/(mean(patch2[, i]))*(clinexters[2, (i-1)])
}
patchre <- rbind(patch1, patch2)
patchre <- rbind(patchre, patch3)

#特征标准化
patch <- read.csv("F:/PRIVATE DATA/外部验证特征 锦医/463353/tumor feature.csv", sep = ",", header = F)
colnames(patch) <- cname2
map <- read.csv("F:/PRIVATE DATA/外部验证特征 锦医/463353/map.csv", sep = ",", header = T)
patch <- merge(patch, map, by = "patch")
patch1 <- subset(patch, cluster == 1)
patch2 <- subset(patch, cluster == 2)
patch3 <- subset(patch, cluster == 3)
patch4 <- subset(patch, cluster == 4)
patch5 <- subset(patch, cluster == 5)
patch6 <- subset(patch, cluster == 6)
patch7 <- subset(patch, cluster == 7)
for (i in 2:2049){
  patch1[, i] <- (patch1[, i])/(mean(patch1[, i]))*(clinexter1s[171, (i-1)])
}
for (i in 2:2049){
  patch2[, i] <- (patch2[, i])/(mean(patch2[, i]))*(clinexter1s[172, (i-1)])
}
for (i in 2:2049){
  patch3[, i] <- (patch3[, i])/(mean(patch3[, i]))*(clinexter1s[173, (i-1)])
}
for (i in 2:2049){
  patch4[, i] <- (patch4[, i])/(mean(patch4[, i]))*(clinexter1s[174, (i-1)])
}
for (i in 2:2049){
  patch5[, i] <- (patch5[, i])/(mean(patch5[, i]))*(clinexter1s[175, (i-1)])
}
for (i in 2:2049){
  patch6[, i] <- (patch6[, i])/(mean(patch6[, i]))*(clinexter1s[176, (i-1)])
}
for (i in 2:2049){
  patch7[, i] <- (patch7[, i])/(mean(patch7[, i]))*(clinexter1s[177, (i-1)])
}
patchre <- rbind(patch1, patch2)
patchre <- rbind(patchre, patch3)
patchre <- rbind(patchre, patch4)
patchre <- rbind(patchre, patch5)
patchre <- rbind(patchre, patch6)
patchre <- rbind(patchre, patch7)
nacol <- c()
for (i in 2:2049){
  if (TRUE %in% is.na(patchre[, i])){
    nacol <- append(nacol, i)
  }
}
patchre <- patchre[, -nacol]

#cellsum
cellnum <- c()
cellsub <- subset(cell, ImageNumber == 1)
cellnum <- append(cellnum, nrow(cellsub))
cellsubmean <- c()
for (j in 12:46){
    cellsubmean <- append(cellsubmean, mean(cellsub[, j]))
}
for (i in 2:100){
    cellsub <- subset(cell, ImageNumber == i)
    cellnum <- append(cellnum, nrow(cellsub))
    cellsubmeans <- c()
    for (j in 12:46){
        cellsubmeans <- append(cellsubmeans, mean(cellsub[, j]))
    }
    cellsubmean <- rbind(cellsubmean, cellsubmeans)
}
cellsubmean <- cbind(cellnum, cellsubmean)
cellsubmean <- as.data.frame(cellsubmean)
colnames(cellsubmean) <- cellname

cellname <- colnames(cell)[12:46]
cellname <- append("cellnum", cellname)
colnames(cellsubmean) <- cellname