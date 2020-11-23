#Project code
#Cory and Blake
#Continuous Clustering Simulations

library(MASS);
library(clValid);
library(mixture);
library(MixGHD);
library(teigen);
library(mixsmsn);
library(meanShiftR);

###########################################################
#Simulate Data#############################################
###########################################################

#parameters
mu1 <- c(0,0)
mu<- list(c(2.5,2.5),c(1,1),c(.5,0))
ss<- list(diag(c(1,1)),matrix(c(1,.75,.75,1),2),matrix(c(1,-.75,-.75,1),2),
          matrix(c(4,4.5,4.5,9),2),matrix(c(4,-4.5,-4.5,9),2) )
skew<- list(c(0,0),c(1,-1),c(1,1))  

#fewer models just for now. Remove later
#mu<- list(c(2.5,2.5),c(1,1))
#ss<- list( diag(c(1,1)), diag(c(2,2))  )
#skew<- list(c(0,0),c(1,1))  

#results matrix
results<-matrix(0,nrow = 675, ncol = 17) #setting, ARI1, Dunn1
colnames(results)<-c('setting','ARI 1','Dunn 1',
                     'ARI2','Dunn2','BIC2','ARI3','Dunn3','BIC3',
                     'ARI4', 'nARI4', 'Dunn4', 
                     'nDunn4','ARI5','Dunn5','ARI6','Dunn6')

i=0
for(mu2 in mu){
    for(s1 in ss){
        for(s2 in ss){
            for(skew1 in skew){
                for (skew2 in skew){
                    
                    ############
                    #Makes data#
                    ############
                    
                    arg1=list(mu=mu1,Sigma=s1,shape=skew1)
                    arg2=list(mu=mu2,Sigma=s2,shape=skew2)
                    args=list(arg1,arg2)
                    stuff=rmmix(n=200,p=c(.5,.5),family = 'Skew.normal',arg=args,
                                cluster = TRUE)
                    data<-stuff$y
                    labels <- stuff$cluster
                    #plot(data,col = labels)
                    i=i+1
                    if (i%%25==0){print(i)}
                    results[i,1] <- i 
                    
                    ########
                    #Models#
                    ########
                    
                    #K-Means (1)
                    #ARI, Dunn
                    model1 = kmeans(data, nstart = 25, centers = 2)
                    results[i,2] <- ARI(model1$cluster,labels)
                    results[i,3]<-dunn(clusters = model1$cluster, Data = data)
                    
                    #GMM (2)
                    #ARI, Dunn, BIC
                    model2 = gpcm(data, G=2,mnames = NULL)
                    results[i,4] <- ARI(model2$map,labels)
                    results[i,5] <- dunn(clusters = model2$map, Data = data)
                    results[i,6] <- model2$bicModel$bic
                    
                    #t-mix (3)
                    #ARI, Dunn, BIC
                    tmix <- teigen(x=data,Gs=2,models="UUUU");
                    class.tmix <- tmix$classification 
                    results[i,7] <- ARI(class.tmix,labels)
                    results[i,8] <- dunn(clusters = class.tmix, Data = data)
                    results[i,9] <- tmix$bic
                    
                    #meanshift
                    #ARI, nARI, Dunn, nDunn
                    ls<-c()
                    lsd<-c()
                    for (j in 2:10){
                        model <- meanShift(data, trainData = data, nNeighbors = j,iterations=25)
                        ls<-c(ls,ARI(labels, model$assignment))
                        ddd <- dunn(clusters = model$assignment,Data=data)
                        if (ddd != Inf){
                            lsd<-c(lsd,dunn(clusters = model$assignment,Data=data))
                        } else {lsd<-c(lsd,0)}
                        #print( dunn(clusters = model$assignment,Data=data) )
                        #lsd<-c(lsd,dunn(clusters = model$assignment,Data=data))
                    }
                    results[i,10]<-max(ls)
                    results[i,11]<-which.max(ls)
                    results[i,12]<-max(lsd)
                    results[i,13]<-which.max(lsd)
                    #print(max(ls))
                    #print(which.max(ls))
                    
                    #Ward
                    #ARI
                    #Dunn
                    modelw<-cutree(agnes(data,method = 'ward'),2)
                    results[i,14] <- ARI(modelw,labels)
                    results[i,15] <- dunn(clusters = modelw, Data = data)              
                    
                    #PAM
                    #ARI Dunn
                    modelp <- pam(data,k=2)
                    results[i,16] <- ARI(modelp$clustering,labels)
                    results[i,17] <- dunn(clusters = modelp$clustering, Data = data)               
                    
                    
                }
            }
        }
    }
}

write.csv(results,'results_df.csv', row.names = FALSE)

###########################################################
#Parameter Table###########################################
###########################################################

#parameters
mu1 <- c(0,0)
mu<- list(c(2.5,2.5),c(1,1),c(.5,0))
ss<- list(diag(c(1,1)),matrix(c(1,.75,.75,1),2),matrix(c(1,-.75,-.75,1),2),matrix(c(4,4.5,4.5,9),2),matrix(c(4,-4.5,-4.5,9),2) )
skew<- list(c(0,0),c(1,-1),c(1,1))  


#results matrix
param_mat<-matrix(0,nrow = 675, ncol = 17) 
colnames(param_mat)<-c('G1_mu_x1','G1_mu_x2','G2_mu_x1','G2_mu2_x2','dist_between_means','G1_var_x1','G1_var_x2','G1_covar','G1_corr','G2_var_x1','G2_var_x2','G2_covar','G2_corr','G1_skew_x1','G1_skew_x2','G2_skew_x1','G2_skew_x2')

i=0
for(mu2 in mu){
    for(s1 in ss){
        for(s2 in ss){
            for(skew1 in skew){
                for (skew2 in skew){
                    i=i+1
                    param_mat[i,3]<-mu2[1]
                    param_mat[i,4]<-mu2[2]
                    param_mat[i,5]<-sqrt( sum((mu1-mu2)^2)        )
                    param_mat[i,6]<-s1[1,1]
                    param_mat[i,7]<-s1[2,2]
                    param_mat[i,8]<-s1[1,2]
                    param_mat[i,9]<-s1[1,2]/(sqrt(s1[1,1]*s1[2,2]))
                    param_mat[i,10]<-s2[1,1]
                    param_mat[i,11]<-s2[2,2]
                    param_mat[i,12]<-s2[1,2]
                    param_mat[i,13]<-s2[1,2]/(sqrt(s2[1,1]*s2[2,2]))
                    param_mat[i,14]<-skew1[1]
                    param_mat[i,15]<-skew1[2]
                    param_mat[i,16]<-skew2[1]
                    param_mat[i,17]<-skew2[2]
                    
                    
                }
            }
        }
    }
}

write.csv(param_mat,'parameter_matrix.csv', row.names = FALSE)

###########################################################
#Interpreting##############################################
###########################################################

params<-read.csv('parameter_matrix.csv',header=TRUE)
results<-read.csv('results_df.csv')

#how does the dist between means affect the scores
cor(params[,'dist_between_means'],results[,-1])
cor(params[,'G1_corr'],results[,-1])
cor(params[,'G2_corr'],results[,-1])

##############
#ARI Boxplots#
##############

boxplot(results[,c('ARI.1','ARI2','ARI3','ARI4','ARI5','ARI6')],
        main = 'ARIs for all simulations',
        names =c('K-means','GMM','t-mix','Meanshift','Ward','PAM'))

boxplot(results[params$G2_mu_x1==2.5,c('ARI.1','ARI2','ARI3','ARI4','ARI5','ARI6')],
        main = 'ARIs for far groups')

boxplot(results[params$G2_mu_x1==1,c('ARI.1','ARI2','ARI3','ARI4','ARI5','ARI6')],
        main = 'ARIs for medium distance groups',
        names =c('K-means','GMM','t-mix','Meanshift','Ward','PAM'))

boxplot(results[params$G2_mu_x1==0.5,c('ARI.1','ARI2','ARI3','ARI4','ARI5','ARI6')],
        main = 'ARIs for Close Groups',
        names =c('K-means','GMM','t-mix','Meanshift','Ward','PAM'))

######
#Skew#
######

rows = params$G1_skew_x1==0 & params$G1_skew_x2==0 & params$G2_skew_x1==0 & params$G2_skew_x2==0
boxplot(results[rows,c('ARI.1','ARI2','ARI3','ARI4','ARI5','ARI6')],main = 'ARIs for no skew')

#Some skew
rows = params$G1_skew_x1==1 & params$G1_skew_x2==1 & params$G2_skew_x1==1 & params$G2_skew_x2==1
boxplot(results[rows,c('ARI.1','ARI2','ARI3','ARI4','ARI5','ARI6')])

rows = params$G1_skew_x1==1 & params$G1_skew_x2==-1 & params$G2_skew_x1==1 & params$G2_skew_x2==-1
boxplot(results[rows,c('ARI.1','ARI2','ARI3','ARI4','ARI5','ARI6')],main = 'ARIs for close groups')

rows = params$G1_skew_x1==0 & params$G1_skew_x2==0 & params$G2_skew_x1==1 & params$G2_skew_x2==1
boxplot(results[rows,c('ARI.1','ARI2','ARI3','ARI4','ARI5','ARI6')],
        main = 'ARI for skew1=(0,0), skew2=(1,1)',
        names =c('K-means','GMM','t-mix','Meanshift','Ward','PAM'))



#crossing groups
rows = params$G1_corr*params$G2_corr < 0
boxplot(results[rows,c('ARI.1','ARI2','ARI3','ARI4','ARI5','ARI6')],
        main = 'Crossing Distributions',
        names =c('K-means','GMM','t-mix','Meanshift','Ward','PAM'))

arg1=list(mu=c(0,0),Sigma=matrix(c(1,.75,.75,1),2),shape=c(0,0))
arg2=list(mu=c(1,1),Sigma=matrix(c(1,-.75,-.75,1),2),shape=c(0,0))
args=list(arg1,arg2)
stuff=rmmix(n=200,p=c(.5,.5),family = 'Skew.normal',arg=args,
            cluster = TRUE)
plot(stuff$y,col = stuff$cluster,xlab = 'x1',ylab = 'x2')

############
#Dunn Index#
############

boxplot(results[,c('Dunn.1','Dunn2','Dunn3','Dunn4','Dunn5','Dunn6')],
        main = 'Dunn Index for all simulations',
        names =c('K-means','GMM','t-mix','Meanshift','Ward','PAM'))


#####
#BIC#
#####
a = results$BIC2
b = abs(results$BIC3)
boxplot(a,b,
        main = 'BIC for all simulations',
        names =c('GMM','t-mix'))

a = results$BIC2[params$G2_mu_x1==2.5]
b = results$BIC2[params$G2_mu_x1==1]
c = results$BIC2[params$G2_mu_x1==0.5]
boxplot(a,b,c,
        main = 'BIC for GMM based on distance between means',
        names =c('Far','Medium','Close'))

a = abs(results$BIC3[params$G2_mu_x1==2.5])
b = abs(results$BIC3[params$G2_mu_x1==1])
c = abs(results$BIC3[params$G2_mu_x1==0.5])
boxplot(a,b,c,
        main = 'BIC for t-mix based on distance between means',
        names =c('Far','Medium','Close'))

###########################################
#Correlation between ARI and other metrics#
###########################################

cor(results[,c('ARI.1','Dunn.1')])
cor(results[,c('ARI2','Dunn2','BIC2')]) #GMM BIC and ARI are basically independent 
cor(results[,c('ARI3','Dunn3','BIC3')])#t-mix is the best
cor(results[,c('ARI4','Dunn4')])
cor(results[,c('ARI5','Dunn5')])
cor(results[,c('ARI6','Dunn6')])



#Meanshift
#number of groups
boxplot(results[,c('nARI4','nDunn4')],
        main = 'Number of Groups for Mean Shift',
        names=c('n based on ARI','n based on Dunn'))

print(sum(results$nARI4 == results$nDunn4)/675)
print(sum(results$nARI4 != results$nDunn4)/675)


########################################################################
#finding the number of groups for meanshift#############################
########################################################################

#parameters
mu1 <- c(0,0)
mu<- list(c(2.5,2.5),c(1,1),c(.5,0))
ss<- list(diag(c(1,1)),matrix(c(1,.75,.75,1),2),matrix(c(1,-.75,-.75,1),2),
          matrix(c(4,4.5,4.5,9),2),matrix(c(4,-4.5,-4.5,9),2) )
skew<- list(c(0,0),c(1,-1),c(1,1))  

#results matrix
results<-matrix(0,nrow = 675, ncol = 17) #setting, ARI1, Dunn1
colnames(results)<-c('setting','ARI 1','Dunn 1',
                     'ARI2','Dunn2','BIC2','ARI3','Dunn3','BIC3',
                     'ARI4', 'nARI4', 'Dunn4', 
                     'nDunn4','ARI5','Dunn5','ARI6','Dunn6')

N_ARI<-c()
N_Dunn<-c()

i=0
for(mu2 in mu){
    for(s1 in ss){
        for(s2 in ss){
            for(skew1 in skew){
                for (skew2 in skew){
                    
                    ############
                    #Makes data#
                    ############
                    
                    arg1=list(mu=mu1,Sigma=s1,shape=skew1)
                    arg2=list(mu=mu2,Sigma=s2,shape=skew2)
                    args=list(arg1,arg2)
                    stuff=rmmix(n=200,p=c(.5,.5),family = 'Skew.normal',arg=args,
                                cluster = TRUE)
                    data<-stuff$y
                    labels <- stuff$cluster
                    #plot(data,col = labels)
                    i=i+1
                    if (i%%25==0){print(i)}
                    results[i,1] <- i 
                    
                    ########
                    #Models#
                    ########
                    
                    #meanshift
                    #ARI, nARI, Dunn, nDunn
                    ls<-c()
                    lsd<-c()
                    n_groups_ARI <-c()
                    n_groups_Dunn <-c()
                    for (j in 2:10){
                        model <- meanShift(data, trainData = data, nNeighbors = j,iterations=25)
                        ls<-c(ls,ARI(labels, model$assignment))
                        n_groups_ARI<-c(n_groups_ARI,max(model$assignment))
                        ddd <- dunn(clusters = model$assignment,Data=data)
                        if (ddd != Inf){
                            lsd<-c(lsd,dunn(clusters = model$assignment,Data=data))
                        } else {
                            lsd<-c(lsd,0)
                        }
                    }
                    N_ARI<-c(N_ARI, n_groups_ARI[which.max(ls)])
                    N_Dunn<-c(N_Dunn, n_groups_ARI[which.max(lsd)])
                    
                    #print(max(ls))
                    #print(which.max(ls))
                    
                }
            }
        }
    }
}

boxplot(N_ARI,N_Dunn,
        main='Number of groups',
        names=c('ARI','Dunn'))

boxplot(N_ARI,N_Dunn,
        main='Number of groups',
        names=c('ARI','Dunn'),
        ylim = c(0,10))

sum(N_ARI == 2)/675
sum(N_Dunn == 2)/675



