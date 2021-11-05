library("hadron")
betas<-c(0.050000,0.100000,0.200000,0.300000,0.400000,0.600000,0.800000,1.000000,1.500000,2.000000,2.500000,3.000000,3.500000,4.000000,4.500000,5.000000,5.500000,1.200000,1.400000,1.600000,1.800000)

write.csv("#0", "resulthadron.csv") #makes sure to write to an empty file
for (beta in betas){
    betaprint<-rep(beta, times=16)
    string<-sprintf("measbeta%fNt10Ns10.txt", beta)
    file<-read.table(string)
    data<-file$V1
    message(beta,"\t",1-mean(data),"\t", sd(data)/sqrt(length(data)))
    result<-bootstrap.analysis(file$V1, boot.R=1000, tsboot.sim=fixed)
    result$beta<-c(betaprint)
#~     result<-cbind(result, c(betaprint))
    write.table(result, "resulthadron.csv", append=TRUE, col.names=FALSE)#sep=" "
#~     print(result)
}
print(warnings())
