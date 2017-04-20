my_simulate<-function(data,num,causal_list,output)
 {
   data=read.table(data,header=T,na.strings="NA")
   col=colnames(data)
   num<-as.integer(num)
   a=runif(num,min=2,max=length(col))
   beta=runif(num,min=0.02,max=0.1)
   a=as.integer(a+0.5)
 #  cat("a is: ",a,"\n")
 #  cat("beta is: ",beta,"\n")
   j=1
   LD_r2<-matrix(NA,nrow=length(a),ncol=length(a))
   for(i in 1:(length(a)))
    {
      for(j in i:length(a))
       {
         LD_r2[i,j]=(cor(data[,a[i]],data[,a[j]],use="complete"))^2
         LD_r2[j,i]=(cor(data[,a[i]],data[,a[j]],use="complete"))^2
       }
    }
   mean_LD_r2=apply(LD_r2, 2, function(x) sum(x)-1)
   mean_LD_r2=mean_LD_r2/length(a) 
   frq=array()
   for(i in 1:length(a))
    {
      arr=data[,a[i]]
   #   cat("arr is: ",arr,"\n")
      f=(2*length(which(arr==0))+length(which(arr==1)))/(2*length(arr))
      if(f>0.5)
       {
         f=1-f
       }
      frq[i]=f  
    }
 #  cat("frq is: ",frq,"\n")
   background=rnorm(length(rownames(data)),mean=0,sd=1)
   phe=background
   allel_eff<-array()
   allel_eff[1]=-sqrt(beta[1]/(2*frq[1]*(1-frq[1])))
   phe=phe+data[,a[1]]*allel_eff[1]
   allel_eff[2]=-sqrt(beta[2]/(2*frq[2]*(1-frq[2]))) 
   phe=phe+data[,a[2]]*allel_eff[2]
   for(i in 3:(length(a)))
    {
      allel_eff[i]=sqrt(beta[i]/(2*frq[i]*(1-frq[i])))
      phe=phe+data[,a[i]]*allel_eff[i]
    #  phe=data[,a[1]]*allel_eff[1]+data[,a[2]]*allel_eff[2]+data[,a[3]]*allel_eff[3]+data[,a[4]]*allel_eff[4]
    }
   causal_SNP=data.frame(col[a],frq,allel_eff)
   colnames(causal_SNP)=c("Causal_SNP","frequency","Allele_effect")
   write.table(causal_SNP,file=causal_list,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
 #  background=rnorm(length(rownames(data)),mean=0,sd=1)
   phe=phe+background
   data1=data.frame(phe,data[,1],data[,a])
   colnames(data1)<-c("phe",colnames(data)[1],colnames(data[,a]))
   output_LD_r2<-paste(output,"_LD_r2",sep="")
   colnames(LD_r2)<-colnames(data)[a]
   LD_r2_result<-cbind(colnames(data)[a],beta,allel_eff,mean_LD_r2)
   colnames(LD_r2_result)<-c("Causal","Effect_size","Allelic","Mean_LD")
   write.table(data1,file=output,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
   write.table(LD_r2_result,file=output_LD_r2,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
 }
Args <- commandArgs()

cat("Useage: R --vanilla --slave --input infile --output outfile  < simulate_phenotype.R\n")

for (i in 1:length(Args))
 {
    if (Args[i] == "--input") file_name = Args[i+1]
    if (Args[i] == "--output") out_file = Args[i+1]
    if (Args[i] == "--num") num = Args[i+1]
    if (Args[i] == "--list")   causal_list     = Args[i+1]
 }
my_simulate(file_name,num,causal_list,out_file)

