my_simulate<-function(data,list,output)
 {
   data=read.table(data,header=T,na.strings="NA")
   col=colnames(data)
   a=runif(4,min=3,max=length(col))
   beta=runif(4,min=0.05,max=0.1)
   a=as.integer(a+0.5)
 #  cat("a is: ",a,"\n")
 #  cat("beta is: ",beta,"\n")
   j=1
   frq=array()
   for(i in a)
    {
      arr=data[,i]
   #   cat("arr is: ",arr,"\n")
      f=(2*length(which(arr==0))+length(which(arr==1)))/(2*length(arr))
      if(f>0.5)
       {
         f=1-f
       }
      frq[j]=f
      j=j+1  
    }
 #  cat("frq is: ",frq,"\n")
   allel_eff=sqrt(beta/(2*frq*(1-frq)))
   phe=data[,a[1]]*allel_eff[1]+data[,a[2]]*allel_eff[2]+data[,a[3]]*allel_eff[3]+data[,a[4]]*allel_eff[4]
   causal_SNP=data.frame(col[a],frq,allel_eff)
   colnames(causal_SNP)=c("Causal_SNP","frequency","Allele_effect")
   write.table(causal_SNP,file=list,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
   background=rnorm(length(rownames(data)),mean=0,sd=1)
   phe=phe+background
   data1=data.frame(phe,data)
   colnames(data1)<-c("phe",colnames(data))
   write.table(data1,file=output,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
 }
Args <- commandArgs()

cat("Useage: R --vanilla --slave --input infile --output outfile  < simulate_phenotype.R\n")

for (i in 1:length(Args))
 {
    if (Args[i] == "--input") file_name = Args[i+1]
    if (Args[i] == "--output") out_file = Args[i+1]
    if (Args[i] == "--list")   list     = Args[i+1]
 }
my_simulate(file_name,list,out_file)

