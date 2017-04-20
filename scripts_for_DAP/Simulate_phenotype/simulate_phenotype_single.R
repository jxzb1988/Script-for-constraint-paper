my_simulate<-function(file,num,causal_list,output,ld_file,weight_file,ld_filter)
 {
      data=read.table(file,header=T,na.strings="NA")
#      data=data[,-2]
      col=colnames(data)
      num<-as.integer(num)
      rand=sample(3:(length(col)-100))
   #   a=sample(rand:(rand+100),num,replace=F)
      a=sample(3:length(col),num,replace=F)
      beta=runif(num,min=0.02,max=0.1)
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
      LD_r2_all<-matrix(NA,nrow=dim(data)[2]-1,ncol=dim(data)[2]-1)
      for(i in 1:(dim(data)[2]-1))
       {
         for(j in 1:(dim(data)[2]-1))
          {
            LD_r2_all[i,j]=(cor(data[,i+1],data[,j+1],use="complete"))
            LD_r2_all[j,i]=(cor(data[,i+1],data[,j+1],use="complete"))
          }
       }
      if(!is.na(ld_filter))
       {
         time=0
         while(min(abs(LD_r2))<ld_filter && time<=1000)
          {
            time=time+1
            a=sample(3:length(col),num,replace=F)
        #    a=sample(rand:(rand+100),num,replace=F)
            for(i in 1:(length(a)))
             {
               for(j in i:length(a))
                {
                  LD_r2[i,j]=(cor(data[,a[i]],data[,a[j]],use="complete"))^2
                  LD_r2[j,i]=(cor(data[,a[i]],data[,a[j]],use="complete"))^2
                }
             }
          }
        if(time>1000)
         {
           return(0)
         }
       }
   write.table(LD_r2_all,file=ld_file,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
   frq=array()
   allel_eff=array()
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
      allel_eff[i]<-sqrt(beta[i]/(2*f*(1-f)))
    }
   cat("frq is: ",frq,"\n")
   background=rnorm(length(rownames(data)),mean=0,sd=1)
   background1=rnorm(length(rownames(data)),mean=0,sd=1)
   phe=background
   phe1=background1
   allel_eff<-array()
   
   for(i in 1:length(a))
    {
     #snp_test<-a[i]
     allel_eff[i]<-sqrt(beta[i]/(2*frq[i]*(1-frq[i])))
    }
  # phe=phe+data[,a[1]]*allel_eff[1]
  # phe1=phe1+data[,a[1]]*allel_eff[1]
   for(i in 1:(length(a)))
    {
    #  allel_eff[i]=0
    #  phe=phe+data[,a[i]]*allel_eff[i]
      phe=phe+data[,a[i]]*allel_eff[i]
      phe1=phe1+data[,a[i]]*allel_eff[i]
    }
   cat("Right now, it is Ok\n")
   causal_SNP=data.frame(col[a],frq,allel_eff)
   colnames(causal_SNP)=c("Causal_SNP","frequency","Allele_effect")
   write.table(causal_SNP,file=causal_list,quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
 #  return()
 #  background=rnorm(length(rownames(data)),mean=0,sd=1)
 #  phe=phe+background
 #  data1=data.frame(phe,data)
   data[,2]=phe
 #  colnames(data1)<-c("phe1","phe2",colnames(data))
   pheno_file=paste(output,"_simulated",sep="")
   write.table(data,file=pheno_file,quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
   return()
   output_LD_r2<-paste(output,"_LD_r2",sep="")
   cat("Right now, it is Ok\n")
   col1<-colnames(data1)
   z_out1<-c("snp","z_score")
   z_out2<-c("snp","z_score")
   for(i in 5:dim(data1)[2])
     {
           # cat("Coming\n")
            string=paste(col1[1],"~",col1[i],sep="")
         #   cat("string is: ",string,"\n")
            fit<-lm(data=data1,formula = as.formula(paste(col1[1],"~",col1[i],sep="")),na.action=na.exclude)
          #  fit1=fit
          #  num<-length(summary(fit)$coefficients[,1])
           # for(q in 2:num)
           #  {
               #p=summary(fit)$coefficients[2,4]
             alle=summary(fit)$coefficients[2,1]
             sd<-sd(data1$phe1,na.rm=TRUE)/(dim(data)[2]-3)
             z_score1=alle/sd
             z_out1<-rbind(z_out1,paste(c(col1[i],z_score1),sep=" "))
           #    snp_test=col1[i]
           #  }
               # cat("Coming\n")
            string=paste(col1[2],"~",col1[i],sep="")
         #   cat("string is: ",string,"\n")
            fit<-lm(data=data1,formula = as.formula(paste(col1[2],"~",col1[i],sep="")),na.action=na.exclude)
            
            sd<-sd(data1$phe2,na.rm=T)/sqrt(dim(data)[2]-3)
          #  fit1=fit
          #  num<-length(summary(fit)$coefficients[,1])
          #  for(q in 2:num)
          #   {
            alle=summary(fit)$coefficients[2,1]
            z_score2=alle/sd
            
            z_out2<-rbind(z_out2,paste(c(col1[i],z_score2),sep=" "))
          #     alle=summary(fit)$coefficients[q,1]
          #     snp_test=col1[i]
          #   }

          
     }
   #command=paste("/nv/hp10/bzeng30/data2/software/caviar-master/CAVIAR-C++/eCAVIAR", " -o ",output," -l ", "-l"  , " -z "," -z")
   output1=paste(output,"_1",sep="")
   output2=paste(output,"_2",sep="")
   z_out1<-z_out1[-1,]
   z_out2<-z_out2[-1,]
   write.table(z_out1,file=output1,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
   write.table(z_out2,file=output2,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
  # colnames(LD_r2)<-colnames(data)[a]
  # LD_r2_result<-cbind(colnames(data)[a],beta,allel_eff,mean_LD_r2)
  # colnames(LD_r2_result)<-c("Causal","Effect_size","Allelic","Mean_LD")
   #write.table(LD_r2_result,file=output_LD_r2,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
   command=paste("/nv/hp10/bzeng30/data2/software/caviar-master/CAVIAR-C++/eCAVIAR", " -o ",output," -l ", ld_file," -l "  ,ld_file, "-c 2  -z ",output1," -z ",output2,sep=" ")
   cat("Command is: ",command,"\n")
   system(command)
   smr_out<-paste(output,"_col",sep="")
   smr_result<-read.table(smr_out)
   cat("Come here1\nCausal variants are: ",col[a],"\n")
   smr_causal_result<-smr_result[which(smr_result[,1] %in% col[a]),]
   cat("Come here2\n")
   output_causal<-paste(output,"_CLPP_causl",sep="")
   write.table(smr_causal_result,output_causal,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
   combined_output<-paste(output_causal,"_combined_with_causal_infor",sep="") 
   command<-paste("perl ~/data/my_script/combined_effect_size_and_CLPP.pl --causal",causal_list," --clpp ",output_causal," --output ",combined_output,sep=" ")
   system(command) 
    
 }
Args <- commandArgs()

cat("Useage: R --vanilla --slave --input infile --output outfile  < simulate_phenotype.R\n")
weight=NA
ld_filter=NA
for (i in 1:length(Args))
 {
    if (Args[i] == "--input") file_name = Args[i+1]
    if (Args[i] == "--output") out_file = Args[i+1]
    if (Args[i] == "--num") num = Args[i+1]
    if (Args[i] == "--list")   causal_list     = Args[i+1]
    if (Args[i] == "--ld_file") ld_file = Args[i+1]
    if (Args[i] == "--ld_filter") ld_filter = Args[i+1]
    if (Args[i] == "--weight") weight = Args[i+1]
 }
my_simulate(file_name,num,causal_list,out_file,ld_file,weight,ld_filter)

