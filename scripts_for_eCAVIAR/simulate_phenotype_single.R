my_simulate<-function(file,num,causal_list,output,ld_file)
 {
   #num=4
   #file<-"combined_ZZZ3_genotype_phenotype"
   #ld_file="ld_file"
   #output<-"output"
   #causal_list<-"causal"
   data=read.table(file,header=T,na.strings="NA")
   data<-data[,c(1:102)]
#   data=data[,-2]
   col=colnames(data)
   num=1
   num<-as.integer(num)
   a=runif(num,min=3,max=length(col))
   beta=runif(num,min=0.05,max=0.1)
   a=as.integer(a+0.5)
#   cat("a is: ",a,"\n")
 #  cat("beta is: ",beta,"\n")
   j=1
   LD_r2<-matrix(NA,nrow=length(a),ncol=length(a))
   
   LD_r2_all<-matrix(NA,nrow=dim(data)[2]-2,ncol=dim(data)[2]-2)
   for(i in 1:(dim(data)[2]-2))
    {
      for(j in 1:(dim(data)[2]-2))
       {
         r<-(cor(data[,i+2],data[,j+2],use="complete"))
         LD_r2_all[i,j]<-r
         LD_r2_all[j,i]<-r
       }
    }
   write.table(LD_r2_all,file=ld_file,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
 #  mean_LD_r2=apply(LD_r2, 2, function(x) sum(x)-1)
 #  mean_LD_r2=mean_LD_r2/length(a)
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
 #  for(i in 2:(length(a)))
 #   {
    #  allel_eff[i]=0
    #  phe=phe+data[,a[i]]*allel_eff[i]
   phe=phe+data[,a[1]]*allel_eff[1]
   phe1=phe1+data[,a[1]]*allel_eff[1]
  #  }
  # cat("Right now, it is Ok\n")
   causal_SNP=data.frame(col[a],frq,allel_eff)
   colnames(causal_SNP)=c("Causal_SNP","frequency","Allele_effect")
   write.table(causal_SNP,file=causal_list,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
 #  background=rnorm(length(rownames(data)),mean=0,sd=1)
 #  phe=phe+background
   data1=data.frame(phe,phe1,data)
   colnames(data1)<-c("phe1","phe2",colnames(data))
   pheno_file=paste(output,"_simulated",sep="")
   write.table(data1,file=pheno_file,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))

   output_LD_r2<-paste(output,"_LD_r2",sep="")
  # cat("Right now, it is Ok\n")
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
           #  sd<-sd(data1$phe1,na.rm=TRUE)/sqrt(dim(data)[2]-3)
             z_score1=summary(fit)$coefficients[2,3]
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
            z_score2=summary(fit)$coefficients[2,3]
            
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
   command=paste("/nv/hp10/bzeng30/data2/software/caviar-master/CAVIAR-C++/eCAVIAR", " -o ",output," -r 0.9 -l ", ld_file," -l "  ,ld_file, "-c 3  -z ",output1," -z ",output2,sep=" ")
  # cat("Command is: ",command,"\n")
   system(command)
   smr_out<-paste(output,"_col",sep="")
   smr_result<-read.table(smr_out)
  # cat("Come here1\nCausal variants are: ",col[a],"\n")
   smr_causal_result<-smr_result[which(smr_result[,1] %in% col[a]),]
  # cat("Come here2\n")
   output_causal<-paste(output,"_CLPP_causl",sep="")
   write.table(smr_causal_result,output_causal,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
   combined_output<-paste(output_causal,"_combined_with_causal_infor",sep="") 
   command<-paste("perl  /nv/hp10/bzeng30/data2/simulation_with_random_picked_genes_to_evaluate_influence_on_CLPP/scripts/calculate_pp_for_CAVIAR2/another/combined_effect_size_and_pp_CLPP.pl  --causal",causal_list," --clpp ",smr_out," --output ",combined_output,sep=" ")
   system(command) 
   command<-paste("rm ",output,"*_set ",output,"*_post ",output,"*_simulated ",output,"*_CLPP_causl ",output,"*_col ",output,"*_? ","*_causal",sep="")
   system(command) 
 }
Args <- commandArgs()

cat("Useage: R --vanilla --slave --input infile --output outfile  < simulate_phenotype.R\n")

for (i in 1:length(Args))
 {
    if (Args[i] == "--input") file_name = Args[i+1]
    if (Args[i] == "--output") out_file = Args[i+1]
    if (Args[i] == "--num") num = Args[i+1]
    if (Args[i] == "--list")   causal_list     = Args[i+1]
    if (Args[i] == "--ld_file") ld_file = Args[i+1]
 }
my_simulate(file_name,num,causal_list,out_file,ld_file)

