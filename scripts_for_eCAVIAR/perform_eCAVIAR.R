my_simulate<-function(file,number,causal_list,output,ld_file)
 {

  # output_LD_r2<-paste(output,"_LD_r2",sep="")
  # output<-"output"
  # number=1
   cat("Right now, it is Ok\n")
  # file="./output_simulated_phenotype_1"
   data1=read.table(file,header=T)
   col1<-colnames(data1)
   z_out1<-c("snp","z_score")
   z_out2<-c("snp","z_score")
   for(i in 4:dim(data1)[2])
     {
           # cat("Coming\n")
            string=paste(col1[1],"~",col1[i],sep="")
            cat("string is: ",string,"\n")
            fit<-lm(data=data1,formula = as.formula(paste(col1[1],"~",col1[i],sep="")),na.action=na.exclude)
          #  fit1=fit
          #  num<-length(summary(fit)$coefficients[,1])
           # for(q in 2:num)
           #  {
               #p=summary(fit)$coefficients[2,4]
             alle=summary(fit)$coefficients[2,1]
             cat("alle is: ",alle,"\n")
             sd<-sd(data1$phe1,na.rm=TRUE)/sqrt(dim(data1)[2]-3)
             z_score1=alle/sd
             cat("col1i is: ",col1[i], " and z_score1 is: ",z_score1,"\n")
             z_out1<-rbind(z_out1,paste(c(col1[i],z_score1),sep=" "))
           #    snp_test=col1[i]
           #  }
               # cat("Coming\n")
            string=paste(col1[2],"~",col1[i],sep="")
         #   cat("string is: ",string,"\n")
            fit<-lm(data=data1,formula = as.formula(paste(col1[2],"~",col1[i],sep="")),na.action=na.exclude)
            
            sd<-sd(data1$phe2,na.rm=T)/sqrt(dim(data1)[2]-3)
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
           cat("z_out1 is: ",z_out1,"\n")
           cat("z_out2 is: ",z_out2,"\n")
          
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
   command=paste("/nv/hp10/bzeng30/data2/software/caviar-master/CAVIAR-C++/eCAVIAR", " -o ",output," -r 0.9 -l ", ld_file," -l "  ,ld_file, "-c  ", number," -z ",output1," -z ",output2,sep=" ")
   cat("Command is: ",command,"\n")
   system(command)
   smr_out<-paste(output,"_col",sep="")
   #smr_result<-read.table(smr_out)
   #cat("Come here1\nCausal variants are: ",col[a],"\n")
   #smr_causal_result<-smr_result[which(smr_result[,1] %in% col[a]),]
   #cat("Come here2\n")
   #output_causal<-paste(output,"_CLPP_causl",sep="")
   #write.table(smr_causal_result,output_causal,quote = FALSE, sep = "\t", eol = "\n", na = "0", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
   combined_output<-paste(output,"_CLPP_combined_with_causal_infor",sep="") 
   command<-paste("perl ~/data/my_script/combined_effect_size_and_CLPP.pl --causal",causal_list," --clpp ",smr_out," --output ",combined_output,sep=" ")
   cat("command is: ",command,"\n")
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

