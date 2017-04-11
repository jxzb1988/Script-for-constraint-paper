
my_conditional<-function(file,output)
 {
   library(fmsb)
   data=read.table(file,header=T,row.names=1,na.strings="NA")
   dim<-dim(data)
   len<-dim[2]
   out_text=""
   col=colnames(data)
   j=0
   for(fea in col)
    {
      if(grepl ("IL",fea)==1)
       {
         j=j+1
       }
      if(grepl ("phe",fea)==1)
       {
         j=j+1
       }
    }
  # cutoff=0.05/(len-j)
   cutoff=0.0001
   for(t in 1:j)
    {
      cat("I am a\n")
      x=0
      bi=0
      left<-array()
      target=0
      p_target=1
      data1=data
   #   data2=data
      for(k in (j+1):len)
       {         
          cat("I am in\n")
      #   cat("*******************************************************Run",k,"***************************","\n")
         col1<-colnames(data1)
         len1<-length(col1)
       #  cat("len1 is:",len1,"\n")
         data1[,t]=data[,t]
         target=0
         p_target=1
         alle_target=0
         residual=0
         data2=data
         if(bi>0)
          {
            cat("I am B\n")
            cat("left is: ",left,"\n")
            for(name in left)
             {
               string=paste(col1[t],"~",paste(name,collapse="+"),sep="")
               cat("string is:",string,"\n")
               fit<-lm(data=data2,formula = as.formula(string),na.action=na.exclude)
               data2[,t]=residuals(fit)
               cat("Am I Ok\n")
             }
            data1[,t]=data2[,t]
          }
         cat("Conditioning is over\n")
         for(i in (j+1):len1)
          {
            string=paste(col1[t],"~",col1[i],sep="")
         #   cat("string is: ",string,"\n")
            fit<-lm(data=data1,formula = as.formula(paste(col1[t],"~",col1[i],sep="")),na.action=na.exclude)
            num<-length(summary(fit)$coefficients[,1])
            for(q in 2:num)
             {
          #     cat("I am g\n")
               p=summary(fit)$coefficients[q,4]
               alle=summary(fit)$coefficients[q,1]
               snp_test=col1[i]
           #    cat("snp is: ",snp_test," and p_value is ",p," allele effect is ",alle,"alle_target is: ",alle_target,"\n")
           #    if(p<p_target)
                if(p<p_target)
                 {
                   p_target=p
                   target=i
                   snp_test=col1[i]
           #        cat("change SNP to ",snp_test,"\n")
                 }          
             }
          }
          cat("p_target is:",p_target,"\n")
        if(p_target>cutoff)
         {
           break()
         }
        if(p_target<=cutoff)
          {
             cat("bi is:",bi,"\n")
             if(bi>0)
              {
         #       cat("Get in\n")
                y=0
         #       cat("left is:",left, " and col1[target] is: ",col1[target],"\n")
                left_test=c(left,col1[target])
                for(z in 1:length(left_test))
                 {
                   col_test<-left_test[-z]
           #        cat("calculate VIF\n")
                   string=paste(left_test[z],"~",paste(col_test,collapse="+"),sep="")
            #       cat("string is:",string,"\n")
                   vif<-VIF(lm(data=data,as.formula(paste(left_test[z],"~",paste(col_test,collapse="+"),sep="")),na.action=na.exclude))
             #      cat("is it Ok\n")
                   if(vif>=10)
                    {
                      y=1
                #      cat("Bad, VIF is too large\n")
                    }
              #     cat("Run next one\n")
                 }
                if(y==1)
                 {
               #    cat("I am n\n")
                   data1<-data1[,-target]
                   next
                 }
               # cat("Here am I\n")
                test=c(left,col1[target])
                x=0
             #   cat("extract effective factor\n")
                fit<-lm(data=data,formula = as.formula(paste(col1[t],"~",paste(test,collapse="+"),sep="")),na.action=na.exclude)
                for(i in 2:length(rownames(summary(fit)$coefficients)))
                 {
              #     cat("I am e\n")
               #    cat("i is:",i,"\n")
                   p=(summary(fit)$coefficients)[i,4]
                   variable=rownames(summary(fit)$coefficients)[i]
                #   cat("p is:", p, "and variable is: ",variable,"\n")
                   if(variable == col1[target]&& p<=cutoff)
                    {
                      x=1
                    }
                 }
               if(x==1)
                { 
                 # cat("I am z\n")
                  data1<-data1[,-target]
                  bi=bi+1
                  left[bi]=col1[target]
                  cat("left is:",left,"\n")
                  len=len-1
                  if(len<=0)
                   {
                     break()
                   }
                }
             cat("second run\n")              
              } else
              {
                cat("first run\n")
                data1<-data1[,-target]
                bi=bi+1
                left[bi]=col1[target]
                len=len-1
                cat("left is:",left,"\n")
               # cat("first run finish\n") 
              }                     
           }
          cat("Why there is problem with my scripts\n")
        }           
   cat("why is here\n")
   cat("left is:",left,"\n")
   if(is.na(left))
    {
      return(0)
    }
   fit<-lm(data=data,formula = as.formula(paste(col1[t],"~",paste(left,collapse="+"),sep="")),na.action=na.exclude)
   out_text=""
   num=length(summary(fit)$coefficients[,1])
   for(i in 2:num)
     {
       variable=rownames(summary(fit)$coefficients)[i]
       p_value=(summary(fit)$coefficients)[i,4]
       effect =(summary(fit)$coefficients)[i,1]
       if(p_value<=cutoff)
         {
           text<-c(col[t]," ",variable," ",p_value," ",effect)
           out_text<-rbind(out_text,text)
         }
     }
   }
   write.table(out_text,file=output,quote = FALSE, sep = "\t", eol = "\n", na = "0",dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))  
 }

Args <- commandArgs()

cat("Useage: R --vanilla --slave --input input_file --output output_file  < Yang_condition_analysis.R\n")
for (i in 1:length(Args))
 {
   if (Args[i] == "--input") file_name = Args[i+1]
   if (Args[i] == "--output") out_file = Args[i+1]
 }
my_conditional(file_name,out_file)




 
