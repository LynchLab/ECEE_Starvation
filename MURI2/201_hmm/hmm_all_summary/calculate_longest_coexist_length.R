rm(list = ls())

library(plyr)  #round_any

summary <- read.table('/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/summary.txt', header=T)

outDF <- data.frame(pop=summary$pop, num_TP=summary$num_TP, pop_type=summary$pop_type,out=summary$out, max_coexist_length=NA, begin=NA, end=NA)

for (i in 1:120){
  pop <- summary$pop[i]
  mc_len <- 0
  myBegin <- 0
  myEnd <- 0
  if(pop!=208 & pop!=224){
    pop_100 <- round(pop/100)
    print(pop)
    if(summary$out[i]!=-1){
      myHaplotypeDF <- read.table(paste0('/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/',pop,'_haplotype_timecourse.txt'), sep=",")
      nTP <- ncol(myHaplotypeDF)
      myGen <- as.numeric(myHaplotypeDF[1,1:nTP])
      myMajorCladeFreq <- as.numeric(myHaplotypeDF[2,1:nTP])
      myMinorCladeFreq <- as.numeric(myHaplotypeDF[3,1:nTP])
      myDay <- round_any(myGen/((pop_100==1)*3.3+(pop_100==2)*3.3*4+(pop_100==3)*3.3*7+(pop_100==4)*0.33+(pop_100==5)*0.25),10)

      strong_coexist_property <- (myMajorCladeFreq>0.20) & (myMajorCladeFreq<0.80) & (myMinorCladeFreq>0.20) & (myMinorCladeFreq<0.80) #time point evident for both clade 
      
      if(sum(strong_coexist_property)>0){
        weak_coexist_property <- ((myMajorCladeFreq>0.01) & (myMajorCladeFreq<0.99)) | ((myMinorCladeFreq>0.01) & (myMinorCladeFreq<0.99)) #time point evident for at least one clade 
     
        #test data
        #myGen <- c(0,297,660,990,1320,1650,2310,2640)
        #weak_coexist_property <- c(F,T,F,T,F,T,F,T,F,T)
        #weak_coexist_property <- c(F,T,T,T,F,T,T,F,F,T)
        #weak_coexist_property <- c(F,T,T,T,F,T,T,T,F,F)
        count_consecutive_weak_coexist_property <- rep(0,length(weak_coexist_property))
      
        for (t in 2:(nTP-1)){
          if(weak_coexist_property[t] == TRUE){
            for(t_end in (t+1):nTP){
              if(weak_coexist_property[t_end] == TRUE){
                count_consecutive_weak_coexist_property[t] <- (myGen[t_end]-myGen[t])
              }else{
                break
              }
            }
          }
        }
      
        mc_len <- max(count_consecutive_weak_coexist_property)
        myBegin <- myGen[which(count_consecutive_weak_coexist_property==mc_len)]
        myEnd <- myBegin+mc_len
      }
    }
    outDF$max_coexist_length[i] <- mc_len
    outDF$begin[i] <- myBegin
    outDF$end[i] <- myEnd
  }
}

write.table(outDF,"/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/calculate_longest_coexist_length.out")