



evolution_zg <- function(pop, pop.obj, obj.value, 
                      pop.size = pop.size, pheno.data,
                      select.rate = select.rate, add.rate = add.rate, 
                      mut.rate = mut.rate){
  
  ##select
  selected.individ <- obj.value[1:round((nrow(pop) * select.rate)),1]
  deleted.individ <- obj.value[-(1:round((nrow(pop) * select.rate))),1]
  
  added.individ <- deleted.individ[runif(length(deleted.individ)) < add.rate]
  selected.individ <- c(selected.individ,added.individ)
  
  selected.pop <- pop.obj[pop.obj$indival %in% selected.individ,]
  selected.obj <- obj.value[obj.value$indival %in% selected.individ,]
  
  ##crossover
  desired.len <- pop.size
  
  children.pop <- vector()
  
  for (i in 1:(nrow(selected.pop)-1)){
    individ.1 <- as.numeric(selected.pop[i,1:(ncol(selected.pop)-2)])
    individ.2 <- as.numeric(selected.pop[i+1,1:(ncol(selected.pop)-2)])
    
    cross.point <- round(length(individ.1)*runif(1,min = 0.1,max = 0.9))
    
    child.1 <- c(individ.1[1:cross.point],individ.2[-(1:cross.point)])
    child.2 <- c(individ.2[1:cross.point],individ.1[-(1:cross.point)])
    child2 <- rbind(child.1,child.2)
    
    ##mutation
    random.rate <- runif(nrow(child2))
    if(length(random.rate[random.rate < mut.rate])>0){
      mut.individ <- which(random.rate<mut.rate)
      for(m in mut.individ){
        pos.mutate <- sample(1:ncol(child2),1)
        child2[m,pos.mutate] <- abs(child2[m,pos.mutate] - 1)
      }
    }
    
    gene.num <- apply(child2,1,sum)
    child2 <- child2[gene.num > 0 & gene.num <= nrow(pheno.data)*0.5,]
    
    children.pop <- rbind(children.pop,child2)
  }
  
  children.pop.len <- nrow(children.pop) 
  
  while(children.pop.len < desired.len){
    selected.obj$prob <- exp(selected.obj$obj.fit)/sum(exp(selected.obj$obj.fit))
    individ.1 <- sample(selected.obj$indival[-1],1,prob = selected.obj$prob[-1])
    individ.2 <- selected.obj$indival[which(selected.obj$indival==individ.1)-1]
    
    individ.1 <- as.numeric(selected.pop[selected.pop$indival==individ.1,
                                         1:(ncol(selected.pop)-2)])
    individ.2 <- as.numeric(selected.pop[selected.pop$indival==individ.2,
                                         1:(ncol(selected.pop)-2)])
    
    cross.point <- round(length(individ.1)*runif(1,min = 0.1,max = 0.9))
    child.1 <- c(individ.1[1:cross.point],individ.2[-(1:cross.point)])
    child.2 <- c(individ.2[1:cross.point],individ.1[-(1:cross.point)])
    child2<-rbind(child.1,child.2) 
    
    ##mutation
    random.rate <- runif(nrow(child2))
    if(length(random.rate[random.rate < mut.rate])>0){
      mut.individ <- which(random.rate<mut.rate)
      for(m in mut.individ){
        pos.mutate <- sample(1:ncol(child2),1)
        child2[m,pos.mutate] <- abs(child2[m,pos.mutate] - 1)
      }
    }
    
    ## new generateion
    gene.num <- apply(child2,1,sum)
    child2 <- child2[gene.num > 0 & gene.num <= nrow(pheno.data)*0.5,]
    
    children.pop <- rbind(children.pop,child2)
    
    children.pop.len <- nrow(children.pop)  
  }
  
  return(children.pop)
}


