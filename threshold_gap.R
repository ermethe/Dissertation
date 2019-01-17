#----------------FILE DESCRIPTION-------------------------------
#---------------------------------------------------------------
#This program computes the threshold gap of a proper sequence
#The definition cam from "Threshold Sequences" by Hammer,
#Ibaraki, and Simeone
#---------------------------------------------------------------

#---------Necessary Functions-----------------------------------

#check EG status
erdos.gallai = function(deg.vals,deg.cts){
  # checks if degree sequence is graphic
  # if deg.cts  = corresponding counts for each value in deg.vals
  #   deg.vals = unique degree values, in descending order
  N = sum(deg.cts)
  n.vals = length(deg.vals)
  check.cts = cumsum(deg.cts)
  
  j = 1
  H = R = 0
  L = n.vals
  
  # check Erdos-Gallai condition
  while ((deg.vals[j]>=check.cts[j])&(H>=0)){
    Q = deg.cts[j]*R
    while (deg.vals[L]<check.cts[j]){
      Q = Q+(check.cts[j]-deg.vals[L])*deg.cts[L]
      R = R+deg.cts[L]
      L = L-1
    } # end while L
    H = H+deg.cts[j]*(N-1-deg.vals[j])-Q
    j = j+1
  } # end while j
  EG.cond = (H>=0)
  
  # check if total degree is even
  odd.vals = round(round(deg.vals,-1)-deg.vals)%in%c(-5,-3,-1,1,3,5)
  odd.cts = round(round(deg.cts,-1)-deg.cts)%in%c(-5,-3,-1,1,3,5)
  odd.both = sum(odd.vals&odd.cts)
  even.sum = (round(odd.both,-1)-odd.both)%in%c(-4,-2,0,2,4)
  
  graphic = even.sum&EG.cond
  results = list(graphic=graphic,even.sum=even.sum,EG.cond=EG.cond)
  return(results)
}

#Compute EGVs
erdos.gallai.value <- function(degseq,k){
  n = length(degseq)
  if(k==1){
    y = degseq[1] - (n-1)
    return(y)
  }else if(k==n){
    y = sum(degseq) - n*(n-1)
    return(y)
  }else{
    t1 = sum(degseq[1:k])
    t2 = k*(k-1)
    indices = which(degseq[(k+1):n]>k)
    t3 = degseq[(k+1):n]
    t3[indices] = k
    y = t1-t2-sum(t3)
  }
  
}


#----------Program----------------------------------------------

delta.calc <- function(egvs){
  m = length(egvs)
  deltas = rep(0,m)
  deltas[1] = - egvs[1]
  for(i in 2:m){
    deltas[i] = egvs[(i-1)] - egvs[i]
  }
  return(m)
}

threshold.gap <- function(degs){
  n = length(degs)
  degs = sort(degs, decreasing = TRUE)
  if(max(degs)<= (n-1) & min(degs)>= 0){  #Ensure that the sequence is proper
    #Find istar
    compvals = seq(0,n-1)
    istar = max(which(degs<compvals))
    kvals = 1:istar
    EGVs = sapply(kvals,erdos.gallai.value,degseq=degs)
    delta.values = delta.calc(EGVs)
    threshold.gap = sum(abs(delta.values))/2
    output = list(threshold.gap = threshold.gap)
    return(output)
  }
  
  
}