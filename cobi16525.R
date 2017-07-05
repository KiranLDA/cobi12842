#########################################################################
#                                                                       #
#  16/09/2016 code written by Kiran DHanjal-Adams as supporting         #
#  information for "Setting conservation priorities for migratory       #
#  networks under uncertainty" (doi:10.1111/cobi.16525). More           #
#  specifically, this code uses satellite tracks of Eastern Curlew      #
#  (Numenius madagascariensis) published in Driscoll and Ueta (2002)    #
#  (doi:10.1046/j.1474-919X.2002.00081.x) to (i) set up a connectivity  #
#  matrix, (ii) use this matrix to calculate the maximum flow of        #
#  birds through the network and (iii) prioritise sites according to    #
#  maximum flow and maximum count.                                      #
#                                                                       #
#########################################################################


#Load the necessary packages
library(Matrix)
library(gurobi) # note that gurobi must be setup before this package will work. See http://www.gurobi.com/products/modeling-languages/r 
library(maptools)
library(rgdal)
library(fields) 
library(plyr)


### Set working directory
setwd("C:\\write\\your\\filepath\\here\\")
load('AppendixS2.RData')


###############################################################
# Create a network of stopover sites specific to the species  #
###############################################################


#read in Bamford data
EC_bamford=bamford

#Create dataset for north and south migration sites
North=bamford[order(bamford$Latdec),]
North=North[c((length(North[,1])-1),1:(length(North[,1])-2),length(North[,1])),]

South=bamford[order(bamford$Latdec,decreasing=TRUE),]
South=South[c(length(South[,1]),1:(length(South[,1])-2),(length(South[,1])-1)),]

#project North South
North_pts=SpatialPoints(coords=North[2:(length(North[,1])-1),3:2], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
South_pts=SpatialPoints(coords=South[2:(length(North[,1])-1),3:2], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))

### calculate distances
North_Dists=rdist.earth(North_pts@coords)
South_Dists=rdist.earth(South_pts@coords)

#name columns and rows
colnames(North_Dists)=rownames(North_Dists)=North[2:(length(North[,1])-1),1]
colnames(South_Dists)=rownames(South_Dists)=South[2:(length(South[,1])-1),1]

# add columns and rows for non-breeding (sitename 0) and breeding (sitename 999)
# so that they can be identified later as different from stopover sites
#breeding supersink node t'
add_row=matrix(0,1,dim(North_Dists)[2])
North_Dists=rbind(add_row,North_Dists,add_row)
add_col=matrix(0,dim(North_Dists)[1],1)
North_Dists=cbind(add_col,North_Dists,add_col)
colnames(North_Dists)[1]=0
colnames(North_Dists)[length(North_Dists[1,])]=999
rownames(North_Dists)[1]=0
rownames(North_Dists)[length(North_Dists[,1])]=999


# non-breeding supersource node s'
add_row=matrix(0,1,dim(South_Dists)[2])
South_Dists=rbind(add_row,South_Dists,add_row)
add_col=matrix(0,dim(South_Dists)[1],1)
South_Dists=cbind(add_col,South_Dists,add_col)
colnames(South_Dists)[1]=999
colnames(South_Dists)[length(South_Dists[1,])]=0
rownames(South_Dists)[1]=999
rownames(South_Dists)[length(South_Dists[,1])]=0

#replace with zeros on one side to make sure it is a directional network
end=length(colnames(South_Dists))
for (i in 1:end){
  North_Dists[i:end,i]=0
  South_Dists[i:end,i]=0
}

#Execute this bit to make North and South different
not_NM=sort(unique(c(which(North$NM==1),which(North$B==1),which(North$NB==1))))
North_Dists=North_Dists[not_NM,not_NM]

not_SM=sort(unique(c(which(South$SM==1),which(South$B==1),which(South$NB==1))))
South_Dists=South_Dists[not_SM,not_SM]

#Link Non-breeding and Breeding sites with connecntivity matrix
North_Dists[1,which(rownames(North_Dists) %in% North$Site.Code[(which(North$NB==1))])]=Inf
North_Dists[1,1]=0
North_Dists[length(colnames(North_Dists))-1,length(colnames(North_Dists))]=Inf

South_Dists[which(rownames(South_Dists) %in% South$Site.Code[(which(South$NB==1))]),length(colnames(South_Dists))]=Inf
South_Dists[length(colnames(South_Dists)),length(colnames(South_Dists))]=0
South_Dists[1,2]=Inf



###########################################################
### Determine the northness/southness probability distribution
### using cosine of the azimuth A[uv] = |cos(Phi[uv])| in ms
###########################################################


#North
North_probability=North_Dists

for (i in 2:length(North_Dists[1,])-1){
  for (j in 2:length(North_Dists[,1])-1){
    if (North_Dists[i,j]>0){
      xLon=bamford$Londec[bamford$Site.Code==as.integer(rownames(North_Dists)[i])]
      xLat=bamford$Latdec[bamford$Site.Code==as.integer(rownames(North_Dists)[i])]
      yLon=bamford$Londec[bamford$Site.Code==as.integer(colnames(North_Dists)[j])]
      yLat=bamford$Latdec[bamford$Site.Code==as.integer(colnames(North_Dists)[j])]
      z=gzAzimuth(c(xLon,xLat),c(yLon,yLat))
      North_probability[i,j]=abs(cos(z*pi/180))
    }
  }
}
North_probability[1,]=North_Dists[1,]
North_probability[North_probability==Inf]<-1
North_probability[North_probability==0]<-NA

#South
South_probability=South_Dists

for (i in 2:length(South_Dists[1,])-1){
  for (j in 2:length(South_Dists[,1])-1){
    if (South_Dists[i,j]>0){
      xLon=bamford$Londec[bamford$Site.Code==as.integer(rownames(South_Dists)[i])]
      xLat=bamford$Latdec[bamford$Site.Code==as.integer(rownames(South_Dists)[i])]
      yLon=bamford$Londec[bamford$Site.Code==as.integer(colnames(South_Dists)[j])]
      yLat=bamford$Latdec[bamford$Site.Code==as.integer(colnames(South_Dists)[j])]
      z=gzAzimuth(c(xLon,xLat),c(yLon,yLat))
      South_probability[i,j]=abs(cos(z*pi/180))
    }
  }
}
South_probability[1,]=South_Dists[1,]
South_probability[South_probability==Inf]<-1
South_probability[South_probability==0]<-NA


###########################################################
### Determine the distance probability distribution
### i.e. P[uv] in ms
###########################################################


#using the density function to fit the probability distribution to tracks
probN=approxfun(density(trax$distance[trax$species=="Eastern Curlew" & trax$NorS=="N"],adjust=2))
probS=approxfun(density(trax$distance[trax$species=="Eastern Curlew" & trax$NorS=="S"],adjust=2))

# apply the distance probability distribution to distance data between the pairs of nodes in Bamford
N_dist_p=North_Dists
S_dist_p=South_Dists

N_dist_p[N_dist_p==0]<-NA
S_dist_p[S_dist_p==0]<-NA

N_dist_p[N_dist_p==Inf]<-5000
S_dist_p[S_dist_p==Inf]<-5000

N_dist_p=apply(N_dist_p,c(1,2),function(N_dist_p) probN(N_dist_p))
S_dist_p=apply(S_dist_p,c(1,2),function(S_dist_p) probS(S_dist_p))

N_dist_p[which(N_dist_p[,length(N_dist_p[,1])]>0), length(N_dist_p[,1])]=1
N_dist_p[1,which(N_dist_p[1,]>0)]=1

S_dist_p[which(S_dist_p[,length(S_dist_p[,1])]>0), length(S_dist_p[,1])]=1
S_dist_p[1,which(S_dist_p[1,]>0)]=1


###########################################################
### Dermine the bird proportion
### i.e. N[v] in manuscript
###########################################################

#north
N_prop_max_K=North_probability
bamford_North=North$Max.Count[not_NM]

for (i in 2:(length(colnames(N_prop_max_K))-2)){
  N_prop_max_K[i,i:(length(colnames(N_prop_max_K))-2)]=bamford_North[i:(length(colnames(N_prop_max_K))-2)]/max(bamford_North[2:(length(colnames(N_prop_max_K))-2)])
}

for (i in 2:(length(colnames(N_prop_max_K))-2)){
  if (!is.na(N_prop_max_K[1,i])){
    N_prop_max_K[1,i]=bamford_North[i]/max(bamford_North[2:(length(colnames(N_prop_max_K))-2)])
  }
}
N_prop_max_K[2:(length(colnames(N_prop_max_K))-2),length(colnames(N_prop_max_K))-1]=1

#South
S_prop_max_K=South_probability
bamford_South=South$Max.Count[not_SM]

for (i in 3:(length(colnames(S_prop_max_K))-1)){
  S_prop_max_K[i,i:(length(colnames(S_prop_max_K))-1)]=bamford_South[i:(length(colnames(S_prop_max_K))-1)]/max(bamford_South[2:(length(colnames(S_prop_max_K))-1)])
}
S_prop_max_K[(length(colnames(S_prop_max_K))-1),(length(colnames(S_prop_max_K))-1)]=NA


#########################################################################
### calculate the overall probability of getting from one site to the next
### i.e. w[uv] in manuscript
#########################################################################

EC_N_overall_probability=North_probability*N_dist_p*N_prop_max_K
EC_S_overall_probability=South_probability*S_dist_p*S_prop_max_K

####################################################################################
### calculate numbers going from one bamford site to the other
### using c[uv]  =  x[u] ( w[uv] / (Sum {v: (u, v) element of E} w[uv] )) in ms
####################################################################################

##############
#### North
bamford=EC_bamford

###calculate the overall probability of getting from one site to the next
overall_probability=EC_N_overall_probability
proportion<-apply(overall_probability,1,function(overall_probability) overall_probability[which(!is.na(overall_probability))]/sum(overall_probability,na.rm=TRUE))

#Calculate the number of birds getting from one site to the next (in=out)
pop=bamford$Max.Count[1]
proportions=numbers=overall_probability
for (i in 1:length(proportion)){
  for (j in 1:length(proportion[[i]])){
    proportions[rownames(proportions)==names(proportion[i]),colnames(proportions)==names(proportion[[i]][j])]=proportion[[i]][j]    
  }
}
numbers[1,]=pop*proportions[1,]
for (i in 2:length(proportion)){
  numbers[i,]=sum(numbers[,i],na.rm=TRUE)*proportions[i,]
}

EC_N_numbers=numbers


#############
### South
bamford=EC_bamford

###calculate the overall probability of getting from one site to the next
overall_probability=EC_S_overall_probability
proportion<-apply(overall_probability,1,function(overall_probability) overall_probability[which(!is.na(overall_probability))]/sum(overall_probability,na.rm=TRUE))


#Calculate the number of birds getting from one site to the next (in=out)
#pop=bamford$Max.Count[1]
proportions=numbers=overall_probability
for (i in 1:length(proportion)){
  for (j in 1:length(proportion[[i]])){
    proportions[rownames(proportions)==names(proportion[i]),colnames(proportions)==names(proportion[[i]][j])]=proportion[[i]][j]    
  }
}
numbers[1,]=pop*proportions[1,]
for (i in 2:length(proportion)){
  numbers[i,]=sum(numbers[,i],na.rm=TRUE)*proportions[i,]
}

EC_S_numbers=numbers

##############################################################
### Combine N and S to make huge proportions matrix of birds going between pairs of sites
############################################################

EC_N_numbers=EC_N_numbers[-length(rownames(EC_N_numbers)),-1]
EC_S_numbers=EC_S_numbers[,-1]

topright=matrix(NA,length(rownames(EC_N_numbers)),length(colnames(EC_S_numbers)))                          
bottomleft=matrix(NA,length(rownames(EC_S_numbers)),length(colnames(EC_N_numbers))) 

top = cbind(EC_N_numbers,topright)
bottom = cbind(bottomleft,EC_S_numbers)
          
EC_matrix=rbind(top,bottom)

colnames(EC_matrix)=c(paste0("N",colnames(EC_N_numbers)),paste0("S",colnames(EC_S_numbers)))
rownames(EC_matrix)=c(paste0("N",rownames(EC_N_numbers)),paste0("S",rownames(EC_S_numbers)))

rownames(EC_matrix)[which(rownames(EC_matrix) %in% "S999")]="N999"
colnames(EC_matrix)[which(colnames(EC_matrix) %in% "S999")]="N999"



#################################
##   GUROBI OPTIMISATION
##   prioritise flow
#################################

#initialise matrix with numbers of birds going to each site from each other site. Note these are not round
bird=apply(EC_matrix,1,function(EC_matrix) EC_matrix[which(!is.na(EC_matrix))])


### set up variables for gurobi


a=0 #initialise the upper bounds matrix 
lp_matrix<-matrix(0,length(rownames(EC_matrix)),1) # initialise connectivity matrix
rownames(lp_matrix)=rownames(EC_matrix) #naming rows and columns
counter=1 #used to track iterations

for (i in 1:(length(rownames(EC_matrix))-1)){
  for (j in 1:length(bird[[i]])){
    a <-c(a, bird[[i]][j]) #basically puts bird data in a format that gurobi can deal with
    lp_matrix[rownames(lp_matrix)==names(bird[i]),counter]<-1 # put a 1 if going from-to, otheriwse 0
    lp_matrix[rownames(lp_matrix)==names(bird[[i]][j]),counter]<--1
    counter=counter+1
    lp_matrix<-cbind(lp_matrix,matrix(0,length(rownames(EC_matrix)),1))
  }
}
a<-a[2:length(a)]
a<-c(a,Inf)
lp_matrix[1,length(lp_matrix[1,])]<--1 #create connection with supersource node
lp_matrix[length(rownames(EC_matrix)),length(lp_matrix[1,])]<-1 #create connection with supersing node
I = lp_matrix #save connectivity matrix
Up_bounds = a # save the upper bounds
Cost=matrix(-1,1,length(a)) #cost is a dummy variable set to -1 to generate profit

rownames(I)=rownames(EC_matrix)
colnames(I)=names(Up_bounds)

#sites at the moment are named N and S according to when birds use them, but we are removing sites during the prioritisation regardless of  whether they are used during N or S migration, so we remove these from the name
site_rows=matrix(rownames(I),1,length(rownames(I)))
for (i in 1:length(site_rows)){
  site_rows[i]=strsplit(site_rows[i],"N|S")[[1]][2]
}
rownames(I)=site_rows

site_cols=matrix(colnames(I),1,length(colnames(I)))
for (i in 1:length(site_cols)){
  site_cols[i]=strsplit(site_cols[i],"N|S")[[1]][2]
}

#continue naming everything correctly
colnames(I)=site_cols
colnames(I)[length(I[1,])]="0"
names(Up_bounds)=colnames(I)


#set up the matrix that will hold the results
prioritisation=matrix(NA,1,6)
colnames(prioritisation)=c("site","Lon","Lat","flow","K","count")
count=0

#initialise loops and save I data as before I as it will get re-written in prioritisation
before_I=I
before_Up_bounds=Up_bounds

site_numbers=unique(as.vector(rownames(before_I)))
never=c("0","999","9997") #these sites are never to be removed
a=which(site_numbers %in% never) #indexes which sites should not be removed 
site_numbers=site_numbers[-a] #removes these sites from the possible site list used in prioritisation

#start loop - when there are no longer any sites, the loop ends
while (length(site_numbers)>0){
  count=count+1 #keeps track
  flow=matrix(NA,2,2)
  colnames(flow)=c("site","flow")
  
  #Remove each site in sequence to see how much to contributes to flow, find the site that contributes least, and removed it from the network, then repeat
  for (R in 1:length(site_numbers)){ # loop through all the sites that can be removed
    I=before_I
    Up_bounds=before_Up_bounds
    
    seq_removal=site_numbers[R]
    
    #index where site is to remove
    rows=NA
    columns=NA
    rows = c(rows,which(rownames(I) %in% seq_removal))
    for (j in 1:length(rows)){
      columns = c(columns,which(I[rows[j],]!=0))
    }
    rows=unique(rows[2:length(rows)])
    columns=unique(columns[2:length(columns)])
    #remove the site
    I=I[-rows,-columns]
    Up_bounds=Up_bounds[-columns]
    
    #do the max flow
    model <- list()
    model$rhs <- as.vector(matrix(0,1,length(rownames(I))))# right hand side vector for linear constraint
    model$sense  <- as.vector(matrix('=',length(rownames(I)),1)) # says that we want the same number coming in and going out of a node
    model$ub <-Up_bounds # upper bound vector - sets an upper limit to the number of birds migrating along an edge
    model$lb <-as.vector(matrix(0,1,length(Up_bounds)))# lower bound vector i.e. can't have less than 0 birds migrating along an edge
    model$A <-as.matrix(I,ncol = ncol(dd), dimnames = NULL) #linear constraint matrix
    model$obj <- matrix(-1,1,length(Up_bounds)) #linear object vector, must specify 1 value for each row of A
    model$modelsense <- "min"
    result <- gurobi(model)
    print(result$objval)
    #print(result$x)
    print(paste0("the optimal solution: ",result$x[length(result$x)]," birds"))
    
    flow=rbind(flow,c(seq_removal,result$x[length(result$x)]))
  }
  
  flow=flow[3:length(flow[,1]),]
  flow=flow[order(flow[,2],decreasing=TRUE),]
  
  # once we have figured out each site contribution to flow, we permanently remove the one that contributes least 
  to_remove=flow[,1][flow[,2]==max(flow[,2])]
  to_remove_flow=flow[,2][flow[,2]==max(flow[,2])]
  
  print(to_remove)
  
  #reset all the variables with the site now deleted from list of available sites
  rows=NA
  columns=NA
  for (i in 1:length(to_remove)){
    rows = c(rows,which(rownames(before_I) %in% to_remove[i]))
    for (j in 1:length(rows)){
      columns = c(columns,which(before_I[rows[j],]!=0))
    }
  }
  rows=unique(rows[2:length(rows)])
  columns=unique(columns[2:length(columns)])
  before_I=before_I[-rows,-columns]
  before_Up_bounds=before_Up_bounds[-columns]
  
  site_numbers=unique(as.vector(rownames(before_I)))
  never=c("0","999","9997")
  a=which(site_numbers %in% never)
  site_numbers=site_numbers[-a]
  
  #save data
  end=matrix(NA,length(to_remove),6)
  for (i in 1:length(to_remove)){
    end[i,1]=to_remove[i]
    end[i,2]=bamford$Lon[bamford$Site.Code==to_remove[i]]
    end[i,3]=bamford$Lat[bamford$Site.Code==to_remove[i]]
    end[i,4]=to_remove_flow[i]
    end[i,5]=bamford$Max.Count[bamford$Site.Code==to_remove[i]]
    end[i,6]=count
  }
  prioritisation=rbind(prioritisation,end)
  
}
end=matrix(NA,1,6)
end[1]=flow[1]
end[2]=bamford$Lon[bamford$Site.Code==flow[1]]
end[3]=bamford$Lat[bamford$Site.Code==flow[1]]
end[4]=flow[2]
end[5]=bamford$Max.Count[bamford$Site.Code==flow[1]]
end[6]=count+1

prioritisation=rbind(prioritisation,end)




### Save the flow prioritsation data
flow_prioritisation=prioritisation




#################################
##   GUROBI OPTIMISATION 
##   Does the same but prioritises sites with big bird counts
#################################

bird=apply(EC_matrix,1,function(EC_matrix) EC_matrix[which(!is.na(EC_matrix))])#/1,na.rm=TRUE)


#set up variables for gurobi
a=0
lp_matrix<-matrix(0,length(rownames(EC_matrix)),1)
rownames(lp_matrix)=rownames(EC_matrix)
counter=1
for (i in 1:(length(rownames(EC_matrix))-1)){
  for (j in 1:length(bird[[i]])){
    a<-c(a, bird[[i]][j])
    lp_matrix[rownames(lp_matrix)==names(bird[i]),counter]<-1
    lp_matrix[rownames(lp_matrix)==names(bird[[i]][j]),counter]<--1
    counter=counter+1
    lp_matrix<-cbind(lp_matrix,matrix(0,length(rownames(EC_matrix)),1))
  }
}
a<-a[2:length(a)]
a<-c(a,Inf)
lp_matrix[1,length(lp_matrix[1,])]<--1
lp_matrix[length(rownames(EC_matrix)),length(lp_matrix[1,])]<-1
I= lp_matrix
Up_bounds=a
Cost=matrix(-1,1,length(a)) #cost of dummy variable set to -1 to generate profit

rownames(I)=rownames(EC_matrix)
colnames(I)=names(Up_bounds)

site_rows=matrix(rownames(I),1,length(rownames(I)))
for (i in 1:length(site_rows)){
  site_rows[i]=strsplit(site_rows[i],"N|S")[[1]][2]
}
rownames(I)=site_rows

site_cols=matrix(colnames(I),1,length(colnames(I)))
for (i in 1:length(site_cols)){
  site_cols[i]=strsplit(site_cols[i],"N|S")[[1]][2]
}
colnames(I)=site_cols
colnames(I)[length(I[1,])]="0"
names(Up_bounds)=colnames(I)


#order the sites accoring to count so they can be removed in this order in the for loop
sites=bamford$Site.Code[order(bamford$Max.Count)]
never=c("0","999","9997")
a=which(sites %in% never)
sites=sites[-a]


prioritisation=matrix(NA,1,6)
colnames(prioritisation)=c("site","Lon","Lat","flow","K","count")
count=0


for (iter in 1:length(sites)){
  to_remove=sites[iter]
  
  rows=NA
  columns=NA
  for (i in 1:length(to_remove)){
    rows = c(rows,which(rownames(I) %in% to_remove[i]))
    for (j in 1:length(rows)){
      columns = c(columns,which(I[rows[j],]!=0))
    }
  }
  rows=unique(rows[2:length(rows)])
  columns=unique(columns[2:length(columns)])
  I=I[-rows,-columns]
  Up_bounds=Up_bounds[-columns]
  
  model <- list()
  model$rhs <- as.vector(matrix(0,1,length(rownames(I))))#c(as.vector(matrix(0,1,length(rownames(matrix)))),CarryingCapacity) #right hand side vector for linear constraint
  model$sense  <- as.vector(matrix('=',length(rownames(I)),1))#c(as.vector(matrix('=',length(rownames(matrix)),1)),as.vector(matrix('<=',length(rownames(matrix)),1)))
  model$ub <-Up_bounds #max #upper bound vector
  model$lb <-as.vector(matrix(0,1,length(Up_bounds)))#lower bound vector
  model$A <-as.matrix(I,ncol = ncol(dd), dimnames = NULL) #linear constraint matrix
  model$obj <- matrix(-1,1,length(Up_bounds)) #linear object vector, must specify 1 value for each row of A
  model$modelsense <- "min"
  result <- gurobi(model)
  print(result$objval)
  print(paste0("the optimal solution: ",result$x[length(result$x)]," birds"))
  
  # put all the data in this matrix called end
  end=matrix(NA,length(to_remove),6)
  for (i in 1:length(to_remove)){
    end[i,1]=to_remove[i]
    end[i,2]=bamford$Lon[bamford$Site.Code==to_remove[i]]
    end[i,3]=bamford$Lat[bamford$Site.Code==to_remove[i]]
    end[i,4]=result$x[length(result$x)]
    end[i,5]=bamford$Max.Count[bamford$Site.Code==to_remove[i]]
    end[i,6]=count
  }
  #add this to the priorisation
  prioritisation=rbind(prioritisation,end)
}

plot(prioritisation[,6],prioritisation[,4])
plot(prioritisation[,4],type="l")



### Save the data
N_prioritisation=prioritisation


#################################
##   PLOT  EVERYTHING
#################################


plot(flow_prioritisation[,4],pch=17,type="o",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, 
     xlab='Number of sites removed',
     ylab='Population migrating through network')
lines(N_prioritisation[,4],pch=16,lty=2, type="o")
legend('bottomleft', 
       legend=c("Prioritise flow",
                "Prioritise count"), 
       lty=c(1,2), pch=c(17,16),bty='n')
