

set.seed(10)

library(dismo)			### Library dismo will be needed for the sensitivity analysis using BRTs

### Here starts the fun! Metapopulation model function. Arguments: n.sim: number of simulations to undertake; dim.sq: vector containing the dimensions of one of the sides of the square defining the landscape; 
### n.year: number of consecutive years simulated 

metapop.model<-function(n.sim, dim.sq, n.year){
	### Load the required libraries: fields for the pairwise euclidean distance and lhs for the latin hypercube sampling
	### igraph for networks
	require(fields)
	require(lhs)
  require(igraph)
  require(qgraph)
  
  #### structure for holding the results
	outsim <- vector("list") 
	
	### Latin hypercube sampling of some parameters
  lh <- randomLHS(n.sim, 4)									### Generate values
	input <- matrix(0, nrow=n.sim, ncol=4)
	colnames(input) <- c("Number of patches", "Alpha", "c", "z")
	
	### Actual parameters values
	input[,1] <- round(qunif(lh[,1], min=2, max=100), digits=0)	### Number of patches; PGD: I've increased it to 100 patches maximum
	input[,2] <- round(qunif(lh[,2], min=0, max=900), digits=2)	### Density per km2 (alpha); from estimates in a report by Warburton et al. (2009).
	######## I converted these units to density per squared kilometer. Double check this (PGD: I've included zero, so some patches can be vacant)
	input[,3] <- round(qunif(lh[,3], min=1/5, max=1/1), digits=2)		### Exponential distance decay dispersal kernel in kilometres (from 1 to 5 km); based on estimates from Etherington et al. (2014): http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0088293
	input[,4] <- round(qunif(lh[,4], min=0.1, max=1), digits=2)		### Power law exponent (for the abundance-area relationship)
			
	### Start simulations
	
	for (s in 1:n.sim){			
    pars <- vector("list") 					### Create empty list to store results
		
    A <- pars$Area<- round(runif(input[s,1], 0.01, 1), digits=2)		
    
    ### PGD: Patch area in hectares drawn from a Negative Binomial distribution and transformed into km2 (not fully convinced of this)
	  ### ML: I agree, we need a more heterogeneous distribution for area size. Can you implement this?
    ### Also, are we sure we cannot have bigger areas?
    
    r <- runif(input[s,1], min=.1, max=3)					### Per capita population growth rate (due to births and death, but not because of migration - migration is modelled directly in raw numbers)
		
		##### I have increased the max r to 3 to allow for richer dynamics in the logistic equation; PGD: yeah, good. 

		### Create the patch centroids and calculate the distance matrix
		x.sim <- sample(dim.sq, input[s,1])					### Simulate x-coordinates for the centroids using sample without replacement	
		y.sim <- sample(dim.sq, input[s,1])					### Simulate y-coordinates for the centroids using sample without replacement	
		coord.patches <- matrix(c(x.sim, y.sim), nrow=input[s,1], ncol=2)		### Put together the patch centroids
		D <- rdist(coord.patches)							### Pairwise euclidean distance matrix calculated using function rdist from package fields
		diag(D) <- NA								### Diagonal values are set to NA - trick to avoid problems with computations and calculations

		##### Initial population size and k (carrying capacity) for each patch

		N <- matrix(0, ncol=input[s,1], nrow=n.year)			### Create empty matrix to store abundance results
		N[1,] <- round(input[s,2]*A^input[s,4], digits=0)		### Define initial population size based on the power-law relationship
		
		##########
		k <- round(runif(input[s,1], min=min(N[1,]), max=max(N[1,])),digits=0)		### Define patch-level k as the inital population size plus a random number of additional animals
		##########  I have changed the way k is defined.... double check this! (Miguel); PGD: yeah, good

		#### Exponential distance-decay dispersal kernel
		dk <- exp(-input[s,3]*D)
		
		#### we normalise dk so the probabilities of moving sum up to one (to make sure all immigrants are distributed)
		### PGD: we don't need to normalise the dispersal kernel twice - I've removed the spatial normalisation in the previous line :)
		dk <- (dk/rowSums(dk,na.rm=T))
		
		### Individual probability of dispersal in each patch
    prob.disp <- round(runif(input[s,1], min=0.01, max=0.2), digits=2)	### Annual probability of dispersal (adult possums are recorded up to a maximum of 0.2); based on estimates from Ramsey & Efford (2010): http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2010.01839.x/abstract
		### If according to Ramsey et al, the maximum dispersal probability should be 0.2, why do you use 0.7 as the max?
    
    ### Changed according to the comment above; PGD: yeah, fine. 
		### Create an array to store the results of the movement patterns (a matrix each year)
		movement <- array(0, dim=c(nrow=input[s,1], ncol=input[s,1], n.year))
		
		#### Simulate the dynamics
		#### Here start the simulations
		
		for (i in 2:n.year){        ### n.year: number of years simulated
		  ##### here we calculate dispersers for each local population and store the results in the
		  ##### movement matrix
		  emigrants <- round((N[(i-1),] * prob.disp), digits=0)
		  current_disp <- round(dk*emigrants, digits=0)
		  
		  emigrants <- rowSums(current_disp, na.rm=T)
		  immigrants <- colSums(current_disp, na.rm=T)
		  
		  movement[,,i] <- current_disp
		  
		  #### We should think about whether dispersing individuals should be incorporated 
		  #### before or after growth (i.e., should be considered for growth in the current iteration)
		  #### if considered then they should be added also next to the k
		  
		  N[i,] <- N[i-1, ] + round( ((r * (N[(i-1),] - emigrants + immigrants ) ) * ((k - (N[(i-1),]- emigrants + immigrants ))/k)), digits=0)
		  N[i,][N[i,]< 0] <- 0
	
    }
		### PGD: Can we take the graphs out of the function? For me, it is easier to visualise and work with them
		
		matplot(N,type='l', lty=1, lwd=2, col=rainbow(14), main='Local pops dynamics', xlab='Time', ylab='Abundance')
		plot(rowSums(N), main='Metapop dynamics', xlab='Time', ylab='Abundance', type='l', lwd=2)
  
		# a possible way to define the network of connections between patches observed from the random process is to set a threshold
    # for example:
		conn_thrs <- mean(dk, na.rm=T)
		
		# ML: I changed the threshold to be the mean of dk, but I think we need a more clever
		# way of deciding the threshold
		
		# this is equivalent to say that no connection is assumed 
		# between patches for which the probability to arrive from 
		#one to the other is smaller than .1
		adj_m <- dk
		adj_m[adj_m > conn_thrs] <- 1
		adj_m[is.na(adj_m)] <- 0
		adj_m[adj_m < conn_thrs] <- 0
		
		### this is how the network of patches looks like
		g <- graph_from_adjacency_matrix(adj_m, mode='undirected')
		plot(g)
		
		# and then we can play around with this and get some stats like for example
		# obtain the number of neighbors of the largest patch (in terms of its area)
		# and maybe then check later if this has anything to do with the observed dynamics?
		
		largest_patch <- which(A == max(A))
		print(igraph::degree(g,largest_patch))
		
		#### The small world index
		print(as.numeric(qgraph::smallworldness(g)[1]))
		
		#### so we can se results a little bit...
		#### I print the dynamics for the first 5 populations
		matplot(N[,1:5], type='l', lwd=2)
		
		#### PARAMETERS TO OUTPUT
		
		##### maybe a better measure for the sensitivity analyses would be some sort of 
		##### weighted mean abundance? (weighted by the area) - so it is more than just 
		##### persistent and non-persistent
		### PGD: we need to change all of this below to adapt it to the updated model; just to reflect what are we doing :)
		
    n.patches <- as.numeric(input[s,1])			### Number of patches
  	pars$N<-N						### Matrix with population size per year in each patch
  	pars$N.end<-N[n.year,]					### Population size in each patch at the last year
  	pars$init<-N[1,]					### Initial population size in each patch
  	
  	### PGD: I've made changes to the parameters to output, so I just changed the whole batch of code (I'm too lazy)
    pars$k <- k						### Patch carrying capacity
  	pars$alpha <- rep(input[s,2], n.patches)			### Density, alpha (repeat the same value for each patch)
  	pars$c <- rep(input[s,3], n.patches)			### Exponential distance decay dispersal kernel (repeat the same value for each patch)
  	pars$z <- rep(input[s,4], n.patches)			### power law exponent (repeat the same value for each patch)
  	pars$prob.disp <- prob.disp				### Annual probability of dispersalfor each patch
  	pars$min.dist <- apply(D, 1, function(x){min(x, na.rm=T)})		### Distance to the nearest neighbour patch (in km)
  	
  	
  	pars$conn <- degree(g)					
  	
  	### Connectivty measured as the number of neighbours; PGD: check this one Miguel
  	### ML: Looks good to me.
  	
  	pars$r <- r						### Per capita growth rate in each patch
  	pars$prop.hab <- rep(sum(A)/(max(dim.sq)^2), input[s,1])		### Proportion of the landscape that is 'good' habitat
		outsim[[s]] <- pars
  }
  outsim
}



#this function incorporates allee effects to the model
metapop.nallee <- function(n.sim, dim.sq, n.year, cat.year, with.allee=F){
  ### Load the required libraries: fields for the pairwise euclidean distance and lhs for the latin hypercube sampling
  require(fields)
  require(lhs)
  require(igraph)
  require(qgraph)
  
  #### structure for holding the results
  outsim <- vector("list") 
  
  ### Latin hypercube sampling of some parameters
  lh <- randomLHS(n.sim, 5)									### Generate values
  input <- matrix(0, nrow=n.sim, ncol=5)
  colnames(input) <- c("Number of patches", "Alpha", "c", "z", "conn.trh")
  
  ### Actual parameters values
  input[,1] <- round(qunif(lh[,1], min=2, max=200), digits=0)				### Number of patches
  input[,2] <- qunif(lh[,2], min=0, max=1000)  						### Density per km2 (alpha); from estimates in a report by Warburton et al. (2009 - per hectare is better
  input[,3] <- qunif(lh[,3], min=1/5, max=1/1)					        ### Exponential distance decay dispersal kernel in metres (from 1 to 5 km); based on estimates from Etherington et al. (2014): http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0088293
  input[,4] <- qunif(lh[,4], min=0.1, max=1)			### Power law exponent (for the abundance-area relationship)
  input[,5] <- qunif(lh[,5], min=1, max=3)					### Connectivity threshold for delineating the networks (in kms)
  
  for (s in 1:n.sim){			### Start simulations
    pars <- vector("list") 								### Create empty list to store results
    A <- pars$Area <- round(runif(input[s,1], min=0.01, max=2), digits=2)			### Patch area varies from 1 to 100 ha
    cat.mag <- round(runif(input[s,1], min=0, max=1), digits=2)				### Magnitude of the catastrophe (proportion remaing after catastrophe)
    
    ### Create patch centroids and calculate the distance matrix
    x.sim <- sample(dim.sq, input[s,1])						### Simulate x-coordinates for the centroids using sample without replacement	
    y.sim <- sample(dim.sq, input[s,1])						### Simulate y-coordinates for the centroids using sample without replacement	
    coord.patches <- matrix(c(x.sim, y.sim), nrow=input[s,1], ncol=2)		### Put together the patch centroids
    D <- rdist(coord.patches)							### Pairwise euclidean distance matrix calculated using function rdist from package fields
    diag(D) <- NA									### Diagonal values are set to NA - trick to avoid problems with computations and calculations
    
    ##### Initial population size and k (carrying capacity) for each patch
    N <- matrix(0, ncol=input[s,1], nrow=n.year)					### Create empty matrix to store results
    N[1,] <- k <- round(input[s,2]*A^input[s,4], digits=0)				### Define initial population size based on the power-law relationship; initialise the local population at k
    
    ### for the allee effect
    if(with.allee) allee <- round((k*.1))
    
    #### Exponential distance-decay dispersal kernel 
    dk <- exp(-input[s,3]*D)
    
    #### we normalise dk so the probabilities of moving sum up to one (to make sure all immigrants are distributed) ### Don't need to normalise twice
    dk <- (dk/rowSums(dk,na.rm=T))
    
    ### Probability of dipsersal in each patch
    prob.disp <- round(runif(input[s,1], min=0.01, max=0.7), digits=2)	### Annual probability of dispersal
    
    ### Create an array to store the results of the movement patterns (a matrix each year)
    movement <- array(0, dim=c(nrow=input[s,1], ncol=input[s,1], n.year))
    
    ### Create a matrix to store the per capita growth rate per year and patch
    r <- matrix(0, ncol=input[s,1], nrow=n.year)
    
    #### Simulate the dynamics
    #### Here start the simulations for each + 1 year
    for (i in 2:n.year){        ### n.year: number of years simulated
      ### Per capita growth rate
      
      #### I don't understand this. why do you want to change r per year, and why round it?
      
      r[i,] <- round(runif(input[s,1], min=0.1, max=2))
      
      ### here we calculate dispersers for each local population and store the results in the
      ### movement matrix
      emigrants <- round((N[(i-1),] * prob.disp), digits=0)
      current_disp <- round(dk*emigrants, digits=0)
      
      emigrants <- rowSums(current_disp, na.rm=T)
      immigrants <- colSums(current_disp, na.rm=T)
      movement[,,i] <- current_disp
      
      #### Dispersal after growth
      # N[i,] <- N[i-1, ] + round( ((r[i,] * (N[(i-1),] - emigrants + immigrants ) ) * ((k - (N[(i-1),]- emigrants + immigrants ))/k)), digits=0)
      # N[i,][N[i,]< 0] <- 0
      
      #### ML: If dispersal has to be after growth, I would do it like this:
      if(with.allee) N[i,] <- N[i-1, ] + round( ((r[i,] * (N[(i-1),])) * ((k - (N[(i-1),] ))/k) * ((N[(i-1),]/allee) - 1 )), digits=0) - emigrants + immigrants
      else N[i,] <- N[i-1, ] + round( ((r[i,] * (N[(i-1),] ) ) * ((k - (N[(i-1),] ))/k) ), digits=0) - emigrants + immigrants
      N[i,][N[i,]< 0] <- 0
      
      
      #### ML: I re-implemented this bit to have it all in the same for loop,
      #### but not sure if you wanted catastrophes to happen every year.
      #### It didn't make sense to me to have it each iteration on the catastrophe year
      #### Your original code is still below...
      if(i <= cat.year){
        N[i, ] <- round(N[i, ]*cat.mag, digits=0)
      }
      
      # ### Popualtion in the year of the catastrophe - the catastrophe happens after growth and dispersal
      # N[cat.year, ] <- round(N[cat.year, ]*cat.mag, digits=0)
      # N[cat.year,][N[cat.year,]< 0] <- 0
    }
    
    # Define the network of connections between patches observed from the random process is to set a threshold 
    # this is equivalent to say that no connection is assumed 
    # between patches for which the probability to arrive from 
    #one to the other is smaller than the threshold
    
    adj_m <- D
    adj_m[adj_m < input[s,5]] <- 1
    adj_m[is.na(adj_m)] <- 0
    adj_m[adj_m > input[s,5]] <- 0
    
    ### this is how the network of patches looks like
    g <- graph_from_adjacency_matrix(adj_m, mode='undirected')
    plot(g)
    
    # and then we can play around with this and get some stats like for example
    # obtain the number of neighbors of the largest patch (in terms of its area)
    # and maybe then check later if this has anything to do with the observed dynamics?
    
    largest_patch <- which(A == max(A))
    print(paste('The largest patch has area ', A, 'and its degree is:', igraph::degree(g,largest_patch)))
#    print(igraph::degree(g,largest_patch))
    
    #### The small world index
    print(paste('small worldness is:', as.numeric(qgraph::smallworldness(g)[1])))
    
    #### so we can se results a little bit...
    matplot(N[,1:5], type='l', lwd=2)

    ### Parameters of interest: ann.bef: mean annual abundance before catastrophe
    #ann.bef<-apply(N[1:(cat.year-1), ], 2, mean)
    
    ### time until recovering to a threshold density of interest of 10/km2 (0.1/ha)
    t.recov<-rep(0, input[s,1])
    for (j in 1:input[s,1]){
      t.recov[j] <- which(N[,j]/A[j] >= 10)[1]					
    }
    
    ### OUTPUT VALUES
    n.patches<-pars$n.patches<-as.numeric(input[s,1])			### Number of patches
    # pars$ann.dens<-ann.bef							### Mean density pre-disturbance
    pars$t.recov<-t.recov							### Time to recovery
    pars$k<-k								### Patch carrying capacity
    pars$alpha<-rep(input[s,2], n.patches)					### Density, alpha (repeat the same value for each patch)
    pars$c<-rep(input[s,3], n.patches)					### Exponential distance decay dispersal kernel in metres (repeat the same value for each patch)
    pars$z<-rep(input[s,4], n.patches)					### power law exponent (repeat the same value for each patch)
    pars$prob.disp<-prob.disp						### Annual probability of dispersal
    pars$thres<-rep(input[s,5], n.patches)					### Threshold distance for network
    pars$min.dist<-apply(D, 1, function(x){min(x, na.rm=T)})		### Distance to the nearest neighbour patch (in km)
    pars$patch.conn<-degree(g)						### Patch connectivty measured as the number of neighbours
    pars$tot.conn<-sum(pars$patch.conn)					### Total landscape connectivity
    pars$prop.hab<-rep(sum(A)/(max(dim.sq)^2),  n.patches)			### Proportion of the landscape that is 'good' habitat
    pars$mag<-cat.mag							### Magnitude of the catastrophe
    outsim[[s]]<-pars
  }
  outsim
}




############ SIMULATIONS
### Number of simulations and sensitivity analysis parameters based on Prowse et al. (2016): http://onlinelibrary.wiley.com/doi/10.1002/ecs2.1238/full 

n.sims <- 10				#### Number of simulations: 20000; PGD: it is currently 10 so it is easy to see if the code works
dimensions <- seq(0, 10, 0.01)		#### Dimensions of the landscape (in kilometres because this is the unit used by the dispersal kernel),total area =10000ha
years <- 13*3				#### Number of years simulated, representing three generations; maximum lifespan of a possum in the wild = 13 years 
cat.years <- 10

##### RUN THE SIMULATIONS

# sim.patch <- metapop.model(n.sim=n.sims, dim.sq=dimensions,  n.year=years)

sim.patch <- metapop.nallee(n.sims, dimensions, years, cat.years, T)
# sim.patch  ### Delete the hastag to examine the list containing the results :)

sim.patch <- metapop.model(n.sims, dimensions, years, T, cat.years, thres.dens = .2)

### Store the simulation data in a more manageable format (data frame); PGD: the output simulation have to have the same length for this to work

out.sim <- do.call(rbind, lapply(sim.patch, data.frame, stringsAsFactors=FALSE))

##################################################### PGD: I've changed this to a sensitivity analysis of the abundance in the last simulated year
############################################# Not quite sure it makes any sense, and it probably better to analyse something else. For example,
############################################# if you run the sa, you'll find that k is the most important factor. Of course it is!
#################################### I've just copied the new thing (too lazy, again), but we may need to think about it. 

################################ SENSITIVITY ANALYSIS OF THE ABUNDANCE IN THE LAST YEAR
library(dismo)

### Create a data frame storing the data

data.sens<-data.frame(N=out.sim$N, Area=out.sim$Area, dist.nearest=out.sim$min.dist, connectivity=out.sim$conn, n0=out.sim$init, k=out.sim$k, 
			alpha=out.sim$alpha, c=out.sim$c, z=out.sim$z, prob.disp=out.sim$prob.disp, r=out.sim$r, prop.hab=out.sim$prop.hab)

### Extract only complete cases - excluding those that produced NAs (because these simulations were exploring a restricted parameter space)

data.sens<-data.sens[complete.cases(data.sens),]


### Fit a boosted regression tree (takes 1 hour in my laptop)

sens<-gbm.step(data.sens, gbm.x=2:12, gbm.y=1, tree.complexity=3,learning.rate=0.001, bag.fraction=0.7, step.size=1000, max.trees=4000,  n.folds=3, family="poisson")

### Examine the results and variable importance
summary(sens)



##### let's test the logistic equation

logistic_growth <- function(N, r, k){

  print(N)
  print(round((r*N)*((1-N)/k), digits=0))
  N1 <- N + round((r*N)*((1-N)/k), digits=0)
  print(N1)
  N1
}

N <- 6461
r <- .1 #-0.7587540  #runif(1, min=-2, max=2)
k <- 6450 #round(N + runif(1, min=2, max=10),digits=0)

dynamics <- c()
for(i in 1:20){

  dynamics <- append(dynamics, N)
  N <- logistic_growth(N,r,k)
}

plot(dynamics, type='l')



##################### MODEL 11/04/2017 ################
#this function incorporates allee effects to the model
metapop.model <- function(n.sim, dim.sq, n.year, with.allee=F, cat.year, thres.dens){
  ### Load the required libraries: fields for the pairwise euclidean distance and lhs for the latin hypercube sampling
  require(fields)
  require(igraph)
  require(qgraph)
  require(msm)
  
  #### structure for holding the results
  outsim <- vector("list") 
  for (s in 1:n.sim){			### Start simulations
    pars <- vector("list") 														### Create empty list to store results
	  ### Actual parameters values
	  n.patches <- round(rtnorm(1, mean=10, sd=3, lower=5, upper=50), digits=0)		### Number of patches drawn from a truncated normal distribution (rounding a normal dist approximates a Poisson dist)
	  density<- rpois(1, lambda=500)  								### Density per km2 (alpha); mean density = 500/km2. From estimates in a report by Warburton et al. (2009r
	  c <- runif(1, min=1/10, max=1/5)					        ### Exponential distance decay dispersal kernel in metres (from 1 to 5 km); based on estimates from Etherington et al. (2014): http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0088293
	  z <- runif(1, min=0.1, max=1)							### Power law exponent (for the abundance-area relationship)
    A <- round(rtnorm(n.patches, mean=0.02, sd=1.5, lower=0.01), digits=2)		### Patch area sampled from a truncated normal distribution (mean = 2ha, sd=1.5, and minumum size of 1 ha)
    cat.mag<-round(rbeta(n.patches, 1, 1), digits=2)							### Magnitude of the catastrophe (proportion remaing after catastrophe) - flat probability from 0 to 1
    
    ### Create patch centroids and calculate the distance matrix
    x.sim <- sample(dim.sq, n.patches)						### Simulate x-coordinates for the centroids using sample without replacement	
    y.sim <- sample(dim.sq, n.patches)						### Simulate y-coordinates for the centroids using sample without replacement	
    coord.patches <- matrix(c(x.sim, y.sim), nrow=n.patches, ncol=2)		### Put together the patch centroids
    D <- rdist(coord.patches)							### Pairwise euclidean distance matrix calculated using function rdist from package fields
    
    ##### Initial population size and k (carrying capacity) for each patch
    N <- matrix(0, ncol=n.patches, nrow=n.year)					### Create empty matrix to store results
    N[1,] <- k <- rpois(n.patches, density*A^z)				### Define initial population size based on the power-law relationship; initialise the local population at k
    
    ### for the allee effect
    if(with.allee) allee <- 20					### This is an example value that seems to work (we can try other things :) 
    
    #### Exponential distance-decay dispersal kernel 
    dk <- exp(-c*D)
  	diag(dk) <- 0
	
    #### I was talking with Dean and a good way of normalising things so they sum up to 1 is to use a Multinomial distribution
	  #### That's what I've done below within the 'year' loop
	
    ### Probability of dipsersal in each patch
    prob.disp <- round(rbeta(n.patches, 2, 5), digits=2)					### Annual probability of dispersal using a Beta distribution with a peak density around ~0.2
    
    ### Create an array to store the results of the movement patterns (a matrix each year)
    movement <- array(0, dim=c(nrow=n.patches, ncol=n.patches, n.year))
    r <- runif(n.patches, min=0.1, max=3)					### Per capita population growth rate (due to births and death, but not because of migration - migration is modelled directly in raw numbers)
    
    #### Simulate the dynamics
    #### Here start the simulations for each + 1 year
    for (i in 2:n.year){        ### n.year: number of years simulated
  	  #### PGD-: If the dispersal is after growth, I think that we need to calculate the growth first and the dispersal based on that. I've changed the allee effects formulat to the one in the post
  	  #### which is the same as model (4) in Table 1 of the paper (it shows interesting dynamics - figures 1b and 1c in paper).
  	  #### Also, I've removed some brackets - it may be because I'm crazy and seeing many brackets makes me nervous
  	 
      if(with.allee) growth <- rpois(n.patches, N[i-1, ] + round((r * N[(i-1),])* (1-(N[(i-1),]/k))*((N[(i-1),]/k)-(allee/k))))
	    else growth <- rpois(n.patches, N[i-1, ] +  round(((r * N[(i-1),]  ) * ((k - N[(i-1),] )/k) )))
      
	    growth[is.na(growth)] <- 0
      ### Movement after growth. Here we calculate dispersers for each local population and store the results in the
      ### movement matrix
	    ### Emigrants per year drawn from a Binomial distribution
      
		  emigrants <- rbinom(n.patches, size=growth, prob=prob.disp)
		
	    ### Estimate the movement matrix (total number of current dispersers) drawn from a multinomial distribution (so the probabilities sum up to 1 row-wise and there is no need to normalise)
	    ### I'm kind of certain that there has to be an easier way of doing this without using a loop, but I couldn't figure it out. Sleep and caffeine deprivation?
	  
	    for (j in 1:n.patches){
	      movement[j,,i] <- rmultinom(1, size=emigrants[j], prob=dk[j,])
      }
	  
      emigrants <- rowSums(movement[,,i], na.rm=T)
      immigrants <- colSums(movement[,,i], na.rm=T)
      
      print('current pop')
      print(N[i-1,])
      print('emigrants')
      print(emigrants)
      print('immigrants')
      print(immigrants)
      
      #### Now the population size with growth and dispersal. ML: If dispersal has to be after growth, I would do it like this. PGD: You're the expert! I'm happy wit that
      if(with.allee) N[i,] <-  rpois(n.patches, growth - emigrants + immigrants)
      else N[i,] <- rpois(n.patches, growth - emigrants + immigrants)
      #N[i,][N[i,]< 0] <- 0
      
      if(i == cat.year){
        N[i, ] <- rbinom(n.patches, N[i-1, ], prob=cat.mag)
	    }
	
    } 
    # Define the network of connections between patches observed from the random process is to set a threshold 
    # this is equivalent to say that no connection is assumed 
    # between patches for which the probability to arrive from 
    # one to the other is smaller than the threshold
	  #### PG-D: Changed this again to km (see email)
	  conn_thrs <- 5
    adj_m <- D
    adj_m[adj_m < conn_thrs] <- 0
    adj_m[is.na(adj_m)] <- 0
    adj_m[adj_m > conn_thrs] <- 1
    
    ### this is how the network of patches looks like
    g <- graph_from_adjacency_matrix(adj_m, mode='undirected')
    
    #### so we can se results a little bit...
    matplot(N[,1:5], type='l', lwd=2)
    
    ### time until recovering to a threshold density of interest
    t.recov <- rep(0, n.patches)
	  for (j in 1:n.patches){
      t.recov[j] <- which(N[cat.year:n.year,j]/A[j] >= thres.dens*100)[1]					
	  }

    ### OUTPUT VALUES- classified by groups
    pars$t.recov <- t.recov									### Time to recovery
    ### Landscape-level parameters
    pars$n.patches <- n.patches									### Number of patches
    pars$smwi <- as.numeric(smallworldness(g)[1])						### Network small world index
	  pars$tot.conn <- ecount(g)/((vcount(g)*(vcount(g)-1))/2)			### Connectance
	  pars$lp <- ecount(g)/vcount(g)									### Mean number of links per patch in the landscape
	  pars$prop.hab <- rep((sum(A))/(max(dim.sq)^2),  n.patches)			### Proportion of the landscape that is 'good' habitat
	
	  ### Patch-level parameters (Miguel - add centrality here, please :)
	  pars$Area<-A														### Patch area, in km
	  pars$min.dist<-apply(D, 1, function(x){sort(x, decreasing=TRUE)[2]})		### Distance to the nearest neighbour patch (in km), discarding the first one that's going to be itself
	  pars$patch.conn<-degree(g)										### Patch connectivty measured as the number of neighbours
	  pars$k<-k														### Patch carrying capacity
	  
	  #### ML: As centrality measure I think betweenness is more appropriate for our 
	  #### research question, since it tells you how many paths cross through the node, 
	  #### or in other words, thorugh which patch it is more likely that possums will 
	  #### have to travel to be able to colonise the whole metapop
	  pars$between <- betweenness(g)      
	  #### I also like this version because it normalises by the degree, so less connected patches
	  #### are more important if they have a high centrality
	  pars$norm_between <- betweenness(g)/degree(g)
	  
	  ### Species-level parameters
	  pars$c<-rep(c, n.patches)								### Exponential distance decay dispersal kernel in metres (repeat the same value for each patch)
    pars$prob.disp<-prob.disp										### Annual probability of dispersal
    pars$alpha<-rep(density, n.patches)							### Density, alpha (repeat the same value for each patch)
    pars$z<-rep(z, n.patches)								### power law exponent (repeat the same value for each patch)

	  ### Catastrophe magnitude (patch-level)
	  pars$mag<-1-cat.mag												### Magnitude of the catastrophe
	  pars$N<-N
    pars$g<-g
	
  	outsim[[s]]<-pars
  }
  outsim
}

