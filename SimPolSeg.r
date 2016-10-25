
# This is the function that runs the MAP theory simulation
# The imputs are "areas" = number of geographical areas
# "blackpop" = minority population of all geographic areas combined
# "sd_blackpop" = Standard deviation of the minority population in all areas
# "whitepop" = White population  of all geographic areas combined
# "sd_whitepop" = Standard de

segsim<-function(areas, blackpop,sd_blackpop, whitepop, sd_whitepop,Id, sims){
  e = 0.001
  library(msm)
  library(micEcon)
  library(sp)
  library(ggplot2)
  
  # Tipping point function
  mstar<-function(eta) {
    eta/2  # Eta is tolerance
  }
  
  # Loss function
  loss<-function(eta,ma1,d) {
    eta*ma1 - ma1^2  - d 
  }
  
  # Output
  countydat<-data.frame()
  move_out_white<-c()
  simulations<-c()
  move_out_black<-c()
  moveout_county<-data.frame()
  
  
  # Generate the initial time period agents
  race.s0<-c()
  countynum.s0<-c()
  Ia.s0<-c()
  rho<-c()
  blackpopup.s0<-c()
  prep<-c()
  countypop.s0<-c()
  blackpop.s0<-c()
  whitepop.s0<-c()
  time_period.s0<-c()
  pctblk.s0<-c()
  blackpopgrow.s0<-c()
  ones.s0<-c()   
  xc.s0<-c() # x-coordinate
  yc.s0<-c() # y-coordinate
  
  # Ideology is between zero and one , set it equal to .5
  Ideology<-rtnorm(
    n = areas, 
    mean = Id, 
    sd= .1, 
    lower = 0 + e, 
    upper= 1 - e)
  
  
  
  if(length(whitepop) == 1) {
    # Generate blackpop for each area
    whitepop<-round(
      rtnorm(n = areas,mean = whitepop,sd = sd_whitepop, lower = 1),
      digits=0)
    
    blackpop<-round(
      rtnorm(n = areas,mean = blackpop, sd = sd_blackpop),
      
      digits=0)
  }
  
  whitepop<-round(whitepop,digits = 0)
  blackpop<-round(blackpop+1, digits = 0)
  
  
  
  countypop<-round(whitepop + blackpop,digits=0)
  
  
  pernum.s0<-1:sum(countypop)
  countySnum<-1:areas
  tipping.s0 <-rep(0,sum(countypop))
  eta_norm_scale.s0 <-rep(0,sum(countypop))
  whitepopgrow.s0 <-rep(0,sum(countypop))
  #Small number to get values above the threshold
  
  # Create spatial coordinates
  # Create the spatial coordinates for the areas
  xc<-round(runif(areas),4)
  yc<-round(runif(areas),4)
  
  
  # Generate the initial dataset , white agents
  for(i in 1:areas) {
    
    time_period.s0<-c(time_period.s0,rep(0,whitepop[i]))
    countynum.s0<-c(countynum.s0,rep( i,whitepop[i]))
    Ia.s0<-c(Ia.s0,rtnorm(n=whitepop[i],mean = Ideology[i],sd=.1,lower = 0 + e, upper = 1-e))
    blackpop.s0 <-c(blackpop.s0,rep(blackpop[i],whitepop[i]))
    whitepop.s0 <-c(whitepop.s0,rep(whitepop[i],whitepop[i]))
    countypop.s0 <-c(countypop.s0,rep(countypop[i],whitepop[i]))
    pctblk.s0<-c(pctblk.s0,rep(blackpop[i]/countypop[i],whitepop[i]))
    race.s0<-c(race.s0, rep(1,whitepop[i]))
    blackpopgrow.s0<- c(blackpopgrow.s0,rep(0,whitepop[i]))
    ones.s0<-c(ones.s0, rep(1,whitepop[i]))
    # Add the area coordinates 
    xc.s0<-c(xc.s0,rep(xc[i],whitepop[i]))
    yc.s0<-c(yc.s0,rep(yc[i],whitepop[i]))
  }
  
  # Generate the initial dataset , black agents
  for(i in 1:areas) {
    time_period.s0<-c(time_period.s0,rep(0,blackpop[i]))
    countynum.s0<-c(countynum.s0,rep( i,blackpop[i]))
    Ia.s0<-c(Ia.s0,rtnorm(n=blackpop[i],mean = Ideology[i] + .3,sd=.1,lower = 0 + e, upper = 1-e))  # Blacks do not move
    blackpop.s0 <-c(blackpop.s0,rep(blackpop[i],blackpop[i]))
    whitepop.s0 <-c(whitepop.s0,rep(whitepop[i],blackpop[i]))
    countypop.s0 <-c(countypop.s0,rep(countypop[i],blackpop[i]))
    pctblk.s0<-c(pctblk.s0,rep(blackpop[i]/countypop[i],blackpop[i]))
    race.s0<-c(race.s0, rep(0,blackpop[i]))
    blackpopgrow.s0<- c(blackpopgrow.s0,rep(0,blackpop[i]))
    ones.s0<-c(ones.s0, rep(1,blackpop[i]))
    # Add the area coordinates 
    xc.s0<-c(xc.s0,rep(xc[i],blackpop[i]))
    yc.s0<-c(yc.s0,rep(yc[i],blackpop[i]))
  }
  
  
  # Generate initial agent data frame t=0
  agentdat<-cbind(pernum.s0,
                  time_period.s0,
                  countynum.s0,
                  Ia.s0,
                  blackpop.s0, #This is blackpop + migration shock
                  whitepop.s0,
                  countypop.s0,
                  pctblk.s0,
                  tipping.s0,
                  eta_norm_scale.s0,
                  whitepopgrow.s0 ,
                  race.s0,
                  blackpopgrow.s0,
                  ones.s0,
                  xc.s0,
                  yc.s0)
  
  agentdat<-agentdat[order(agentdat[,3]),]
  
  finagent<-agentdat                                     
  # Generate initial county data frame t= 0
  countydat<-data.frame(pernum.s0,
                        time_period.s0,
                        countynum.s0,
                        Ia.s0,
                        blackpop.s0, #This is blackpop + migration shock
                        whitepop.s0,
                        countypop.s0,
                        pctblk.s0,
                        tipping.s0,
                        eta_norm_scale.s0,
                        whitepopgrow.s0,
                        race.s0,
                        blackpopgrow.s0,
                        ones.s0,
                        xc.s0,
                        yc.s0)
  
  
  countydat<-aggregate(countydat,list(countydat$countynum.s0),FUN=mean)
  
  
  
  # The Spatial Dataset
  xy.sp<-SpatialPoints(cbind(xc,yc))
  xy.countydat<-SpatialPointsDataFrame(xy.sp, countydat) 
  
  
  
  
  
  
  #countydat<-cbind(countydat,xc,yc)
  countydat<-as.matrix(countydat[,2:dim(countydat)[2]]) # This makes rbinding easier
  
  
  #par(mfrow=c(1,2))
  print("Running simulation for time")
  
  for(j in 1:sims) {
    
    print(paste("...",j))
    
    agentdat<-agentdat[order(agentdat[,3]),]
    
    pernum.s1<-agentdat[,1]
    time_period_sim.s1<-agentdat[,2]
    countynum_sim.s1<-agentdat[,3]
    Ia_sim.s1<-agentdat[,4]
    blackpop_sim.s1<-round(agentdat[,5] ,digits=0)
    whitepop_sim.s1<-round(agentdat[,6],digits=0)
    countypop_sim.s1<-round(agentdat[,7],digits=0)
    pctblk_sim.s1<-agentdat[,8]
    whitepopgrow.s1<- agentdat[,11]
    blackpopgrow.s1<- agentdat[,13]
    race.s1<-agentdat[,12]
    ones.s1<-agentdat[,14]
    
    ###########################################################
    # MOving Decisions of people in agentdat                   #
    ###########################################################
    
    
    D<-.01
    eta_norm_scale_white<-Ia_sim.s1 - D
    
    eta_norm_scale_white<-ifelse(eta_norm_scale_white <=0, 0 +e, eta_norm_scale_white)
    
    eta_norm_scale_white<- eta_norm_scale_white * race.s1
    
    eta_norm_scale_black<-(Ia_sim.s1) * abs(1-race.s1)
    eta_norm_scale<-eta_norm_scale_white + eta_norm_scale_black
    
    # Generate Tipping Points
    tipping<-mstar(eta_norm_scale)  # Calculate tipping points
    tipping<-tipping*race.s1 + abs(1-race.s1)
    
    # Moving: Rule 1  Move if tipping point is above the threshold
    move<-ifelse(pctblk_sim.s1 > tipping,1,0)
    
    # How many people moved?
    # Number of people for all areas
    
    whitemovers<-sum(move[race.s1 == 1])
    blackmovers<-sum(move[race.s1 == 0])
    totalpop<-length(move)
    pctmovers<-whitemovers/totalpop
    
    print(paste("whites:",whitemovers, "blacks:", blackmovers, "pctwhitemover:",pctmovers))
    
    moverow<-cbind(whitemovers, blackmovers,pctmovers)
    moveout_county<-rbind(moveout_county,moverow)
    
    if(pctmovers <0.01) {
      break
    }
    
    
    
    simulations<-c(simulations,j)
    
    
    
    #  Agent level data for time t-1
    
    agentdat[,9]<-tipping
    agentdat[,10]<-eta_norm_scale
    
    
    tc<-aggregate(agentdat,list(agentdat[,3]),FUN = mean)
    
    tc<-as.matrix(tc[,2:dim(tc)[2]]) # First row is unnecessary
    
    tc<-tc[order(tc[,3]),]
    # This applies to the ORIGINAL Agentdat without the newly moved in people
    counties<-tc[,3] #County numbers
    pctblkcounty<-tc[,8]  # Percent Black for each county
    id_county<-tc[,4]
    xco<- tc[,15]
    yco<-tc[,16]
    tc2<-cbind(counties,pctblkcounty,id_county,xco,yco)
    
    
    # Give a list of eligible counties to relocate to
    newcount<-lapply(tipping, function (x) tc2[ pctblkcounty < x, ]  )
    
    
    for(i in 1:length(newcount)) {
      if(length(newcount[[i]]) == 0) {
        newcount[[i]] <- tc2[tc2[,2] == min(pctblkcounty),]
      }
      
    }
    
    # Add the spatial coordinates and eta to calculate distance between current location and prospective locations
    # Get the county that has the highest utility 
    
    
    fincount<-c() # Final counties for movers
    
    for(i in 1:length(newcount)) {
      
      if(length(newcount[[i]]) == 5) { # This accounts for when tipping points happen to not be satisfied
        fincount<-c(fincount,newcount[[i]][1])
      }
      else if(length(newcount[[i]]) > 5) {
        
        
        xnow<-as.vector(
          rep(agentdat[i,15],dim(newcount[[i]])[1])
        )
        ynow<-as.vector(
          rep(agentdat[i,16],dim(newcount[[i]])[1])
        )
        xother<-newcount[[i]][,4]
        yother<-newcount[[i]][,5]
        d<-sqrt(
          (xnow - xother)^2 + (ynow-yother)^2
        )
        
        # Eta
        eta.agent<-rep(agentdat[i,10],dim(newcount[[i]])[1])
        
        # Minority pop
        mpct<-newcount[[i]][,2]
        
        # Utility of moving to a certain place
        move.utility<- eta.agent*mpct - mpct^2 - d
        
        movematrix<-data.frame(newcount[[i]],d,eta.agent,move.utility)
        
        fincount<-c(fincount, movematrix$counties[movematrix$move.utility == max(movematrix$move.utility)])
        
      }
    }
    
    
    
    
    
    # Generates the counties that are moved to
    newcount<- fincount *move  + abs(1-move)*countynum_sim.s1 # Isolates only those counties where people decided to move
    newcount<-as.vector(newcount)
    
    ######################################################################
    # County level data for time t  , population growth                  #
    ######################################################################
    # New County Data for t=0
    
    # New counties for people currently in agentdat
    countynum_sim.s2<-newcount  # After Moving counties
    
    # This updates the new dataset with the correct demographics
    s2agent<-agentdat
    
    for(i in 1:length(countynum_sim.s2)) {
      if(countynum_sim.s2[i] != countynum_sim.s1[i]){
        poprow<-tc[tc[,3] == countynum_sim.s2[i],5:8]
        s2agent[i,5:8]<-poprow
        # Also need to change the coordinates for those that moved
        coordrow<-tc[tc[,3] == countynum_sim.s2[i],15:16]
        s2agent[i,15:16]<-coordrow
      } 
    }
    
    
    
    
    # Update the countynums
    s2agent[,3]<-countynum_sim.s2
    s2agent<-s2agent[order(s2agent[,3]),]
    blacks<-abs(1-race.s1)
    
    whitepop_n<-rep(0,areas)
    blackpop_n<-rep(0,areas)
    
    for(i in 1:areas) {
      whitepop_n[i]<-sum(race.s1[countynum_sim.s2 == i])
      blackpop_n[i]<-sum(blacks[countynum_sim.s2 == i])
    }
    
    whitepop_sim.s2<-c()
    blackpop_sim.s2<-c()
    
    
    countypop_n<-whitepop_n + blackpop_n
    
    for(i in 1:length(countypop_n)){
      whitepop_sim.s2<-c(whitepop_sim.s2,rep(round(whitepop_n[i],digits=0),countypop_n[i]))
      blackpop_sim.s2<-c(blackpop_sim.s2,rep(round(blackpop_n[i],digits=0),countypop_n[i]))
    }
    
    
    # Update the whitepop
    s2agent[,6]<-whitepop_sim.s2
    s2agent[,5]<-blackpop_sim.s2
    
    # New Variables within Each County
    time_period_sim.s2<- s2agent[,2] + 1
    Ia_sim.s2<- s2agent[,4]
    blackpop_sim.s2<-s2agent[,5]
    blackpop_sim.s2<-ifelse((is.na(blackpop_sim.s2) == TRUE | blackpop_sim.s2 <1), 1, blackpop_sim.s2)
    blackpop_sim.s2<-round(blackpop_sim.s2,digits = 0)
    whitepop_sim.s2<-s2agent[,6]
    whitepop_sim.s2<-ifelse((is.na(whitepop_sim.s2) == TRUE | whitepop_sim.s2 < 1), 1, whitepop_sim.s2)
    whitepop_sim.s2<-round(whitepop_sim.s2,digits = 0)
    countypop_sim.s2<-whitepop_sim.s2 + blackpop_sim.s2
    pctblk_sim.s2<-blackpop_sim.s2/countypop_sim.s2
    whitepopgrow.s2<- 0.01* s2agent[,6]
    whitepopgrow.s2<-ifelse((is.na(whitepopgrow.s2) == TRUE | whitepopgrow.s2 < 1), 1, whitepopgrow.s2)
    whitepopgrow.s2<-round(whitepopgrow.s2,digits = 0 )
    blackpopgrow.s2<-1/j^2* s2agent[,5]
    blackpopgrow.s2<-ifelse((is.na(blackpopgrow.s2) == TRUE | blackpopgrow.s2 <1), 1, blackpopgrow.s2)
    blackpopgrow.s2<-round(blackpopgrow.s2,digits = 0)
    race.s2<-s2agent[,12]
    ones.s2<-s2agent[,14]
    xc.s2<-s2agent[,15]
    yc.s2<-s2agent[,16]
    
    # Final Agent Dataset
    finagent<-rbind(finagent,as.matrix(s2agent))
    
    # New Agent County Data AFTER Moving
    move_agent1<-data.frame(pernum.s1,
                            time_period_sim.s2,
                            countynum_sim.s2,
                            Ia_sim.s2,
                            blackpop_sim.s2,
                            whitepop_sim.s2,
                            countypop_sim.s2,
                            pctblk_sim.s2,
                            tipping,
                            eta_norm_scale,
                            whitepopgrow.s2,
                            race.s2 ,
                            blackpopgrow.s2,
                            ones.s2,
                            xc.s2,
                            yc.s2)
    
    move_agent1<-move_agent1[order(move_agent1$countynum_sim.s2),]
    
    move_county<-aggregate(move_agent1,list(move_agent1$countynum_sim.s2), FUN=mean)
    move_county<-data.frame(move_county[,2:dim(move_county)[2]]) # First row is unnecessary
    
    move_county<-move_county[order(move_county$countynum_sim.s2),]
    
    sd_county<-as.vector(by(move_agent1[,4], move_agent1[,3], sd))
    
    # Final County Dataset
    countydat<-rbind(countydat,as.matrix(move_county))
    #countydatplot<-data.frame(countydat)
    
    #p<-ggplot(countydatplot,aes(countydatplot[,2],countydatplot[,4]))
    #p + geom_point(colour=countydatplot[,4])
    
    plot(countydat[,2], countydat[,4],xlab = "Time", ylab = "Avg. Area Ideology") 
    
    
    
    # New Agent County Data AFTER Moving
    blackpop_sim.s2<-blackpop_sim.s2 + blackpopgrow.s2
    whitepop_sim.s2<-whitepop_sim.s2 + whitepopgrow.s2
    countypop_sim.s2<-whitepop_sim.s2 + blackpop_sim.s2
    
    move_agent2<-data.frame(pernum.s1,
                            time_period_sim.s2,
                            countynum_sim.s2,
                            Ia_sim.s2,
                            blackpop_sim.s2,
                            whitepop_sim.s2,
                            countypop_sim.s2,
                            pctblk_sim.s2,
                            tipping,
                            eta_norm_scale,
                            whitepopgrow.s2,
                            race.s2 ,
                            blackpopgrow.s2,
                            ones.s2,
                            xc.s2,
                            yc.s2)
    
    move_county2<-aggregate(move_agent2,list(move_agent2$countynum_sim.s2), FUN=mean)
    move_county2<-data.frame(move_county2[,2:dim(move_county2)[2]]) # First row is unnecessary
    
    move_county2<-move_county2[order(move_county2$countynum_sim.s2),]
    
    sd_county2<-as.vector(by(move_agent2[,4], move_agent2[,3], sd))
    
    ######################################################################
    # Agent level data for time t                                        #
    ######################################################################
    # First: Generate New Black and White Agents with updated preferences
    
    pernum.s3<-c()
    time_period_sim.s3<-c()
    countynum_sim.s3<-c()
    Ia_sim.s3<-c()
    blackpop_sim.s3<-c()
    whitepop_sim.s3<-c()
    countypop_sim.s3<-c()
    pctblk_sim.s3<-c()
    tipping.s3<-c()
    eta_norm_scale.s3<-c()
    whitepopgrow.s3<-c()
    num.s3<-j
    blackpopgrow.s3<-c()
    race.s3<-c()
    ones.s3<-c()
    xc.s3<-c()
    yc.s3<-c()
    
    # This loop generates data for the new white agents
    
    for(i in 1:areas){
      
      if( round(move_county2$whitepopgrow.s2[i],digits=0) > 0) {
        time_period_sim.s3<-c(time_period_sim.s3,rep(num.s3,round(move_county2$whitepopgrow.s2[i],digits=0)))
        countynum_sim.s3<-c(countynum_sim.s3,rep(move_county2$countynum_sim.s2[i],round(move_county2$whitepopgrow.s2[i],digits=0)))
        Ia_sim.s3<-c(Ia_sim.s3,rtnorm(n=round(move_county2$whitepopgrow.s2[i],digits=0),mean = move_county2$Ia_sim.s2[i],sd = .1 ,lower = 0+e, upper = 1-e))
        blackpop_sim.s3 <-c(blackpop_sim.s3,rep(round(move_county2$blackpop_sim.s2[i],digits=0),round(move_county2$whitepopgrow.s2[i],digits=0)))
        whitepop_sim.s3 <-c(whitepop_sim.s3,rep(round(move_county2$whitepop_sim.s2[i],digits=0),round(move_county2$whitepopgrow.s2[i],digits=0)))
        countypop_sim.s3 <-c(countypop_sim.s3,rep(round(move_county2$countypop_sim.s2[i],digits=0),round(move_county2$whitepopgrow.s2[i],digits=0)))
        pctblk_sim.s3<-c(pctblk_sim.s3,rep(move_county2$pctblk_sim.s2[i],round(move_county2$whitepopgrow.s2[i],digits=0)))
        tipping.s3<-c(tipping.s3,rep(0,round(move_county2$whitepopgrow.s2[i],digits=0)))
        eta_norm_scale.s3<-c(eta_norm_scale.s3,rep(0,round(move_county2$whitepopgrow.s2[i],digits=0)))
        whitepopgrow.s3<-c(whitepopgrow.s3,rep(0,round(move_county2$whitepopgrow.s2[i],digits=0)))
        pernum.s3<-c(pernum.s3,rep(0,round(move_county2$whitepopgrow.s2[i],digits=0)))
        blackpopgrow.s3<-c(blackpopgrow.s3,rep(0,round(move_county2$whitepopgrow.s2[i],digits=0)))
        race.s3<-c(race.s3,rep(1,round(move_county2$whitepopgrow.s2[i],digits=0)))
        ones.s3<-c(ones.s3,rep(1,round(move_county2$whitepopgrow.s2[i],digits=0)))
        xc.s3<-c(xc.s3,
                 rep(move_county2$xc.s2[i], round(move_county2$whitepopgrow.s2[i],digits=0))
        )
        yc.s3<-c(yc.s3,
                 rep(move_county2$yc.s2[i], round(move_county2$whitepopgrow.s2[i],digits=0))
        )
      }
      else {
        next
      }
      
      
    }
    
    for(i in 1:areas){
      if( round(move_county2$blackpopgrow.s2[i],digits=0) > 0) {
        time_period_sim.s3<-c(time_period_sim.s3,rep(num.s3,round(move_county2$blackpopgrow.s2[i],digits=0)))
        countynum_sim.s3<-c(countynum_sim.s3,rep(move_county2$countynum_sim.s2[i],round(move_county2$blackpopgrow.s2[i],digits=0)))
        Ia_sim.s3<-c(Ia_sim.s3,rtnorm(n=round(move_county2$blackpopgrow.s2[i],digits=0),mean = move_county2$Ia_sim.s2[i] + .3 ,sd = .1,lower = 0+ e, upper = 1 - e))
        blackpop_sim.s3 <-c(blackpop_sim.s3,rep(round(move_county2$blackpop_sim.s2[i],digits=0),round(move_county2$blackpopgrow.s2[i],digits=0)))
        whitepop_sim.s3 <-c(whitepop_sim.s3,rep(round(move_county2$whitepop_sim.s2[i],digits=0),round(move_county2$blackpopgrow.s2[i],digits=0)))
        countypop_sim.s3 <-c(countypop_sim.s3,rep(round(move_county2$countypop_sim.s2[i],digits=0),round(move_county2$blackpopgrow.s2[i],digits=0)))
        pctblk_sim.s3<-c(pctblk_sim.s3,rep(move_county2$pctblk_sim.s2[i],round(move_county2$blackpopgrow.s2[i],digits=0)))
        tipping.s3<-c(tipping.s3,rep(0,round(move_county2$blackpopgrow.s2[i],digits=0)))
        eta_norm_scale.s3<-c(eta_norm_scale.s3,rep(0,round(move_county2$blackpopgrow.s2[i],digits=0)))
        whitepopgrow.s3<-c(whitepopgrow.s3,rep(0,round(move_county2$blackpopgrow.s2[i],digits=0)))
        pernum.s3<-c(pernum.s3,rep(0,round(move_county2$blackpopgrow.s2[i],digits=0)))
        blackpopgrow.s3<-c(blackpopgrow.s3,rep(0,round(move_county2$blackpopgrow.s2[i],digits=0)))
        race.s3<-c(race.s3,rep(0,round(move_county2$blackpopgrow.s2[i],digits=0)))
        ones.s3<-c(ones.s3,rep(1,round(move_county2$blackpopgrow.s2[i],digits=0)))
        xc.s3<-c(xc.s3,
                 rep(move_county2$xc.s2[i], round(move_county2$blackpopgrow.s2[i],digits=0))
        )
        yc.s3<-c(yc.s3,
                 rep(move_county2$yc.s2[i], round(move_county2$blackpopgrow.s2[i],digits=0))
        )
      }
      else {
        next
      }
      
    }
    
    
    newpop<-cbind(pernum.s3,
                  time_period_sim.s3,
                  countynum_sim.s3,
                  Ia_sim.s3,
                  blackpop_sim.s3,
                  whitepop_sim.s3,
                  countypop_sim.s3,
                  pctblk_sim.s3,
                  tipping.s3,
                  eta_norm_scale.s3,
                  whitepopgrow.s3,
                  race.s3,
                  blackpopgrow.s3,
                  ones.s3,
                  xc.s3,
                  yc.s3)
    
    # Attaches the new agents                                             
    agentdat_t<-rbind(as.matrix(move_agent2),newpop)
    agentdat_t<-agentdat_t[order(agentdat_t[,3]),]
    
    
    
    
    
    # To use for the next loop
    agentdat<-agentdat_t
    
    
    
    
    
  }
  
  countydat<-data.frame(countydat[,-c(1,11,13,14)])
  finagent<-data.frame(finagent[,-c(1,11,13,14)])
  names(countydat)<- c("time", "area"  , "Ia", "blackpop", "whitepop","countypop", "pctblk", "tipping", "eta_norm_scale","race")
  names(finagent)<-c("time", "area" ,  "Ia", "blackpop", "whitepop","countypop", "pctblk", "tipping", "eta_norm_scale","race")
  
  return(list(finagent=finagent, countydat =countydat, moveout_county = moveout_county))
}

#######################################################################
#######################################END#############################
#######################################################################











