compute.gE <- function(phi){
  
  gE <- matrix(NA,12,1);
  tmp <- daylength.factor.from.lat(phi,TRUE);
  L <- tmp$L;
  ndl <- tmp$ndl
  cdays <- tmp$cdays
  #
  for (t in 1:12){
    gE[t] <-  mean(ndl[(cdays[t]+1):cdays[t+1]]);
  }
  return(gE)
}

daylength.factor.from.lat <- function(phi,return.ndl.and.cdays=FALSE){
  latr <- phi*pi/180;  # change to radians
  ndays <- cbind(0,31,28,31,30,31,30,31,31,30,31,30,31);
  cdays <- cumsum(ndays);
  sd <- t(asin(sin(pi*23.5/180) * sin(pi * (((1:365) - 80)/180))));   # solar declination
  y <- -tan(matrix(latr,365,1)) * t(tan(sd));
  # bound y within (-1,1):
  y[y >= 1] <- 1;
  y[y <= -1] <- -1;
  
  hdl <- acos(y);
  dtsi <- hdl * sin(matrix(latr,365,1)) * t(sin(sd)) +
    (matrix(cos(latr),365,1)) * t(cos(sd)) * sin(hdl);
  ndl <- dtsi/max(dtsi); # normalized day length
  
  # calculate mean monthly daylength (used for evapotranspiration in soil moisture calcs)
  jday <- cdays[1:12] +.5*ndays[2:13];
  m.star <- 1-tan(phi*pi/180)*tan(23.439*pi/180*cos(jday*pi/182.625));
  # bound m.star between 0 and 2:
  m.star[m.star < 0] <- 0
  m.star[m.star > 2] <- 2;
  
  nhrs <- 24*acos(1-m.star)/pi; # the number of hours in the day in the middle of the month
  # mean normalized daylength factor:
  L <- (ndays[2:13]/30) * (nhrs/12)
  if(return.ndl.and.cdays){
    out <- list(L,ndl,cdays)
    names(out) <- c("L","ndl","cdays")
    return(out)
  }else{
    return(L); 
  }
}

leakybucket.monthly <- function(syear,eyear,phi,T,P,Mmax = 0.76,Mmin = 0.01,alph = 0.093,
                                m.th = 4.886,mu.th = 5.8,rootd = 1000,M0 = .2){
  
  iyear <- syear:eyear;
  nyrs <- length(iyear);
  
  # Storage for growth response output variables (size [12 x Nyears]):
  M <- potEv <- matrix(NA,12,nyrs);
  
  if(M0 < 0.){M0 <- 200/rootd;}
  
  # Compute normalized daylength (neglecting small difference in calculation for leap-years)
  L <- daylength.factor.from.lat(phi); 
  
  # Pre-calculation of istar and I, using input T to compute the climatology:
  Tm <- rowMeans(T);
  istar <- (Tm/5)^1.514; 
  istar[Tm < 0] <- 0;
  I <- sum(istar);
  
  # precalculation of the exponent alpha in the Thornwaite (1948) equation:
  a <- (6.75e-7)*I^3 - (7.71e-5)*I^2 + (1.79e-2)*I + 0.49;
  
  #########################################################################################
  #### -- year cycle -- ####
  # syear = start (first) year of simulation
  # eyear = end (last) year of simulation
  # cyear = year the model is currently working on
  # iyear = index of simulation year
  
  for (cyear in 1:nyrs){     # begin cycling over years
    #########################################################################################
    for (t in 1:12){  # begin cycling over months in a year
      
      ##### Compute potential evapotranspiration for current month after Thornthwaite:
      if ( T[t,cyear] < 0 ){Ep = 0;}
      if ( T[t,cyear] >= 0 && T[t,cyear] < 26.5 ){Ep <- 16*L[t]*(10*T[t,cyear]/I)^a;}
      if ( T[t,cyear] >= 26.5 ){Ep <- -415.85 + 32.25*T[t,cyear] - .43* T[t,cyear]^2;}
      potEv[t,cyear] <- Ep;
      
      ##### Now calculate soil moisture according to the CPC Leaky Bucket model
      ##### (see J. Huang et al, 1996).
      
      if (t > 1){
        # evapotranspiration:
        Etrans <- Ep*M[t-1,cyear]*rootd/(Mmax*rootd);
        # groundwater loss via percolation:
        G <- mu.th*alph/(1+mu.th)*M[t-1,cyear]*rootd;
        # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
        R <- P[t,cyear]*(M[t-1,cyear]*rootd/(Mmax*rootd))^m.th +
          (alph/(1+mu.th))*M[t-1,cyear]*rootd;
        dWdt <- P[t,cyear] - Etrans - R - G;
        M[t,cyear] <- M[t-1,cyear] + dWdt/rootd;
      }
      
      if( t == 1 && cyear > 1){
        # evapotranspiration:
        Etrans <- Ep*M[12,cyear-1]*rootd/(Mmax*rootd);
        # groundwater loss via percolation:
        G <- mu.th*alph/(1+mu.th)*M[12,cyear-1]*rootd;
        # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
        R <- P[t,cyear]*(M[12,cyear-1]*rootd/(Mmax*rootd))^m.th +
          (alph/(1+mu.th))*M[12,cyear-1]*rootd;
        dWdt <- P[t,cyear] - Etrans - R - G;
        M[t,cyear] <- M[12,cyear-1] + dWdt/rootd;
      }
      
      if (t == 1 && cyear == 1){
        if (M0 < 0){ M0 <- .20;}
        # evapotranspiration (take initial soil moisture value to be 200 mm)
        Etrans <- Ep*M0*rootd/(Mmax*rootd);
        # groundwater loss via percolation:
        G <- mu.th*alph/(1+mu.th)*(M0*rootd);
        # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
        R <- P[t,cyear]*(M0*rootd/(Mmax*rootd))^m.th + (alph/(1+mu.th))*M0*rootd;
        dWdt <- P[t,cyear] - Etrans - R - G;
        M[t,cyear] <- M0 + dWdt/rootd;
      }
      
      # error-catching:
      if (M[t,cyear] <= Mmin) {M[t,cyear] <- Mmin;}
      if (M[t,cyear] >= Mmax) {M[t,cyear] <- Mmax;}
      if (is.na(M[t,cyear])==1){ M[t,cyear] <- Mmin;}
    } # end month (t) cycle
    #########################################################################################
  } # end year cycle
  
  return(list(M=M,potEv=potEv))
  
}

std.ramp <- function(x,x1,x2){return(
  apply(as.matrix(apply((x-x1)/(x2-x1),1:length(dim(x)),min,1)),
        1:length(dim(x)),max,0)
)}

VSLite <- function(syear,eyear,phi,T,P,
                   T1 = 8, T2 = 23, M1 = .01, M2 = .05,
                   Mmax = 0.76,Mmin = 0.01,alph = 0.093,
                   m.th = 4.886,mu.th = 5.8,rootd = 1000,M0 = .2,
                   substep = 0,I_0 = 1,I_f = 12,hydroclim = "P"){
  #############################################################################
  nyrs <- length(syear:eyear)
  Gr <- gT <- gM <- M <- potEv <- matrix(NA,12,nyrs);
  #############################################################################
  
  ## Load in soil moisture, or estimate it with the Leaky Bucket model:
  if(hydroclim == "M"){
    ## Read in soil moisture:
    M = P;
  }else{# Compute soil moisture:
    if(substep == 1){
      M <- leakybucket.submonthly(syear,eyear,phi,T,P,
                                  Mmax,Mmin,alph,m.th,mu.th,rootd,M0);
    }else{
      M <- leakybucket.monthly(syear,eyear,phi,T,P,
                               Mmax,Mmin,alph,m.th,mu.th,rootd,M0);
    }
    if(substep !=1 && substep != 0){
      cat("'substep' param must either be set to 1 or 0.");
      return
    }
  }
  
  # Compute gE, the scaled monthly proxy for insolation:
  gE <- compute.gE(phi);
  
  #############################################################################
  ### Calculate Growth Response functions gT and gM
  
  # Temperature growth response:
  gT <- std.ramp(T,T1,T2)
  
  # Soil moisture growth response:
  gM <- std.ramp(M$M,M1,M2)
  
  # Compute overall growth rate:
  Gr <- kronecker(matrix(1,1,nyrs),gE)*pmin(gT,gM)
  
  ############## Compute proxy quantity from growth responses #################
  width <- matrix(NA,nyrs,1);
  if (phi>0){ # Site in Northern Hemisphere:
    if (I_0<0){ # if we include part of the previous year in each year's modeled growth:
      startmo <- 13+I_0;
      endmo <- I_f;
      # use average of growth data across modeled years to estimate first year's growth due
      # to previous year:
      width[1] <- sum(Gr[1:endmo,1]) + sum(rowMeans(Gr[startmo:12,]));
      for(cyear in 2:nyrs){
        width[cyear] <- sum(Gr[startmo:12,cyear-1]) + sum(Gr[1:endmo,cyear]);
      }
    }else{ # no inclusion of last year's growth conditions in estimates of this year's growth:
      startmo <- I_0+1;
      endmo <- I_f;
      width <- colSums(Gr[startmo:endmo,])
    }
  }
  if(phi<0){ # if site is in the Southern Hemisphere:
    # (Note: in the Southern Hemisphere, ring widths are dated to the year in which growth began!)
    startmo <- 7+I_0; # (eg. I_0 = -4 in SH corresponds to starting integration in March of cyear)
    endmo <- I_f-6; # (eg. I_f = 12 in SH corresponds to ending integraion in June of next year)
    for (cyear in 1:(nyrs-1)){
      width(cyear) <- sum(Gr[startmo:12,cyear]) + sum(Gr[1:endmo,cyear+1]);
    }
    # use average of growth data across modeled years to estimate last year's growth due
    # to the next year:
    width[nyrs] <- sum(Gr[startmo:12,nyrs])+sum(rowMeans(Gr[1:endmo,]));
  }
  
  # Simulated proxy series standardized width:
  trw <- t((width-mean(width))/sd(width)); 
  
  #############################################################################
  # Return output:
  out <- list(trw = trw, gT = gT, gM = gM, gE = gE, Gr=Gr, M = M$M, potEv = M$potEv,
              sample.mean.width = mean(width), sample.std.width = sd(width))
  return(out)
  
}

