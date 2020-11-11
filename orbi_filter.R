

is.between <- function(x,a,b) {
  
  if(!is.na(a) & !is.na(b)) return(x>=a & x<=b)
  if(is.na(a) & is.na(b))   return(rep(TRUE,length(x))) 
  if(is.na(a)) return(x<=b)
  if(is.na(b)) return(x>=a)
  
}






orbifilter <- function(x,windows_width=0.1,max_rel_int = 0.2,keep_isotopes=TRUE,max_charge=5,isotope_mz_tol = 0.005){

  neutron   <- 1.0033548378 # really C12, C13 difference
  iso_dist  <- neutron / seq(from=1,by=1,to=max_charge)
  
  to_rem    <- rep(FALSE,  nrow(x)  )
  done      <- rep(FALSE,  nrow(x)  )
  int_order <- order(x[,"intensity"], decreasing = TRUE)
  
  
  # To save time in the loop. Cuts about 20%
  find_isotopes <- keep_isotopes & any(iso_dist<(windows_width/2)) 
  windows_width_half <- windows_width/2
  mz <- x[,"mz"]
  int <- x[,"intensity"]
  
  
  while( !all(done | to_rem) ){
    target <- int_order[   !(int_order %in% which(done | to_rem))   ][1] # highest intensity mass that has not already been remove or assessed
    
    rem_candidate <- which(is.between(mz     ,mz[target]-windows_width_half    ,mz[target]+windows_width_half    )) # the index of mz values within windows_width  
    rem_candidate <- rem_candidate[   int[rem_candidate] / int[target] < max_rel_int  ] # only remove peaks that have lower relative intensity than max_rel_int
    
    
    if(  find_isotopes   ){ # Keep masses if they could be isotopes. Also multicharged. Don't do the calculations if the combination of max_charge and windows_width makes isotopes impossible
    target_dist      <- abs(mz[rem_candidate]-mz[target])
    dist_matrix      <- outer(target_dist, iso_dist, "-")
    dist_matrix      <- abs(dist_matrix)
    dist_to_iso_hypo <- apply(dist_matrix,1,min)
    rem_candidate    <- rem_candidate[dist_to_iso_hypo>isotope_mz_tol]
    }
    
    
    to_rem[rem_candidate] <- TRUE
    done[target]          <- TRUE
  }
  
  
  new_scan <- x[!to_rem,,drop=FALSE]
  return(new_scan)  
  
}






xcmsRaw.orbifilter <- function(xraw,windows_width=0.1,max_rel_int = 0.2,keep_isotopes=TRUE,max_charge=5,isotope_mz_tol = 0.005){
  
   
  # Get all scans seperately
  scans <- lapply(1:length(xraw@scanindex),function(i) getScan(xraw,i))
  
  # remove masses from scans
  scans_new <- lapply(scans,function(x) orbifilter(x,windows_width,max_rel_int,keep_isotopes,max_charge,isotope_mz_tol))
  
  # Create new scanindex
  scanindex_new=numeric(length=length(scans_new))
  
  for (i in 1:length(scans_new)){
    if (i==1){
      scanindex_new[i]=0
    }else{
      scanindex_new[i] = scanindex_new[i-1] + nrow(scans_new[[(i-1)]])
    }
  }
  
  
  
  # Put data into new object
  ob<-new("xcmsRaw")
  ob@env = new.env(parent=.GlobalEnv)
  
  scans_new=do.call(rbind,scans_new)
  ob@env$mz<-scans_new[,'mz']
  ob@env$intensity<-scans_new[,'intensity']
  
  ob@scantime <-xraw@scantime
  ob@scanindex <-as.integer(scanindex_new)
  ob@polarity   <- xraw@polarity
  ob@acquisitionNum <- 1:length(ob@scanindex)
  ob@profmethod     <- xraw@profmethod
  ob@mzrange        <- range(ob@env$mz)
  ob@filepath       <- xraw@filepath
  profStep(ob) =   profStep(xraw)
  ob=xcms:::remakeTIC(ob)
  
  
  # and we're done
  return(ob)
  
  
}
