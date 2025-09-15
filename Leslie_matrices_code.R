#########
# Stable-population model
# for the forecasting of obesity prevalence
# across generations 
# WITH 4 BMI CLASSES (< 18.5, [18.5-25), [25-30), >= 30)
# (Aldea, Garcia-Aguirre, Beltran-Sanchez, Daza, & Palloni)
####################

#clean etc.
rm(list = ls())
gc()
memory.limit(size=10^10)

#how many generations and simulations
# per combination of parameters
gen <- 10
sim <- 20

#initial generation (derived from HRS couples, individuals
# with PRS and BMI, restricted data)
data <- readRDS("~/ECHO/Obesity/Submission HH/Revision/2nd/example_of_initial_generation.rds")

#add BMI category
data$cat <- ifelse(data$bmi<18.5, 1,
                 ifelse(data$bmi<25, 2,
                        ifelse(data$bmi < 30,3,4)))

#number of omegas (AM parameter)
oms <- 2
#number of phis (DF paremeter)
phs <- 2
#number of ps (increased genetic penetrance)
ps <- 3

# array for stocking BMI distributions
dis <- array(0:0, c(gen, 4, sim, oms, phs, ps))

for (i in 1:oms){
  for (j in 1:phs){
    for (k in 1:ps){
      for (sm in 1:sim){
        dis [1,,sm, i, j, k] <- table(data$cat)/nrow(data)
      }
    }
  }
}


#heritability matrix
her <- array(0:0, c(16,4))

#loop through parameter values

for (pp in 1:ps){
  #set p
  p <- (pp-1)*0.1
  
  for (w in 1:oms){
    #set omega
    omega <- seq(0,1, by = 1/(oms-1))[w]
    
    for (ph in 1:phs){
      #set phi
      phi <- seq(0,1, by = 1/(phs-1))[ph]
      
      #loop through simulations
      for (i in 1:sim){
        
        dt <- data
        
        #loop through generations
        for (t in 2:gen){
          
          #1--------------
          # Mating matrix
          bmis <- dis [t-1, , i, w, ph, pp]
          
          mat <- c(bmis[1]*(bmis[1]+(1-bmis[1])*omega), 
                   bmis[1]*bmis[2]*(1-omega), 
                   bmis[1]*bmis[3]*(1-omega), 
                   bmis[1]*bmis[4]*(1-omega),
                   bmis[2]*bmis[1]*(1-omega), 
                   bmis[2]*(bmis[2]+(1-bmis[2])*omega), 
                   bmis[2]*bmis[3]*(1-omega), 
                   bmis[2]*bmis[4]*(1-omega),
                   bmis[3]*bmis[1]*(1-omega), 
                   bmis[3]*bmis[2]*(1-omega), 
                   bmis[3]*(bmis[3]+(1-bmis[3])*omega), 
                   bmis[3]*bmis[4]*(1-omega),
                   bmis[4]*bmis[1]*(1-omega), 
                   bmis[4]*bmis[2]*(1-omega), 
                   bmis[4]*bmis[3]*(1-omega), 
                   bmis[4]*(bmis[4]+(1-bmis[4])*omega))
          
          #2--------------
          # Fertility matrix (values for diagonal)
          f <- array(0:0, c(16))
          for (a in 1:4){
            for (b in 1:4){
              f[(a-1)*4+b] <- 2 + (a+b-6)*phi/2
            }
          }
          
          #3. Children
          #--------------
          child <- mat*f
          
          #4--------------
          # Microsimulation to feed H matrix
          #--------------
          
          #4A. Making couples
          
          #establish target couples of same BMI category
          tar <- array(0:0, c(4)) 
          for (c in 1:4){
            tar [c] <- floor(length(which(dt$cat==c)) * omega / 2)
          }
          
          #couples dataset
          coup <- data.frame(matrix(nrow = 0, ncol = 6))
          rest <- data.frame(matrix(nrow = 0, ncol = 5))
          
          #make same BMI category couples
          if (omega > 0){
            for (c in 1:4){
              df <- subset(dt, cat == c)
              who <- sample(1:nrow(df), tar[c]*2, replace = F)
              for (ss in 1:tar [c]){
                coup <- rbind(coup, c(unlist(df[who[ss], c("bmi","pgs", "cat")]),
                                      unlist(df[who[tar[c] + ss], c("bmi","pgs", "cat")])))
              }
              rest <- rbind(rest, df[-who,])
            }
          }else { rest <- dt }
          
          #make random couples
          howm <- floor(nrow(rest)/2)
          
          if (howm>0){
            who <- sample(1:nrow(rest), howm*2, replace = F)
            for (ss in 1:howm){
              coup <- rbind(coup, c(unlist(rest[who [ss], c("bmi","pgs", "cat")]),
                                    unlist(rest[who [howm + ss], c("bmi","pgs", "cat")])))
            }
          }
          
          colnames(coup) <- c("bmi1", "pgs1", "cat1", "bmi2", "pgs2", "cat2")
          coup$type <- paste(coup$cat1, coup$cat2, sep="")
          
          # 4B. Giving children to couples
          
          coup$ch <- rpois(nrow(coup), 2 + (coup$cat1+coup$cat2-6)*phi/2)
          
          # 4C. Give children PGS & CRS
          
          new_gen <- data.frame(matrix(nrow = 0, ncol = 3))
          
          for (c in 1:nrow(coup)){
            if (coup$ch [c] > 0){
              new_gen <- rbind(new_gen, cbind(cbind(rep(coup$type[c],coup$ch[c])), 
                                              cbind(0.5*coup$pgs1[c]+0.5*coup$pgs2[c] + rnorm(n = coup$ch[c], 0, sd = sqrt(0.5)))))
            }
          }
          
          colnames(new_gen) <- c("type", "pgs")
          new_gen$pgs <- as.numeric(new_gen$pgs)
          
          #compute BMI and category
          
          new_gen$bmi <- 27.796 + 1.0715 * new_gen$pgs * (1+p) + 
            rnorm (n = nrow(new_gen), 0, sd = 4.55) 
          
          new_gen$cat <- ifelse(new_gen$bmi<18.5, 1,
                                ifelse(new_gen$bmi<25, 2,
                                       ifelse(new_gen$bmi < 30,3,4)))
          
          # 5. Compute heritability matrix 
          # (with parental and child BMI categories)
          #--------------------------------
          
          
          for (z in 1:4){
            for (y in 1:4){
              ww <- subset(new_gen, type == paste(z,y,sep=""))
              if (nrow(ww) != 0){
                her [(z-1)*4 + y, 1] <- length(which(ww$cat == 1))/nrow(ww)
                her [(z-1)*4 + y, 2] <- length(which(ww$cat == 2))/nrow(ww)
                her [(z-1)*4 + y, 3] <- length(which(ww$cat == 3))/nrow(ww)
                her [(z-1)*4 + y, 4] <- length(which(ww$cat == 4))/nrow(ww)
              }else{
                # if there are no individuals in this case and it's the first time
                # step, children inherit the BMI class of their parents
                if (t == 2){
                  if (z == y){
                    for (x in 1:4){
                      if (x == z){
                        her [(z-1)*4 + y, x] <- 1
                      }else{
                        her [(z-1)*4 + y, x] <- 0
                      }
                    }
                  }else{
                    for (x in 1:4){
                      if (x %in% c(y,z)){
                        her [(z-1)*4 + y, x] <- 0.5
                      }else{
                        her [(z-1)*4 + y, x] <- 0
                      }
                    }
                  }
                }
              }
              
            }
          }
          
          #stock result
          obs <- array(0:0, c(4))
          for (j in 1:4){
            obs [j] <- sum (child * her [,j])  
          }
          
          dis[t, , i, w, ph, pp] <- obs / sum(obs)
          
          #prepare for next time step
          dt <- new_gen [, c("pgs", "bmi", "cat")]
          
          #In order to reduce computation time, one can set a population limit 
          # and randomly remove individuals over it.
          if (nrow(dt) > 10000){
            dt <- dt [sample(1:nrow(dt), 10000, replace = F),]
          }
          
        }
      }
    }
  }
}

