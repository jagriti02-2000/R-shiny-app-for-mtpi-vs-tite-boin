#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  observeEvent(input$simulate, {
    
    pava <- function(x, wt = rep(1, length(x))) {
      n <- length(x)
      if (n <= 1) {
        return(x)
      }
      if (any(is.na(x)) || any(is.na(wt))) {
        stop("Missing values in 'x' or 'wt' not allowed")
      }
      lvlsets <- (1:n)
      repeat {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) {
          break
        }
        i <- min((1:(n - 1))[viol])
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i + 1]
        ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
        x[ilvl] <- sum(x[ilvl] * wt[ilvl]) / sum(wt[ilvl])
        lvlsets[ilvl] <- lvl1
      }
      x
    }
    # betavar computes variances of beta distributions
    betavar <- function(a, b) {
      resp <- a * b / ((a + b)^2 * (a + b + 1))
      return(resp)
    }
    
    
    tite_boin <- function(true_dlt_prob = c(), no_dose, No_of_simulations, taget_tox, Cohort_size, dlt_window, max_sample, max_dose_sample, enroll_time) {
      set.seed(0205)
      startdose <- 1
      d <- startdose
      simN <- No_of_simulations
      pT <- taget_tox ## Target toxicity probability
      true_dlt_prob <- true_dlt_prob
      
      D <- no_dose
      DLTwindow <- dlt_window
      phi <- taget_tox
      phi_1 <- 0.6 * phi
      phi_2 <- 1.4 * phi
      lambda_e <- log((1 - phi_1) / (1 - phi)) / log((phi * (1 - phi_1)) / (phi_1 * (1 - phi)))
      lambda_d <- log((1 - phi) / (1 - phi_2)) / log((phi_2 * (1 - phi)) / (phi * (1 - phi_2)))
      a <- 0.5 * phi # alpha in paper
      b <- 1 - a # beta in paper
      Totaltime <- rep(0, simN)
      ##  No. of doses
      Dsize <- max_dose_sample
      N.max <- max_sample
      
      datan <- matrix(rep(0, simN * D), ncol = D)
      datax <- matrix(rep(0, simN * D), ncol = D)
      res <- rep(0, simN)
      alldosesafe <- 0
      alldosetoxic <- 0
      
      # csize=3
      for (sim in 1:simN) {
        all_data <- data.frame(matrix(ncol = 6, nrow = 0))
        colnames(all_data) <- c("Cohort", "Patient", "Enrolment time", "Dose", "DLT Happens", "Total assessment Time")
        st <- 0
        nodose <- 0
        seldose <- 0
        complete_data <- data.frame(matrix(ncol = 6, nrow = 0))
        colnames(complete_data) <- c("Cohort", "Patient", "Enrolment time", "Dose", "DLT Happens", "Total assessment Time")
        pending_data <- data.frame(matrix(ncol = 6, nrow = 0))
        colnames(pending_data) <- c("Cohort", "Patient", "Enrolment time", "Dose", "DLT Happens", "Total assessment Time")
        particular_dose <- data.frame(matrix(ncol = 6, nrow = 0))
        colnames(particular_dose) <- c("Cohort", "Patient", "Enrolment time", "Dose", "DLT Happens", "Total assessment Time")
        
        totaltime <- 0
        cohort <- 1
        num <- rep(0, D)
        d <- startdose
        st <- 0
        nodose <- 0
        toxdose <- D + 1
        seldose <- 0
        totaltime <- 0
        x <- rep(0, D)
        
        while (st == 0) {
          newvar <- c()
          
          csize <- Cohort_size
          position <- rep(0, csize)
          position2 <- c()
          # st = 1 indicates the trial must be terminated
          # nodose = 1 indicates no dose is selected as MTD
          # toxdose indicates unacceptably toxic dose
          enrol_time <- runif(n = csize, min = 0, max = enroll_time)
          # enrol_time=rexp(n=3,rate=1/3)
          can <- sort(round(enrol_time, 1)) + totaltime
          for (h in 1:csize) {
            if (anyDuplicated(can) > 0) {
              can[anyDuplicated(can)] <- can[anyDuplicated(can)] + 0.1
            }
          }
          data <- data.frame(cohort, Patient = c(1:csize), can, dose = d)
          names(data)[3] <- "Enrolment time"
          xx <- rbinom(1, csize, true_dlt_prob[d])
          
          if (xx > 0) {
            timetoDLT <- round(runif(n = xx, min = 0.1, max = DLTwindow), 1) # Calculation time to persons DLT
            personwithDLT <- sample(can, xx) # Finding Which person had DLT
            for (k in 1:xx) {
              position[k] <- which(can == personwithDLT[k]) # Finding the positon of the person who had DLT
              newvar[position[k]] <- timetoDLT[k] + personwithDLT[k] # #Calculating the total assessment time for the person who had DLT
            }
            
            newcan <- can[-c(position)] # Finding the position of the person who did not had DLT
            for (j in 1:length(newcan)) {
              position2[j] <- which(can == newcan[j])
              if (xx < csize & xx > 0) {
                newvar[position2[j]] <- can[position2[j]] + DLTwindow
              }
            } # Calculating the total assessment time who did not have DLT
            position <- sort(position)
            data["DLT_Happens"] <- 0
            for (i in 1:csize) {
              for (j in 1:length(position)) {
                if (data$Patient[i] == position[j]) {
                  data$DLT_Happens[i] <- 1 # Adding a indicator if the following person had DLT
                }
              }
            }
            data$totalassessmenttime <- newvar # Adding the total time column
          }
          
          # if there is no DLT
          
          if (xx == 0) { # Calculating the total assessment time when there is no DLT
            data["DLT_Happens"] <- 0
            data$totalassessmenttime <- can + DLTwindow
          }
          
          # Sorting the dataset w.r.t time
          
          datasorted <- data[order(data$totalassessmenttime), ]
          # data=datasorted[ceiling(length(datasorted$totalassessmenttime)/2),]
          all_data <- rbind(all_data, datasorted)
          number <- cohort
          # getting all the information of the current doe
          particular_dose <- all_data %>% filter(cohort <= number, dose == d)
          particular_dose_copy <- particular_dose
          # sorting he data set
          particular_dose <- particular_dose[order(particular_dose$totalassessmenttime), ]
          TFT <- particular_dose$totalassessmenttime[ceiling(nrow(particular_dose) / 2)] - particular_dose$`Enrolment time`
          TFT <- TFT[ceiling(nrow(particular_dose) / 2) + 1:nrow(particular_dose)]
          TFT <- TFT[!is.na(TFT)]
          
          extra_time <- 0
          DLT_half_count <- particular_dose$DLT_Happens[1:ceiling(nrow(particular_dose) / 2)]
          if (min(TFT) < (DLTwindow)/2) {
            extra_time <- ((DLTwindow)/2 - min(TFT))
          }
          # TFT=TFT+extra_time
          # STFT=sum(TFT)/DLTwindow
          particular_dose$totalassessmenttime <- particular_dose$totalassessmenttime + extra_time
          # dataused$totalassessmenttime=dataused$totalassessmenttime+3.8
          # write.csv(dataused,"~/Codes")
          counter <- 0
          for (j in 1:length(TFT)) {
            if (particular_dose$totalassessmenttime[ceiling(nrow(particular_dose) / 2)] >= particular_dose_copy$totalassessmenttime[ceiling(nrow(particular_dose_copy) / 2) + j]) {
              counter <- counter + 1
            }
          }
          BOIN <- 0
          if (particular_dose$totalassessmenttime[ceiling(nrow(particular_dose) / 2)] + counter == nrow(particular_dose)) {
            BOIN <- 1
          }
          complete_data <- head(particular_dose, ceiling(nrow(particular_dose) / 2) + counter)
          # Taking the information from the pending patients
          pending_data <- tail(particular_dose, nrow(particular_dose) - ceiling(nrow(particular_dose) / 2) - counter)
          # TFT=dataused1$totalassessmenttime[nrow(dataused1)]-dataused2$`Enrolment time`
          TFT <- particular_dose$totalassessmenttime[ceiling(nrow(particular_dose) / 2)] - pending_data$`Enrolment time`
          # TFT[TFT<0]=0
          STFT <- sum(TFT) / DLTwindow
          
          
          # estimating Phat
          c <- sum(TFT != 0)
          n <- nrow(complete_data) + c # number of subjects
          r <- nrow(complete_data) # subjects with Completed data
          
          # c = n - r # subjects who have not completed DLT assessment
          s <- sum(complete_data$DLT_Happens) # DLT count
          p_mean <- (s + a) / (r + a + b) # posterior mean
          ######### BOIN############################
          # p_hat=0
          if (BOIN == 1 || c == 0) {
            p_hat <- s / n
          } else {
            p_hat <- (s + (p_mean / (1 - p_mean)) * (c - STFT)) / n
          }
          p_hat[p_hat > 1] <- 1
          num[d] <- n # storing no. of patients at a particular dose
          x[d] <- s
          
          ############ Comparing###############################
          if (d == 1) {
            ## if the first dose is too toxic, the trial will be terminated and no dose is selected as MTD
            if ((1 - pbeta(q = phi, shape1 = 1 + x[d], shape2 = num[d] - x[d] + 1)) > 0.95 && num[d] >= 3) {
              st <- 1
              nodose <- 1
              alldosetoxic <- alldosetoxic + 1
            } else {
              if (p_hat >= lambda_d) {
                st <- 1
                nodose <- 1
                alldosetoxic <- alldosetoxic + 1
              }
              if ((p_hat > lambda_e) && (p_hat < lambda_d) && (num[d] < Dsize)) {
                d <- d
              }
              if (p_hat <= lambda_e) {
                d <- d + 1
              }
            }
          } else {
            if (d == D) {
              if ((1 - pbeta(q = phi, shape1 = 1 + x[d], shape2 = num[d] - x[d] + 1)) > 0.95 && num[d] >= 3) {
                toxdose <- d
                d <- d - 1
              } else {
                if (p_hat >= lambda_d) {
                  d <- d - 1
                }
                if ((p_hat > lambda_e) && (p_hat < lambda_d) && (num[d] < Dsize)) {
                  d <- d
                }
                if ((p_hat <= lambda_e) && (num[d] >= Dsize)) {
                  st <- 1
                  nodose <- 1
                  alldosesafe <- alldosesafe + 1
                } #### all dose levels are cleared and nodose is selected as MTD
                if ((p_hat <= lambda_e) && (num[d] < Dsize)) {
                  d <- d
                }
              }
            } else {
              if ((d > 1) && (d < D)) {
                if ((1 - pbeta(q = phi, shape1 = 1 + x[d], shape2 = num[d] - x[d] + 1)) > 0.95 && num[d] >= 3) {
                  toxdose <- d
                  d <- d - 1
                } else {
                  if (p_hat >= lambda_d) {
                    d <- d - 1
                  }
                  if ((p_hat > lambda_e) && (p_hat < lambda_d) && (num[d] < Dsize)) {
                    d <- d
                  }
                  if (p_hat <= lambda_e) {
                    d <- d + 1
                  }
                }
              }
            }
          }
          
          total <- sum(num)
          totaltime <- particular_dose$totalassessmenttime[ceiling(nrow(particular_dose) / 2)]
          if ((total >= N.max)) {
            st <- 1
            nodose <- 0
          } #### stop the trial when maximum sample size is reached
          if ((d >= 1) && (d <= D) && (p_hat < lambda_d && p_hat > lambda_e) && (num[d] >= Dsize)) {
            st <- 1
          }
          if (num[d] >= Dsize) {
            st <- 1
          }
          cohort <- cohort + 1
        }
        Totaltime[sim] <- totaltime
        if (nodose == 0) {
          tdose <- toxdose - 1
          pp <- rep(-100, tdose)
          pp.var <- rep(0, tdose)
          for (i in 1:tdose) {
            pp[i] <- (x[i] + .005) / (num[i] + .01)
            pp.var[i] <- betavar(x[i] + 0.005, num[i] - x[i] + 0.005) ### here adding 0.005 is to use beta(0.005, 0.005) for estimating the MTD, which is different from the dose-finding.
          }
          pp <- pava(pp, wt = 1 / pp.var) # perform the isotonic transformation using PAVA with weights being posterior variances
          for (i in 2:tdose) {
            pp[i] <- pp[i] + i * 1E-10 # by adding an increasingly small number to toxicity probability at higher doses, it will break the ties and make the lower dose level the MTD if the ties were larger than pT or make the higher dose level the MTD if the ties are smaller than pT
          }
          seldose <- c(1:tdose)[sort(abs(pp[1:tdose] - pT), index.return = T)$ix[1]]
          # seldose is the final MTD that is selected based on the order-transformed posterior means
        }
        res[sim] <- seldose
        for (i in 1:D) {
          datan[sim, i] <- num[i]
          datax[sim, i] <- x[i]
        }
      }
      
      abc <- rep(0, D)
      for (i in 1:D) {
        abc[i] <- sum(res == i) / simN # abc is the proportion of each dose is selected as the MTD
      }
      
      
      tite_boin_data <- data.frame(matrix(ncol = 2, nrow = D+5))
      colnames(tite_boin_data) <- c("True Toxicity Probability", "TITE BOIN")
      
      tite_boin_data$`True Toxicity Probability` <-c(true_dlt_prob,"Trials Selecting All Doses as Safe (%)",
                                                     "Trials Selecting All Doses as Toxic (%)",
                                                     "Average time taken (In weeks)",
                                                     
                                                     "Average Sample size","Average Overdosed patients (%)")
      
      patients_at_each_dose <- round(apply(datan,2,mean),3)
      taget_position <- which(min(abs(true_dlt_prob - pT))==(true_dlt_prob - pT))
      over_dose_per  <- (sum(patients_at_each_dose[taget_position+1:D],na.rm = TRUE)*100)/round(sum(apply(datan, 2, mean)), 3)
      
      new_row <- c( round((abc * 100), 2),round(((alldosesafe / simN) * 100), 2), round(((alldosetoxic / simN) * 100), 2), round(sum(Totaltime / simN), 2), round(sum(apply(datan, 2, mean)), 3),round(over_dose_per,2))
      tite_boin_data$`TITE BOIN` <- new_row
      
      return(tite_boin_data)
    }
    
    mtpi_simulation_8weeks <- function(true_dlt_prob = c(), no_dose, No_of_simulations, taget_tox, Cohort_size, dlt_window, max_sample, max_dose_sample, enroll_time) {
      set.seed(0205)
      true_dlt_prob <- true_dlt_prob
      # simulation setup
      simN <- No_of_simulations ## No. of simulations
      pT <- taget_tox ## Target toxicity probability
      D <- length(true_dlt_prob) ##  No. of doses
      Dsize <- 9
      N.max <- 50 ## Maximum sample size
      startdose <- 1 ## Starting dose
      d <- startdose
      eps1 <- 0.05
      eps2 <- 0.05 ## (pT-eps1, pT+eps2) is the equivalence interval
      a <- 1
      b <- 1 ##  Prior is Beta(1,1).
      xi <- 0.95 ##  Cutoff probability for excessive toxicity
      datan <- matrix(rep(0, simN * D), ncol = D)
      datax <- matrix(rep(0, simN * D), ncol = D)
      res <- rep(0, simN)
      alldosesafe <- 0
      alldosetoxic <- 0
      
      DLTwindow <- 8
      st <- 0
      time <- 1:8
      Totaltime <- c()
      # Simulations
      for (sim in 1:simN) {
        x <- rep(0, D)
        n <- rep(0, D)
        pa <- rep(0, D)
        pb <- rep(0, D)
        q <- rep(0, 3)
        d <- startdose
        st <- 0
        nodose <- 0
        toxdose <- D + 1
        seldose <- 0
        totaltime <- 0
        # st = 1 indicates the trial must be terminated
        # nodose = 1 indicates no dose is selected as MTD
        # toxdose indicates unacceptably toxic dose
        while (st == 0) {
          csize <- Cohort_size
          newvar <- c()
          position <- c()
          xx <- rbinom(1, csize, true_dlt_prob[d])
          
          enrol_time <- runif(n = csize, min = 0, max = enroll_time)
          can <- sort(round(enrol_time, 1))
          for (h in 1:csize) {
            if (anyDuplicated(can) > 0) {
              can[anyDuplicated(can)] <- can[anyDuplicated(can)] + 0.1
            }
          }
          
          if (xx > 0) {
            timetoDLT <- round(runif(n = xx, min = 0.1, max = 8), 1)
            personwithDLT <- sample(can, xx)
            for (k in 1:xx) {
              position[k] <- which(can == personwithDLT[k])
            }
          }
          if (xx == 0) {
            can <- can + DLTwindow
            
            assessmenttime <- can
          }
          if (xx > 0) {
            for (j in 1:xx) {
              newvar[j] <- personwithDLT[j] + timetoDLT[j]
            }
            newcan <- can[-c(position)]
            if (xx < 3 & xx > 0) {
              newcan <- newcan + 8
            }
            assessmenttime <- c(newcan, newvar)
          }
          x[d] <- x[d] + xx
          n[d] <- n[d] + csize
          # Update prior beta distribution to obtain posterior distribution
          pa[d] <- x[d] + a
          pb[d] <- n[d] - x[d] + b
          pa[d + 1] <- x[d + 1] + a
          pb[d + 1] <- n[d + 1] - x[d + 1] + b
          # Computation to see if the next dose is too toxic
          if (d < D) {
            temp <- 1 - pbeta(pT, a + x[d + 1], b + n[d + 1])
            if (temp > xi) {
              tt <- 1
              toxdose <- d + 1
            } else {
              tt <- 0
            }
          }
          # Compute the UPM for three intervals
          q[1] <- (1 - pbeta(eps2 + pT, pa[d], pb[d])) / (1 - eps2 - pT) ### overdosing interval
          q[2] <- (pbeta(eps2 + pT, pa[d], pb[d]) - pbeta(pT - eps1, pa[d], pb[d])) / (eps2 + eps1) ### Proper dosing interval
          q[3] <- (pbeta(pT - eps1, pa[d], pb[d]) / (pT - eps1)) * (1 - tt) ### underdosing interval
          # implement the dose-assignment rules based on the UPM
          if (d == 1) {
            ## if the first dose is too toxic, the trial will be terminated and no dose is selected as MTD
            if ((1 - pbeta(pT, a + x[d], b + n[d] - x[d])) > xi) {
              st <- 1
              nodose <- 1
              alldosetoxic <- alldosetoxic + 1
            } else {
              if ((q[1] > q[2]) && (q[1] > q[3])) {
                st <- 1
                nodose <- 1
                alldosetoxic <- alldosetoxic + 1
              }
              if ((q[2] > q[1]) && (q[2] > q[3]) && (n[d] < Dsize)) {
                d <- d
              }
              if ((q[3] > q[1]) && (q[3] > q[2])) {
                d <- d + 1
              }
            }
          } else {
            if (d == D) {
              if ((q[1] > q[2]) && (q[1] > q[3])) {
                d <- d - 1
              }
              if ((q[2] > q[1]) && (q[2] > q[3]) && (n[d] < Dsize)) {
                d <- d
              }
              if ((q[3] > q[1]) && (q[3] > q[2]) && (n[d] >= Dsize)) {
                st <- 1
                nodose <- 1
                alldosesafe <- alldosesafe + 1
              } #### all dose levels are cleared and nodose is selected as MTD
              if ((q[3] > q[1]) && (q[3] > q[2]) && (n[d] < Dsize)) {
                d <- d
              }
            } else {
              if ((d > 1) && (d < D)) {
                if ((q[1] > q[2]) && (q[1] > q[3])) {
                  d <- d - 1
                }
                if ((q[2] > q[1]) && (q[2] > q[3]) && (n[d] < Dsize)) {
                  d <- d
                }
                if ((q[3] > q[1]) && (q[3] > q[2])) {
                  d <- d + 1
                }
              }
            }
          }
          total <- sum(n)
          totaltime <- totaltime + max(assessmenttime)
          if ((total >= N.max)) {
            st <- 1
            nodose <- 0
          } #### stop the trial when maximum sample size is reached
          if ((d >= 1) && (d <= D) && (q[2] > q[1]) && (q[2] > q[3]) && (n[d] >= Dsize)) {
            st <- 1
          }
          if (n[d] >= Dsize) {
            st <- 1
          }
        }
        # Storing the time value
        Totaltime[sim] <- totaltime
        # compute the posterior mean from the beta distribution
        if (nodose == 0) {
          tdose <- toxdose - 1
          pp <- rep(-100, tdose)
          pp.var <- rep(0, tdose)
          for (i in 1:tdose) {
            pp[i] <- (x[i] + .005) / (n[i] + .01)
            pp.var[i] <- betavar(x[i] + 0.005, n[i] - x[i] + 0.005) ### here adding 0.005 is to use beta(0.005, 0.005) for estimating the MTD, which is different from the dose-finding.
          }
          pp <- pava(pp, wt = 1 / pp.var) # perform the isotonic transformation using PAVA with weights being posterior variances
          for (i in 2:tdose) {
            pp[i] <- pp[i] + i * 1E-10 # by adding an increasingly small number to toxicity probability at higher doses, it will break the ties and make the lower dose level the MTD if the ties were larger than pT or make the higher dose level the MTD if the ties are smaller than pT
          }
          seldose <- c(1:tdose)[sort(abs(pp[1:tdose] - pT), index.return = T)$ix[1]]
          # seldose is the final MTD that is selected based on the order-transformed posterior means
        }
        res[sim] <- seldose
        for (i in 1:D) {
          datan[sim, i] <- n[i]
          datax[sim, i] <- x[i]
        }
      }
      abc <- rep(0, D)
      for (i in 1:D) {
        abc[i] <- sum(res == i) / simN # abc is the proportion of each dose is selected as the MTD
      }
      
      
      mtpi_8weeks <- data.frame(matrix(ncol = 2, nrow = D+5))
      colnames(mtpi_8weeks) <- c("True Toxicity Probability", "MTPI")
      
      mtpi_8weeks$`True Toxicity Probability` <-c(true_dlt_prob,"Trials Selecting All Doses as Safe (%)",
                                                  "Trials Selecting All Doses as Toxic (%)",
                                                  "Average time taken (In weeks)",
                                                  
                                                  "Average Sample size","Average Overdosed patients (%)")
      
      patients_at_each_dose <- round(apply(datan,2,mean),3)
      
      taget_position <- which(min(abs(true_dlt_prob - pT))==(true_dlt_prob - pT))
      over_dose_per  <- (sum(patients_at_each_dose[taget_position+1:D],na.rm = TRUE)*100)/round(sum(apply(datan, 2, mean)), 3)
      
      new_row <- c( round((abc * 100), 2),round(((alldosesafe / simN) * 100), 2), round(((alldosetoxic / simN) * 100), 2), round(sum(Totaltime / simN), 2), round(sum(apply(datan, 2, mean)), 3),round(over_dose_per,2))
      mtpi_8weeks$MTPI <- new_row
      
      return(mtpi_8weeks)
    }
    
    
    No_of_simulations <- reactive(input$sim)
    
    true_dlt_prob <- reactive(as.numeric(unlist(strsplit(input$dall,","))))
    
    
    
    taget_tox <- reactive(input$pT)
    max_sample <- reactive(input$max_sample)
    max_dose_sample <- reactive(input$max_sample_dose)
    dlt_window <- reactive(input$dtl_window)
    Cohort_size <- reactive(input$cohort)
    enroll_time <- reactive(input$ac)
    no_dose <- reactive(input$no_dose)
    
    
    tite_boin_data <- tite_boin(
      true_dlt_prob = true_dlt_prob(), no_dose = no_dose(), No_of_simulations = No_of_simulations(),
      taget_tox = taget_tox(), Cohort_size = Cohort_size(), dlt_window = dlt_window(),
      max_sample = max_sample(), max_dose_sample = max_dose_sample(), enroll_time = enroll_time()
    )
    
    mtpi_8weeks_data <- mtpi_simulation_8weeks(
      true_dlt_prob = true_dlt_prob(), no_dose = no_dose(), No_of_simulations = No_of_simulations(),
      taget_tox = taget_tox(), Cohort_size = Cohort_size(), dlt_window = dlt_window(),
      max_sample = max_sample(), max_dose_sample = max_dose_sample(), enroll_time = enroll_time()
    )
    
    all_design_data <- merge(mtpi_8weeks_data, tite_boin_data, by = "True Toxicity Probability", sort = FALSE)
    
    mtpi_8weeks_plotdata <- data.frame(`True DLT rate` = head(mtpi_8weeks_data$`True Toxicity Probability`, no_dose()), Design = "MTPI", `MTD Selection` = head(mtpi_8weeks_data$MTPI, no_dose()))
    tite_boin_plotdata <- data.frame(`True DLT rate` = head(tite_boin_data$`True Toxicity Probability`, no_dose()), Design = "TITE BOIN", `MTD Selection` = head(tite_boin_data$`TITE BOIN`, no_dose()))
    plotdata <- rbind(mtpi_8weeks_plotdata, tite_boin_plotdata)
    
    # DATA FOR PLOT
    
    mtpi_8weeks_timedata <- data.frame(Design = "MTPI", mtpi_8weeks_data %>% filter(`True Toxicity Probability` == "Average time taken (In weeks)"))
    tite_boin_timedata <- data.frame(Design = "TITE BOIN", tite_boin_data %>% filter(`True Toxicity Probability` == "Average time taken (In weeks)"))
    names(tite_boin_timedata)[names(tite_boin_timedata) == "TITE.BOIN"] <- "Time"
    names(mtpi_8weeks_timedata)[names(mtpi_8weeks_timedata) == "MTPI"] <- "Time"
    time_plot_data <- rbind(mtpi_8weeks_timedata, tite_boin_timedata)
    
    
    
    high_slow <- plotdata %>%
      hchart("column", hcaes(x = Design, y = MTD.Selection, group = True.DLT.rate)) %>%
      hc_yAxis_multiples(
        list(title = list(text = "OC'S(%)"), opposite = FALSE),
        list(showLastLabel = T, opposite = TRUE, title = list(text = "Average time taken (In weeks)"))
      ) %>%
      hc_add_series(time_plot_data, "line", hcaes(Design, Time, group = True.Toxicity.Probability), yAxis = 1) %>%
      hc_title(text = "<b>Scenario: MTPI TITE-BOIN Simulation Comparision Plot <b>", align = "center") %>%
      hc_add_theme(hc_theme_ffx())
    
    
    
    output$table <- renderTable(all_design_data)
    output$plot <- renderHighchart(high_slow)
  })
  
  output$image_dose <- renderUI({
    imageOutput("image_plot")
  })
  
  output$image_plot <- renderImage({
    list(src = "C:/Users/babla/Documents/www/3 injection.png", height = "700px", width = "2000px",
         contentType = 'image/png')
  }, deleteFile = FALSE) # Do not forget this option
  
  output$text <- renderUI({
    span(textOutput("image_text"),style = "color:black; font-size:30px")
  })
  
  output$image_text <- renderText(
    "Should we change the dose regime?"
    
  )
  
  
  observeEvent(input$Analyze, {     
    
    mtpi_calculation = function(no_pat,no_dlt,e1,e2,alpha,beta,pT){
      x = no_dlt
      n = no_dlt
      # Update prior beta distribution to obtain posterior distribution
      pa <- no_dlt + alpha
      pb <- no_pat - no_dlt + beta
      q = rep(0,3)
      q[1] <- (1 - pbeta(e2 + pT, pa, pb)) / (1 - e1 - pT) ### overdosing interval
      
      q[2] <- (pbeta(e2 + pT, pa, pb) - pbeta(pT - e1, pa, pb)) / (e2 + e1) ### Proper dosing interval
      q[3] <- (pbeta(pT - e1, pa, pb) / (pT - e1))  ### underdosing interval
      if(q[1]>q[2] & q[1]>q[3]){
        results = print("Go to next lower dose level")
      }else  if(q[2]>q[1] & q[2]>q[3]){
        results = print("Stay at the Current dose level")
      }else{
        results = print("Go to next higher dose level")
      }
      
    }
    
    
    design_type = reactive(input$design_type)
    
    pT <- reactive(input$pT)
    no_pat=reactive(input$no_pat)
    dlt_win <- reactive(input$dlt_win)
    e1 <- reactive(input$e1)
    e2 <- reactive(input$e2)
    alpha <- reactive(input$alpha)
    beta <- reactive(input$beta)
    no_dlt <- reactive(input$no_dlt)
    
    results=mtpi_calculation(no_pat = no_pat(),
                             no_dlt = no_dlt(),e1 = e1(),e2 = e2(),alpha = alpha(),beta = beta(),pT = pT())
    
    output$image_text <- renderText(
      results
    )
  
    output$image_plot <- renderImage({
      
      if(results == "Go to next lower dose level") {
        list(src = "C:/Users/babla/Documents/www/lower.png",height = "500px",
             contentType = 'image/png')
      }else  if(results == "Go to next higher dose level") {
        list(src = "C:/Users/babla/Documents/www/high.png",height = "500px",
             contentType = 'image/png')
      }  else {
        list(src = "C:/Users/babla/Documents/www/stay.png",height = "500px",
             contentType = 'image/png')
      }
      
    }, deleteFile = FALSE) # Do not forget this option
    
    
    
    
  })
  observeEvent(input$Run, {     
    
    
    titeboin_calculation = function(no_pat1,no_dlt1,pend,comp,time,dlt_win,pT)
    {
      DLTwindow <- dlt_win
      phi <- pT
      phi_1 <- 0.6 * phi
      phi_2 <- 1.4 * phi
      lambda_e <- log((1 - phi_1) / (1 - phi)) / log((phi * (1 - phi_1)) / (phi_1 * (1 - phi)))
      lambda_d <- log((1 - phi) / (1 - phi_2)) / log((phi_2 * (1 - phi)) / (phi * (1 - phi_2)))
      a <- 0.5 * phi # alpha in paper
      b <- 1 - a # beta in paper
      # estimating Phat
      c <- pend
      n <- no_pat1 # number of subjects
      r <- comp # subjects with Completed data
      TFT = time
      s <- no_dlt1 # DLT count
      p_mean <- (s + a) / (r + a + b) # posterior mean
      STFT = sum(TFT)/dlt_win
      
      # p_hat=0
      if (pend == 0) {
        p_hat <- s / n
      } else {
        p_hat <- (s + (p_mean / (1 - p_mean)) * (c - STFT)) / n
      }
      
      if (p_hat >= lambda_d) {
        results = print("Go to next lower dose level")
      }else if ((p_hat > lambda_e) && (p_hat < lambda_d) && (num[d] < Dsize)) {
        results = print("Stay at the current dose level")
      }else  {
        results = print("Go to next higher dose level")
        
      } 
   
      
    }
    
    
    
    time <- reactive(as.numeric(unlist(strsplit(input$time,","))))
    pT <- reactive(input$pT)
    no_pat1 <- reactive(input$no_pat1)
    pending <- reactive(input$pending)
    no_dlt1 <- reactive(input$no_dlt1)
    comp <- reactive(input$comp)
    dlt_win <- reactive(input$dlt_win)
    
    results=titeboin_calculation(no_pat1 = no_pat1(),
                                 no_dlt1 = no_dlt1(),pend = pending(),comp = comp(),time = time(),dlt_win = dlt_win(),pT = pT())
   
    output$image_text <- renderText(
      results
    )
    
    output$image_plot <- renderImage({
      
      if(results == "Go to next lower dose level") {
        list(src = "C:/Users/babla/Documents/www/lower.png",height = "500px",
             contentType = 'image/png')
      }else  if(results == "Go to next higher dose level") {
        list(src = "C:/Users/babla/Documents/www/high.png",height = "500px",
             contentType = 'image/png')
      }  else {
        list(src = "C:/Users/babla/Documents/www/stay.png",height = "500px",
             contentType = 'image/png')
      }
      
    }, deleteFile = FALSE) # Do not forget this option
    
    
    
    
  })
  
  
})
