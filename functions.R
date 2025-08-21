# the main function to get the cumulative probability of late-stage cancer
get.prob <-
  function(test.t = c(0:3), # test time
           fu.tt = c(4:20), # follow-up analysis time
           se.e = 0.8, # sensitivity of detecting early-stage cancer
           se.l = 0.9, # sensitivity of detecting late-stage cancer
           w = 0.00231, # incident rate of preclinical cancer
           la = c(0.152, 0.164, 0.764) # lambda_1, lambda_2, lambda_3,
           )
    
  {
    
    # vector of t with tt(-1)=-inf; tt(0)=0, etc. 
    tt <- function(i) {
      c(-Inf, test.t, fu.t)[i + 2]
    }
    
    # late preclinical cancer prob. 'PQ';
    # when m=n=r, it's P_r(t_r);
    # when n=i;r=r, it's P_i(t_i)Q_i(t_r-t_i)
    late.PQ <- function(n = r, #input for EP-LP, EP-EC
                        m = r #input for Q^{LP-LC}(m-x)
    ) {
      # r=j
      term.A <- NA # term A 
      
      term.A <-  w * la[2] / (la[1] + la[2]) *
        exp(-la[3] * tt(m))  * #Q^{LP,LC}(t(m))
        ((exp(la[3] * tt(n)) - exp(la[3] * tt(n - 1))) / la[3] -
           exp((la[1] + la[2]) * tt(n - 1)) /
           (la[3] - la[1] - la[2]) *
           (exp((la[3] - la[2] - la[1]) * tt(n)) -
              exp((la[3] - la[2] - la[1]) * tt(n - 1))))
      
      # r>0; j<r
      term.B <- 0 # term B 
      if (n > 0) {
        term.B.j <- 0
        for (j in 0:(n - 1))
        {
          term.B.j <- term.B.j +
            (exp((la[1] + la[2]) * tt(j)) -
               exp((la[1] + la[2]) * tt(j - 1))) * (1 - se.e) ^ (n - j)
        }
        
        
        term.B <- w * la[2] / (la[1] + la[2]) *
          exp(-la[3] * tt(m)) / #Q^{LP,LC}(t(m))
          (la[3] - la[2] - la[1]) *
          (exp((la[3] - la[2] - la[1]) * tt(n)) -
             exp((la[3] - la[2] - la[1]) * tt(n - 1))) * term.B.j
        
      }
      
      out <- term.A + term.B
      return(out)
    }
    
    # prob. of late-stage cancer detected at the r^th test: D^{L}_r
    # r starts from 0.
    prob.l.d.r <- function(r) {
      term.1 <- late.PQ(n = r, m = r)
      term.2 <- 0
      if (r > 0) {
        for (i in 0:(r - 1)) {
          term.2 <- term.2 + late.PQ(n = i, m = r) * (1 - se.l) ^ {
            r - i
          }
        }
      }
      out <- se.l * (term.1 + term.2)
      return(out)
    }
    
    
    # prob. of late-stage interval cancer in r^th interval: 
    # int_{t_{r-1}}^{t_{r}} I^{L}_r
    # r=1,...,n+1
    prob.l.inter.r <-
      function(r) {
        term.C <- 0
        term.D.1 <- 0
        term.D.2 <- 0
        
        
        # term C
        for (i in 0:(r - 1)) {
          term.C <- term.C +
            (1 - se.l) ^ (r - i) * (exp(la[3] * (tt(i) - tt(r - 1))) - 
                                      exp(la[3] * (tt(i) - tt(r)))) * 
            late.PQ(n = i, m = i)
        }
        
        
        # term D_1
        
        term.D.1 <- w * la[2] / (la[1] + la[2]) *
          (tt(r) - tt(r - 1) +
             (exp(-la[3] * (tt(
               r
             ) - tt(
               r - 1
             ))) - 1) / la[3] +
             la[3] / (la[3] - la[2] - la[1]) * ((exp(
               -(la[2] + la[1]) * (tt(r) - tt(r - 1))
             ) - 1) / (la[1] + la[2]) -
               (exp(-la[3] * (
                 tt(r) - tt(r - 1)
               )) - 1) / la[3]))
        
        
        # term D_2
        for (j in 0:(r - 1)) {
          term.D.2 <- term.D.2 +
            (1 - se.e) ^ (r - j) *
            w * la[2] / (la[1] + la[2]) *
            (exp((la[1] + la[2]) * tt(j)) - exp((la[1] + la[2]) * tt(j - 1))) *
            la[3] / (la[3] - la[2] - la[1]) *
            ((exp(-(la[2] + la[1]) * tt(r - 1)) -
                exp(-(la[2] + la[1]) * tt(r))) / (la[2] + la[1]) +
               1 / la[3] * (exp(-la[3] * tt(r)) - exp(-la[3] * tt(r - 1))) * 
               exp((la[3] - la[2] - la[1]) * tt(r - 1)))
          
        }
        out <- term.C + term.D.1 + term.D.2
        
        return(out)
      }
    


  # prob. of early-stage cancer detected at r^th test; r starts from 0.
  prob.e.d.r <-
    function(r) {

      
      term.1 <-
        w / (la[1] + la[2]) * (1 - exp(-(la[1] + la[2]) * (tt(r) - tt(r - 1))))
      term.2 <- 0
      if (r > 0) {
        for (i in 0:(r - 1)) {
          term.2 <- term.2 +
            w / (la[1] + la[2]) * exp(-(la[1] + la[2]) * tt(r)) *
            (exp((la[1] + la[2]) * tt(i)) - exp((la[1] + la[2]) * tt(i - 1))) *
            (1 - se.e) ^ {
              r - i
            }
        }
      }
      out <- se.e * (term.1 + term.2)
      return(out)
    }
  
  #prob. of early-stage interval cancer in r^th interval: t_r-1, t_r; r=1,...,n+1
  prob.e.inter.r <-
    function(r) {
     
      term.A <- 0
      for (i in 0:(r - 1)) {
        term.A <- term.A + (1 - se.e) ^ {
          r - i
        } *
          w * la[1] / (la[1] + la[2]) ^ 2 *
          (exp(-(la[1] + la[2]) * tt(r - 1)) - exp(-(la[1] + la[2]) * tt(r))) *
          (exp((la[1] + la[2]) * tt(i)) - exp((la[1] + la[2]) * tt(i - 1)))
      }
      
      
      term.B <- w * la[1] / (la[1] + la[2]) *
        (tt(r) - tt(r - 1) - 1 / (la[1] + la[2]) * (1 - exp(-(la[1] + la[2]) *
                                                              (tt(
                                                                r
                                                              ) - tt(
                                                                r - 1
                                                              )))))
      
      out <- term.A + term.B
      return(out)
    }

  test.fu.t <- c(test.t, fu.tt)
  
  # cumulative screen-detected
  d.late.p <- rep(NA, length(test.fu.t))
  d.early.p <- rep(NA, length(test.fu.t))
  
  # cumulative interval cancer
  int.late.p <- c(0, rep(NA, length(test.fu.t) - 1))
  int.early.p <- c(0, rep(NA, length(test.fu.t) - 1))
  
  # prevalence screening
  fu.t = NA # dummy
  
  d.late.p[1] <- prob.l.d.r(r = 0) # the index starts at 1 for time 0
  d.early.p[1] <- prob.e.d.r(r = 0) # early
  
  
  
  if (length(test.t) > 1) {
    for (j in 1:(length(test.t) - 1)) {
      d.late.p[j + 1] <- d.late.p[j] + # cumulative sum
        prob.l.d.r(r = j) # j^th screening at t_j
      d.early.p[j + 1] <- d.early.p[j] + prob.e.d.r(r = j) # early
      
      int.late.p[j + 1] <- int.late.p[j] +
        prob.l.inter.r(r = j) # interval cancer between t_{j-1} and t_{j}
      
      int.early.p[j + 1] <- int.early.p[j] +
        prob.e.inter.r(r = j)
    }
  }
  
  # follow-up phase
  for (i in (length(test.t) + 1):length(test.fu.t)) {
    fu.t <- test.fu.t[i]
    
    d.late.p[i] <- d.late.p[length(test.t)] # just carry over
    d.early.p[i] <- d.early.p[length(test.t)]
    
    int.late.p[i] <-
      int.late.p[length(test.t)] + # cumulative interval cancer incidence at last screening
      prob.l.inter.r(r = length(test.t)) # since last screening
    
    int.early.p[i] <- int.early.p[length(test.t)] +
      prob.e.inter.r(r = length(test.t))
    
  }
  
  test.time <- c(rep(1, length(test.t)), rep(0, length(fu.tt)))
  dat <- data.frame(time = test.fu.t,
                    test.time,
                    d.early.p,
                    d.late.p,
                    int.early.p,
                    int.late.p)
  
  dat$all.late <- dat$d.late.p + dat$int.late.p
  dat$all.early <- dat$d.early.p + dat$int.early.p
  
  dat$all.detect <- dat$d.late.p + dat$d.early.p
  dat$all.interval <- dat$int.late.p + dat$int.early.p
  
  return(dat)
}





# with step-size s, 
# fill in all the time points between each test/follow-up time points
get.prob.s <- function(test.fu.tt = c(0:3, 20),
                       s = 0.01, # step-size
                       se.e = 0.8,
                       se.l = 0.9,
                       w = 0.00231,
                       la = c(0.152, 0.164, 0.764))
{
  dat <-
    get.prob(
      test.t = test.fu.tt[1],
      fu.tt = seq(
        from = test.fu.tt[1] + s,
        to = test.fu.tt[1 +
                          1] - s,
        by = s
      ),
      se.e = se.e,
      se.l = se.l,
      w = w,
      la = la
    )
  if (length(test.fu.tt) > 2) {
    for (i in 2:c(length(test.fu.tt) - 1))
    {
      dat <- rbind(
        dat,
        get.prob(
          test.t = test.fu.tt[1:i],
          fu.tt = seq(
            from = test.fu.tt[i] + s,
            to = test.fu.tt[i + 1] - s,
            by = s
          ),
          se.e = se.e,
          se.l = se.l,
          w = w,
          la = la
        )
      )
    }
  }
  
  dat <- dat[!duplicated(dat),]
  dat <- dat[order(dat$time),]
  
  return(dat)
  
}


# averaged probability due to staggered entry
acurral.ave <- function(prob, a.yr, s=0.01) {
  # assume uniform enrollment; costumed weight may be modified here
  acc.w <- rep(1, a.yr / s) 
  acc.w <- acc.w / sum(acc.w)
  # padded with 0s for 'filter' type input for convolution
  x <- c(rep(0, length(acc.w) - 1), 
         prob)
  # weighted average as convolution
  w.prob <- convolve(x, acc.w, type = 'f') 
  return(w.prob)
}

# power calculation based on the normal approximation
power.appr <- function(p1, p0, n1, n0, alpha = 0.05) {
  test.stats <-
    (p0 - p1) / sqrt(p0 * (1 - p0) / n0 + p1 * (1 - p1) / n1)
  pnorm(q = qnorm(1 - alpha) - test.stats, lower.tail = F)
}