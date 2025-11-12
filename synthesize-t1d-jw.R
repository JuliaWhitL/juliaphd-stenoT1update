library(igraph)
library(TruncatedNormal)


pa <- function(dag, v) {
  neighbors(dag, v, mode = "in")
}

ch <- function(dag, v) {
  neighbors(dag, v, mode = "out")
}

exog <- function(dag) {
  lapply(adjacent_vertices(dag, V(dag), mode = "in"), \(d) length(d) == 0)
}

## EXAMPLE 1 ---
# 3 timepoints w/o competing risks: keeping sex, diabetes duration & statin use

dag2 <- graph_from_literal(
  sex_male,
  diabetes_duration +- sex_male,
  time_cens +- sex_male:diabetes_duration,
  LDL_0 +- sex_male:diabetes_duration,
  LDL_1 +- sex_male:LDL_0:diabetes_duration:statins_0,
  LDL_2 +- sex_male:LDL_1:diabetes_duration:statins_1,
  statins_0 +- LDL_0:diabetes_duration,
  statins_1 +- statins_0:diabetes_duration:LDL_1,
  statins_2 +- statins_1:diabetes_duration:LDL_2,
  time_cvd +- sex_male:diabetes_duration:LDL_0:LDL_1:LDL_2
  # discretize here? (see Keogh, et al) - FIXIT...violation of temporal ordering below (11/10/25)
  
  
)

plot(dag2)
exog(dag2)

# coefficients <- list(
#   sex_male = c(), 
#   statins_0 = c(diabetes_duration = 1),
#   LDL_0 = c(statins_0 = -2, sex_male = 1.25),
#   diabetes_duration = c(sex_male = 0.5),
#   CVD_1 = c(sex_male = 0.5, LDL_0 = 0.3, statins_0 = 1, diabetes_duration = 0.5),
#   statins_1 = c(diabetes_duration = 1.1, statins_0 = 0.5),
#   LDL_1 = c(statins_1 = -2, sex_male = 1.25),
#   CVD_2 = c(sex_male = 0.25, LDL_1 = 0.4, statins_1 = 1, diabetes_duration = 0.2),
#   statins_2 = c(diabetes_duration = 1.2, statins_1 = 0.5),
#   LDL_2 = c(statins_2 = -2, sex_male = 1.25),
#   CVD_3 = c(sex_male = 0.25, LDL_2 = 0.5, statins_2 = 1, diabetes_duration = 0.2)
#   
# )

#ref below
coefficients <- list(
  sex_male = c(),
  diabetes_duration = c(sex_male = 0.5),
  time_cens  = c( sex_male = -0.3, diabetes_duration = -0.05),
  LDL_0 = c(sex_male = 1.25, diabetes_duration = 0.15),
  LDL_1 = c(sex_male = 1.25, LDL_0 = 1, diabetes_duration = 0.15, statins_0 = -2), # clinical point: unsure how we should treat (potential) cumulative effect of statins, DM duration on LDL
  LDL_2 = c(sex_male = 1.25, LDL_1 = 1.1, diabetes_duration = 0.17, statins_1 = -2), # clinical point: assuming statin effect is stable
  statins_0 = c(LDL_0 = 1.3, diabetes_duration = 1), # Q: should we be looking at statins at every timepoint? (see Keogh, et al where they did)
  statins_1 = c(statins_0 = 0.5, LDL_1 = 1.3, diabetes_duration = 1.1),
  statins_2 = c(statins_1 = 0.5, LDL_2 = 1.3, diabetes_duration = 1.2),
  time_cvd = c(sex_male = 0.25, diabetes_duration = 0.2, LDL_0 = 0.3, LDL_1 = 0.4, LDL_2 = 0.5)
)


generate_data <- function(n = 3000, coefficients, dag, intervene = NULL) {
  
  dag_order <- topo_sort(dag, mode = "out")
  thisdata <- data.frame(dummy = rep(1, n))
  
  generator <- list(

    sex_male = function(thisdata) {
      rbinom(n, 1, .45)
    },
    
    diabetes_duration = function(thisdata) {
      
      lpred <- c(thisdata[,pa(dag, "diabetes_duration")$name] * coefficients[["diabetes_duration"]][pa(dag, "diabetes_duration")$name])
      sapply(lpred, \(lp) rtnorm(1, lp + 42, 7.3, 18, 100))
      
    },
    
    time_cens = function(thisdata) {
      pas <- pa(dag, "time_cens")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$time_cens[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(exp(gamma + alpha)) - 15
      }, interval = c(-50, 50))$root
      
      rexp(n, 1 / exp(alpha0 + gamma))
      
    },
    
    LDL_0 = function(thisdata) {
      pas <- pa(dag, "LDL_0")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$LDL_0[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 2.5
      }, interval = c(-500, 500))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 1.25, 0, 15))
      
    },
    
    LDL_1 = function(thisdata) {
      pas <- pa(dag, "LDL_1")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$LDL_1[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 3.15
      }, interval = c(-500, 500))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 1.25, 0, 15))
      
    },
    
    LDL_2 = function(thisdata) {
      pas <- pa(dag, "LDL_2")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$LDL_2[pas])
      alpha0 <- uniroot(function(alpha) {
        mean((gamma + alpha)) - 3.15
      }, interval = c(-500, 500))$root
      sapply(alpha0 + gamma, \(lp) rtnorm(1, lp, 1.25, 0, 15))
    },
    
    statins_0 = function(thisdata) {
      pas <- pa(dag, "statins_0")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$statins_0[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(plogis(gamma + alpha)) - .25
      }, interval = c(-500, 500))$root
      rbinom(n, 1, plogis(alpha0 + gamma))
      
    },

    statins_1 = function(thisdata) {
      pas <- pa(dag, "statins_1")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$statins_1[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(plogis(gamma + alpha)) - .25
      }, interval = c(-500, 500))$root
      rbinom(n, 1, plogis(alpha0 + gamma))
    },
    
    statins_2 = function(thisdata) {
      pas <- pa(dag, "statins_2")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$statins_2[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(plogis(gamma + alpha)) - .25
      }, interval = c(-500, 500))$root
      rbinom(n, 1, plogis(alpha0 + gamma))
    },
    

    
    time_cvd = function(thisdata) {
      pas <- pa(dag, "time_cvd")$name
      gamma <- c(as.matrix(thisdata[, pas]) %*% coefficients$time_cvd[pas] )
      alpha0 <- uniroot(function(alpha) {
        mean(exp(gamma + alpha)) - 12
      }, interval = c(-5000, 5000))$root
      
      rweibull(n, 5, 1 / exp(alpha0 + gamma))
      
    }
    
  )
  
  for(vs in dag_order$name) {
    if(!is.null(intervene)) {
      if(vs %in% names(intervene)) {
        thisdata[[vs]] <- intervene[[vs]]
        next
      }
    }
    
    thisdata[[vs]] <- generator[[vs]](thisdata)
  }
  
  thisdata$dummy <- NULL
  thisdata
  
}

set.seed(120825)
dta <- generate_data(1000, coefficients = coefficients, dag = dag2)

summary(dta)
hist(dta$time_cvd)
mean(dta$time_cvd < .1)
require(Hmisc)
html(describe(dta))

# DAG to use for sequential prediction
seq_dag <- graph_from_literal(
  LDL_0 -+ LDL_1, LDL_1 -+ LDL_2, 
  LDL_0 -+ statins_0, LDL_1 -+ statins_1,
  LDL_2 -+ statins_2, statins_0 -+ CVD_1,
  statins_1 -+ CVD_2, statins_2 -+ CVD_3,
  sex -+ diabetes_duration, sex -+ CVD_1, 
  sex -+ CVD_2, sex -+ CVD_3,
  sex -+ LDL_0, sex -+ LDL_1, sex -+ LDL_2,
  statins_0 -+ statins_1, statins_1 -+ statins_2,
  statins_0 -+ LDL_1, statins_1 -+ LDL_2, 
  CVD_1 -+ statins_1, CVD_2 -+ statins_2,
  diabetes_duration -+ CVD_1, diabetes_duration -+ CVD_2, diabetes_duration -+ CVD_3,
  simplify = FALSE
)

layout_matrix <- matrix(c(
  # x, y coordinates
  -2, -2,    # sex
  -2, -3,    # diabetes_duration
  -1, 2,    # LDL_0
  0, 2,     # LDL_1
  1, 2,     # LDL_2
  -1, 0,    # statins_0
  0, 0,     # statins_1
  1, 0,     # statins_2
  -1, -2,   # CVD_1
  0, -2,    # CVD_2
  1, -2     # CVD_3
), ncol = 2, byrow = TRUE)

rownames(layout_matrix) <- c("sex", "diabetes_duration", "LDL_0", "LDL_1", 
                             "LDL_2", "statins_0", "statins_1", "statins_2",
                             "CVD_1", "CVD_2", "CVD_3")

layout_matrix <- layout_matrix[V(seq_dag)$name, ]

plot(seq_dag, layout = layout_matrix, 
     vertex.size = 30,
     vertex.label.cex = 0.8,
     edge.arrow.size = 0.5)
