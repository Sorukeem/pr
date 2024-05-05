#Royce & Fu replication
install.packages("deSolve")

#define SIR Model function
library(deSolve)
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Differential equations
    dSw <- bw - beta_w * Sw * Iw - mw * Sw
    dIw <- beta_w * Sw * Iw - gamma_w * Iw - mw * Iw
    dRw <- gamma_w * Iw - mw * Rw
    
    dSd <- bd - beta_d * Sd * Id - pd * Sd * Iw - beta_d * Sd * Td - md * Sd
    dId <- beta_d * Sd * Id + pd * Sd * Iw - mu * Id - gamma_d * Id - md * Id
    dTd <- mu * Id + beta_d * Sd * Td - gamma_d * Td - md * Td
    dRd <- gamma_d * Id + gamma_d * Td - md * Rd
    
    dSh <- bh - beta_h * Sh * Ih - ph * Sh * Td - mh * Sh
    dIh <- beta_h * Sh * Ih + ph * Sh * Td - gamma_h * Ih - mh * Ih
    dRh <- gamma_h * Ih - mh * Rh
    
    return(list(c(dSw, dIw, dRw, dSd, dId, dTd, dRd, dSh, dIh, dRh)))
  })
}


#set parameters
initial_state <- c(Sw = 0.5, Iw = 0.5, Rw = 0, Sd = 1, Id = 0, Td = 0, Rd = 0, Sh = 1, Ih = 0, Rh = 0)
parameters <- c(beta_w = 0.3, gamma_w = 0.1, mw = 0.01, 
                beta_d = 0.5, gamma_d = 0.1, md = 0.01, pd = 0.02, mu = 0.02,
                beta_h = 0.4, gamma_h = 0.1, mh = 0.01, ph = 0.02,
                bw = 0.02, bd = 0.02, bh = 0.02)

#solve the differential equations
time <- seq(0, 200, by = 1)  # from 0 to 200 days
output <- ode(y = initial_state, times = time, func = sir_model, parms = parameters)


#plot

matplot(time, output[, c("Sh", "Ih", "Rh")], type = "l", col = c("green", "red", "blue"), lty = 1,
        xlab = "Time (days)", ylab = "Population Proportion",
        main = "SIR Model Dynamics in Human Population")
legend("right", legend = c("Susceptible", "Infected", "Recovered"), col = c("green", "red", "blue"), lty = 1)

library(ggplot2)
time <- seq(0, 200, by = 1)  # from 0 to 200 days
output <- ode(y = initial_state, times = time, func = sir_model, parms = parameters)
output_df <- as.data.frame(output)
output_df$time <- output_df$time
