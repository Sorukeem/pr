#Royce & Fu (2020) replication
#April 24, 2024 (v3)

library(deSolve)
library(ggplot2)

sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS_w <- b_w - (beta_w * S_w * I_w) - (m_w * S_w)
    dI_w <- (beta_w * S_w * I_w) - (gamma_w * I_w) - (m_w * I_w)
    dR_w <- (gamma_w * I_w) - (m_w * R_w)
    
    dS_d <- b_d - (beta_d * S_d * I_d) - (p_d * S_d * I_w) - (beta_d * S_d * T_d) - (m_d * S_d)
    dI_d <- (beta_d * S_d * I_d) + (p_d * S_d * I_w) - (mu * I_d) - (gamma_d * I_d) - (m_d * I_d)
    dT_d <- (mu * I_d) + (beta_d * S_d * T_d) - (gamma_d * T_d) - (m_d * T_d)
    dR_d <- (gamma_d * I_d) + (gamma_d * T_d) - (m_d * R_d)
    
    dS_h <- b_h - (beta_h * S_h * I_h) - (p_h * S_h * T_d) - (m_h * S_h)
    dI_h <- (beta_h * S_h * I_h) + (p_h * S_h * T_d) - (gamma_h * I_h) - (m_h * I_h)
    dR_h <- (gamma_h * I_h) - (m_h * R_h) 
    
    #return the rate of change
    return(list(c(dS_w, dI_w, dR_w, dS_d, dI_d, dT_d, dR_d, dS_h, dI_h, dR_h)))
    })
}

# Parameters as described in the paper
parameters <- c(b_w = 1, # wild birth rate
                m_w = 1, # wild background death rate
                b_d = 1, # domestic birth rate
                m_d = 1, # domestic background death rate
                b_h = 0.00003252, # human birth rate # used daily rate
                m_h = 0.00002477,# human mortality rate # used daily rate
                beta_w = 0.89, # wild transmission rate 
                gamma_w = 0.981, # wild recovery rate 
                beta_d = 0.89, # domestic transmission rate 
                gamma_d = 0.981, # domestic recovery rate
                p_d = 0.51, # wild-domestic transmission rate
                beta_h = 0.078, # human transmission rate
                gamma_h = 0.091, # human recovery rate
                p_h = 0.207, # human-domestic transmission rate
                mu = 0.499 # mutation rate in domestic
                )

# Initial state values
initial_state <- c(S_w = 0.5, # From paper
                   I_w = 0.5, # From paper
                   R_w = 0, # Assumed, not specified in the paper, assuming outbreak onset
                   S_d = 1, # Assumed, not specified
                   I_d = 0, # Assumed, not specified
                   T_d = 0, # Assumed, not specified, no initial transmissible cases assumed (intermediate hosts infected with human-transmissible strain)
                   R_d = 0, # Assumed, not specified, assuming outbreak onset
                   S_h = 1, # Assumed, all humans initially susceptible
                   I_h = 0, # Assumed, no initially infected humans
                   R_h = 0 # Assumed, no initially recovered humans
                   )
# time
times <- seq(0, 1000, by = 1)

# Solve the differential equations
out <- ode(y = initial_state, times = times, func = sir_model, parms = parameters)

# Create a data frame for plotting
out_df <- as.data.frame(out)
out_df$time <- out[, "time"]

#reshape
library(reshape2)
out_long <- melt(out_df, id.vars = "time", variable.name = "Compartment", value.name = "Population")

# Plot
#plot3 - combine ======
library(gridExtra)

#add species group
out_long$SpeciesGroup <- ifelse(grepl("w", out_long$Compartment), "Wild",
                                    ifelse(grepl("d", out_long$Compartment), "Domestic", "Human"))

#add state
out_long$State <- ifelse(grepl("T", out_long$Compartment), "T",
                             substr(out_long$Compartment, 1, 1))

# Define color palette for S, I, T, R
colors <- c("S" = "blue", "I" = "red", "T" = "pink", "R" = "green")

plot_wild <- ggplot(data = subset(out_long, SpeciesGroup == "Wild"), aes(x = time, y = Population, color = State)) +
  geom_line() +
  labs(title = "Wild Hosts", x = "Time (days)", y = "Population Proportion") +
  scale_color_manual(values = colors) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(limits = c(-0.2, 1.2))

plot_domestic <- ggplot(data = subset(out_long, SpeciesGroup == "Domestic"), aes(x = time, y = Population, color = State)) +
  geom_line() +
  labs(title = "Domestic Animals", x = "Time (days)", y = "Population Proportion") +
  scale_color_manual(values = colors) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(limits = c(-0.2, 1.2))

plot_human <- ggplot(data = subset(out_long, SpeciesGroup == "Human"), aes(x = time, y = Population, color = State)) +
  geom_line() +
  labs(title = "Humans", x = "Time (days)", y = "Population Proportion") +
  scale_color_manual(values = colors) +
  theme_minimal() #+
  #scale_x_continuous(limits = c(0, 1000)) +
  #scale_y_continuous(limits = c(0, 1))

# Combine the plots into one figure
combined_plot <- grid.arrange(plot_wild, plot_domestic, plot_human, ncol = 3)

