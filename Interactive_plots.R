library(plotly)

# Example
p <- plot_ly(midwest, x = ~percollege, color = ~state, type = "box")
p

# Manhattan (requires output of Intake_correlation_new.R)
p <- plot_ly(
  df, x = ~RT, y = ~abs(Pcor), 
  text = ~paste("m/z: ", Mass, "RT: ", RT),
  color = df$signif2, colors = "Set2")

