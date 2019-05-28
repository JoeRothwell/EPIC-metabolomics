# From https://mikeyharper.uk/flowcharts-in-r-using-diagrammer/
# A minimal plot
grViz("digraph {
  
  graph[layout = dot, rankdir = LR]

  this
  is
  flow

  this -> is -> flow
}"
)

# Define some sample data
data <- list(a=nrow(ints), b=800, c=600, d=400)


DiagrammeR::grViz("
digraph graph2 {

graph [layout = dot]

# node definitions with substituted label text
node [shape = rectangle, width = 4, fillcolor = Biege]
a [label = '@@1']
b [label = '@@2']
c [label = '@@3']
d [label = '@@4']

a -> b -> c -> d

}

[1]:  paste0('All samples (n = ', data$a, ')')
[2]: paste0('Remove Errors (n = ', data$b, ')')
[3]: paste0('Identify Potential Customers (n = ', data$c, ')')
[4]: paste0('Select Top Priorities (n = ', data$d, ')')
")