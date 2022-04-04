rm(list = ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)
library(stargazer)
library(entropy)
library(cowplot)
library(reshape2)
library(networkD3) # for sankey diagram
library(dplyr) # for sankey diagram

# Sankey diagram based on estimated parameters ----

## (Figure 9) Merger simulation across different subsidy threshold and amount. ----

# Sankey diagrams need to be exported from a clipboard.
# For the draft, I added information of x-axis (subsidy threshold and amount) 
# manually to the figures generated here. 
# The data input of simulated merger configurations 
# is computed on ship_merger_counterfactual.jl.

## Configuration of the pre-merger market ----

# A connection data frame is a list of flows with intensity for each flow
links <- 
  data.frame(
    source=c("group_A", "group_A", "group_B",
             "group_C", "group_C", "group_E"), 
    target=c("group_C", "group_D", "group_E",
             "group_F", "group_G", "group_H"), 
    value=c(2, 1, 2, 3, 1, 3)
  )

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <-
  data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- 
  match(links$source, nodes$name)-1 
links$IDtarget <- 
  match(links$target, nodes$name)-1

# Make the Network
p <- 
  networkD3::sankeyNetwork(Links = links,
                           Nodes = nodes,
                           Source = "IDsource",
                           Target = "IDtarget",
                           Value = "value",
                           NodeID = "name", 
                           sinksRight=FALSE)
# save the widget
# library(htmlwidgets)
# saveWidget(p, file=paste0( getwd(), "/HtmlWidget/sankeyBasiGroup12.html"))


# A connection data frame is a list of flows with intensity for each flow
links <- 
  data.frame(
    source=c("Nippon Yusen",
             "Mitsubishi Kaiun",
             "Osaka Syosen", 
             "Mitsui Senpaku",
             "Nitto Syosen",
             "Daido Kaiun",
             "Kawasaki Kisen",
             "Iino Kaiun",
             "Yamashita Kisen",
             "Shin-nihon Kisen",
             "Nihon Yusosen", 
             "Nissan Kisen"), 
    target=c("Nippon Yusen ",
             "Nippon Yusen ", 
             "Mitsui O.S.K. Line",
             "Mitsui O.S.K. Line",
             "Japan Lines",
             "Japan Lines",
             "Kawasaki Kisen Kaisha", 
             "Kawasaki Kisen Kaisha",
             "Yamashita-Shinnihon Kisen",
             "Yamashita-Shinnihon Kisen",
             "Showa Line", "Showa Line"), 
    value=c(1, 1, 1, 1, 1, 1)
  )

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <-
  data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <-
  match(links$source, nodes$name)-1 
links$IDtarget <-
  match(links$target, nodes$name)-1


my_color <- 'd3.scaleOrdinal() .domain(["Nippon Yusen (id 1)","Mitsubishi Kaiun (id 2)",
           "Osaka Syosen (id 3-9)", "Mitsui Senpaku (id 4-10)",
           "Nitto Syosen (id 5-13)", "Daido Kaiun (id 6-14)",
           "Kawasaki Kisen (id 7-15)", "Iino Kaiun (id 8-16)",
           "Yamashita Kisen (id 9-17)", "Shin-nihon Kisen (id 10-18)",
           "Nihon Yusosen (id 11-23)", "Nissan Kisen (id 12-24)"]) .range(["grey", "grey" , "grey", "grey", "grey" , "grey", "grey", "grey" , "grey", "grey", "grey" , "grey"])'


# Make the Network
p <- 
  networkD3::sankeyNetwork(Links = links,
                           Nodes = nodes,
                           Source = "IDsource",
                           Target = "IDtarget",
                           Value = "value",
                           NodeID = "name",
                           fontSize = 15,
                           colourScale = my_color,
                           sinksRight = FALSE)
p


## Changing subsidy amount ----

# A connection data frame is a list of flows with intensity for each flow
links <- 
  data.frame(
    source = c("Nippon Yusen","Mitsubishi Kaiun",
               "Osaka Syosen)", "Mitsui Senpaku",
               "Nitto Syosen", "Daido Kaiun",
               "Kawasaki Kisen", "Iino Kaiun",
               "Yamashita Kisen", "Shin-nihon Kisen",
               "Nihon Yusosen", "Nissan Kisen",
               # 1st layer
               "U1 (0.1)", "U2 (0.1)", 
               "U3 (0.1)", "U4 (0.1)",
               "U5 (0.1)", "U6 (0.1)",
               "U7 (0.1)", "U8 (0.1)",
               "U9 (0.1)", "U10 (0.1)",
               "U11 (0.1)", "U12 (0.1)",
               # 2nd layer
               "U1 (0.25)", "U2 (0.25)", 
               "U3 (0.25)", "U4 (0.25)",
               "U5 (0.25)", "U6 (0.25)",
               "U7 (0.25)", "U8 (0.25)",
               "U9 (0.25)", "U10 (0.25)",
               "U11 (0.25)", "U12 (0.25)",
               # 3rd layer
               "U1 (0.5)", "Group1 (0.5)", 
               "U3 (0.5)", "U4 (0.5)",
               "Group1 (0.5)", "Group2 (0.5)",
               "Group2 (0.5)", "Group4 (0.5)",
               "Group3 (0.5)", "Group3 (0.5)",
               "Group3 (0.5)", "Group4 (0.5)",
               #4th layer
               "U1 (0.6)", "Group1 (0.6)", 
               "U3 (0.6)", "U4 (0.6)",
               "Group1 (0.6)", "Group2 (0.6)",
               "Group2 (0.6)", "Group4 (0.6)",
               "Group3 (0.6)", "Group3 (0.6)",
               "Group3 (0.6)", "Group4 (0.6)",
               #5th layer
               "U1 (0.7)", "Group1 (0.7)", 
               "U3 (0.7)", "U4 (0.7)",
               "Group1 (0.7)", "Group2 (0.7)",
               "Group2 (0.7)", "Group4 (0.7)",
               "Group3 (0.7)", "Group3 (0.7)",
               "Group3 (0.7)", "Group4 (0.7)",
               #6th layer
               "Group5 (0.8)", "Group6 (0.8)", 
               "Group5 (0.8)", "Group8 (0.8)",
               "Group7 (0.8)", "Group8 (0.8)",
               "Group9 (0.8)", "Group10 (0.8)",
               "Group6 (0.8)", "Group7 (0.8)",
               "Group9 (0.8)", "Group10 (0.8)",
               #7th layer
               "Group5 (0.9)", "Group6 (0.9)", 
               "Group5 (0.9)", "Group8 (0.9)",
               "Group7 (0.9)", "Group8 (0.9)",
               "Group9 (0.9)", "Group10 (0.9)",
               "Group6 (0.9)", "Group7 (0.9)",
               "Group9 (0.9)", "Group10 (0.9)",
               #8th layer
               "Group5 (1.0)", "Group6 (1.0)", 
               "Group5 (1.0)", "Group8 (1.0)",
               "Group7 (1.0)", "Group8 (1.0)",
               "Group9 (1.0)", "Group10 (1.0)",
               "Group6 (1.0)", "Group7 (1.0)",
               "Group9 (1.0)", "Group10 (1.0)"), 
    target = c(
      # 1st layer
      "U1 (0.1)", "U2 (0.1)", 
      "U3 (0.1)", "U4 (0.1)",
      "U5 (0.1)", "U6 (0.1)",
      "U7 (0.1)", "U8 (0.1)",
      "U9 (0.1)", "U10 (0.1)",
      "U11 (0.1)", "U12 (0.1)",
      # 2nd layer
      "U1 (0.25)", "U2 (0.25)", 
      "U3 (0.25)", "U4 (0.25)",
      "U5 (0.25)", "U6 (0.25)",
      "U7 (0.25)", "U8 (0.25)",
      "U9 (0.25)", "U10 (0.25)",
      "U11 (0.25)", "U12 (0.25)",
      # 3rd layer
      "U1 (0.5)", "Group1 (0.5)", 
      "U3 (0.5)", "U4 (0.5)",
      "Group1 (0.5)", "Group2 (0.5)",
      "Group2 (0.5)", "Group4 (0.5)",
      "Group3 (0.5)", "Group3 (0.5)",
      "Group3 (0.5)", "Group4 (0.5)",
      #4th layer
      "U1 (0.6)", "Group1 (0.6)", 
      "U3 (0.6)", "U4 (0.6)",
      "Group1 (0.6)", "Group2 (0.6)",
      "Group2 (0.6)", "Group4 (0.6)",
      "Group3 (0.6)", "Group3 (0.6)",
      "Group3 (0.6)", "Group4 (0.6)",
      #5th layer
      "U1 (0.7)", "Group1 (0.7)", 
      "U3 (0.7)", "U4 (0.7)",
      "Group1 (0.7)", "Group2 (0.7)",
      "Group2 (0.7)", "Group4 (0.7)",
      "Group3 (0.7)", "Group3 (0.7)",
      "Group3 (0.7)", "Group4 (0.7)",
      #6th layer
      "Group5 (0.8)", "Group6 (0.8)", 
      "Group5 (0.8)", "Group8 (0.8)",
      "Group7 (0.8)", "Group8 (0.8)",
      "Group9 (0.8)", "Group10 (0.8)",
      "Group6 (0.8)", "Group7 (0.8)",
      "Group9 (0.8)", "Group10 (0.8)",
      #7th layer
      "Group5 (0.9)", "Group6 (0.9)", 
      "Group5 (0.9)", "Group8 (0.9)",
      "Group7 (0.9)", "Group8 (0.9)",
      "Group9 (0.9)", "Group10 (0.9)",
      "Group6 (0.9)", "Group7 (0.9)",
      "Group9 (0.9)", "Group10 (0.9)",
      #8th layer
      "Group5 (1.0)", "Group6 (1.0)", 
      "Group5 (1.0)", "Group8 (1.0)",
      "Group7 (1.0)", "Group8 (1.0)",
      "Group9 (1.0)", "Group10 (1.0)",
      "Group6 (1.0)", "Group7 (1.0)",
      "Group9 (1.0)", "Group10 (1.0)",
      #9th layer
      "Group5 (2.0)", "Group6 (2.0)", 
      "Group5 (2.0)", "Group8 (2.0)",
      "Group7 (2.0)", "Group8 (2.0)",
      "Group9 (2.0)", "Group10 (2.0)",
      "Group6 (2.0)", "Group7 (2.0)",
      "Group9 (2.0)", "Group10 (2.0)"
    ), 
    value = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 1st layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 2nd layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 3rd layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 4th layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 5th layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 6th layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 7th layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 8th layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  )

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- 
  data.frame(
    name = c(as.character(links$source), 
             as.character(links$target)) %>% unique()
  )

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <-
  match(links$source, nodes$name)-1 
links$IDtarget <-
  match(links$target, nodes$name)-1

# Make the Network
p <- 
  networkD3::sankeyNetwork(Links = links,
                           Nodes = nodes,
                           Source = "IDsource",
                           Target = "IDtarget",
                           Value = "value",
                           NodeID = "name", 
                           fontSize = 12,
                           colourScale = my_color,
                           sinksRight = FALSE)
p# size = 1600 * 699

## Changing subsidy threshold ----

# A connection data frame is a list of flows with intensity for each flow
links <-
  data.frame(
    source = c("Nippon Yusen","Mitsubishi Kaiun",
               "Osaka Syosen", "Mitsui Senpaku",
               "Nitto Syosen", "Daido Kaiun",
               "Kawasaki Kisen", "Iino Kaiun",
               "Yamashita Kisen", "Shin-nihon Kisen",
               "Nihon Yusosen", "Nissan Kisen",
               # 1st layer
               "Group1 (1mil)", "Group2 (1mil)", 
               "Group1 (1mil)", "Group4 (1mil)",
               "Group3 (1mil)", "Group4 (1mil)",
               "Group5 (1mil)", "Group6 (1mil)",
               "Group2 (1mil)", "Group3 (1mil)",
               "Group5 (1mil)", "Group6 (1mil)",
               # 2nd layer
               "U (2mil)", "Group8 (2mil)", 
               "Group8 (2mil)", "Group9 (2mil)",
               "Group7 (2mil)", "Group8 (2mil)",
               "Group7 (2mil)", "Group9 (2mil)",
               "Group7 (2mil)", "Group9 (2mil)",
               "Group8 (2mil)", "Group9 (2mil)",
               # 3rd layer
               "Group11 (3mil)","Group10 (3mil)",
               "Group10 (3mil)","Group10 (3mil)",
               "Group10 (3mil)","Group11 (3mil)",
               "Group10 (3mil)", "Group11 (3mil)",
               "Group11 (3mil)", "Group11 (3mil)",
               "Group11 (3mil)", "Group11 (3mil)",
               #4th layer
               "Group13 (4mil)","Group12 (4mil)",
               "Group12 (4mil)","Group13 (4mil)",
               "Group12 (4mil)","Group13 (4mil)",
               "Group13 (4mil)", "Group12 (4mil)",
               "Group13 (4mil)", "Group12 (4mil)",
               "Group13 (4mil)", "Group12 (4mil)",
               #5th layer
               "Group14 (5mil)","Group14 (5mil)",
               "Group14 (5mil)","Group14 (5mil)",
               "Group14 (5mil)","Group14 (5mil)",
               "Group14 (5mil)", "Group14 (5mil)",
               "Group14 (5mil)", "Group14 (5mil)",
               "Group14 (5mil)", "Group14 (5mil)"
    ), 
    target = c(# 1st layer
      # 1st layer
      "Group1 (1mil)", "Group2 (1mil)", 
      "Group1 (1mil)", "Group4 (1mil)",
      "Group3 (1mil)", "Group4 (1mil)",
      "Group5 (1mil)", "Group6 (1mil)",
      "Group2 (1mil)", "Group3 (1mil)",
      "Group5 (1mil)", "Group6 (1mil)",
      # 2nd layer
      "U (2mil)", "Group8 (2mil)", 
      "Group8 (2mil)", "Group9 (2mil)",
      "Group7 (2mil)", "Group8 (2mil)",
      "Group7 (2mil)", "Group9 (2mil)",
      "Group7 (2mil)", "Group9 (2mil)",
      "Group8 (2mil)", "Group9 (2mil)",
      # 3rd layer
      "Group11 (3mil)","Group10 (3mil)",
      "Group10 (3mil)","Group10 (3mil)",
      "Group10 (3mil)","Group11 (3mil)",
      "Group10 (3mil)", "Group11 (3mil)",
      "Group11 (3mil)", "Group11 (3mil)",
      "Group11 (3mil)", "Group11 (3mil)",
      #4th layer
      "Group13 (4mil)","Group12 (4mil)",
      "Group12 (4mil)","Group13 (4mil)",
      "Group12 (4mil)","Group13 (4mil)",
      "Group13 (4mil)", "Group12 (4mil)",
      "Group13 (4mil)", "Group12 (4mil)",
      "Group13 (4mil)", "Group12 (4mil)",
      #5th layer
      "Group14 (5mil)","Group14 (5mil)",
      "Group14 (5mil)","Group14 (5mil)",
      "Group14 (5mil)","Group14 (5mil)",
      "Group14 (5mil)", "Group14 (5mil)",
      "Group14 (5mil)", "Group14 (5mil)",
      "Group14 (5mil)", "Group14 (5mil)",
      #6th layer
      "Group14 (7.5mil)","Group14 (7.5mil)",
      "Group14 (7.5mil)","Group14 (7.5mil)",
      "Group14 (7.5mil)","Group14 (7.5mil)",
      "Group14 (7.5mil)", "Group14 (7.5mil)",
      "Group14 (7.5mil)", "Group14 (7.5mil)",
      "Group14 (7.5mil)", "Group14 (7.5mil)"
    ), 
    value = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 1st layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 2nd layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 3rd layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 4th layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              # 5th layer
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  )

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- 
  data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <-
  match(links$source, nodes$name)-1 
links$IDtarget <-
  match(links$target, nodes$name)-1

# Make the Network
p <- 
  networkD3::sankeyNetwork(Links = links,
                           Nodes = nodes,
                           Source = "IDsource",
                           Target = "IDtarget",
                           Value = "value",
                           NodeID = "name", 
                           fontSize = 12,
                           colourScale = my_color,
                           sinksRight = FALSE)
p # size = 1600 * 699

