rm(list = ls())

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
    source=c("Nippon Yusen (id 1)",
             "Mitsubishi Kaiun (id 2)",
             "Osaka Syosen (id 3-9)", 
             "Mitsui Senpaku (id 4-10)",
             "Nitto Syosen (id 5-13)",
             "Daido Kaiun (id 6-14)",
             "Kawasaki Kisen (id 7-15)",
             "Iino Kaiun (id 8-16)",
             "Yamashita Kisen (id 9-17)",
             "Shin-nihon Kisen (id 10-18)",
             "Nihon Yusosen (id 11-23)", 
             "Nissan Kisen (id 12-24)"), 
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
                           colourScale = my_color,
                           sinksRight = FALSE)
p


## Changing subsidy amount ----

# A connection data frame is a list of flows with intensity for each flow
links <- 
  data.frame(
    source = c("Nippon Yusen (id 1)","Mitsubishi Kaiun (id 2)",
               "Osaka Syosen (id 3-9)", "Mitsui Senpaku (id 4-10)",
               "Nitto Syosen (id 5-13)", "Daido Kaiun (id 6-14)",
               "Kawasaki Kisen (id 7-15)", "Iino Kaiun (id 8-16)",
               "Yamashita Kisen (id 9-17)", "Shin-nihon Kisen (id 10-18)",
               "Nihon Yusosen (id 11-23)", "Nissan Kisen (id 12-24)",
               # 1st layer
               "Group1 (amount=0.1)", "Group1 (amount=0.1)", 
               "Group1 (amount=0.1)", "Group1 (amount=0.1)",
               "Group1 (amount=0.1)", "Group1 (amount=0.1)",
               "Group1 (amount=0.1)", "Group1 (amount=0.1)",
               "Group1 (amount=0.1)", "Group1 (amount=0.1)",
               "Group2 (amount=0.1)", "Group1 (amount=0.1)",
               # 2nd layer
               "Group1 (amount=0.25)", "Group1 (amount=0.25)", 
               "Group1 (amount=0.25)", "Group1 (amount=0.25)",
               "Group1 (amount=0.25)", "Group1 (amount=0.25)",
               "Group1 (amount=0.25)", "Group1 (amount=0.25)",
               "Group1 (amount=0.25)", "Group1 (amount=0.25)",
               "Group2 (amount=0.25)", "Group1 (amount=0.25)",
               # 3rd layer
               "unmatched1 (amount=0.5)", "Group1 (amount=0.5)", 
               "Group2 (amount=0.5)", "unmatched4 (amount=0.5)",
               "Group3 (amount=0.5)", "Group2 (amount=0.5)",
               "Group4 (amount=0.5)", "Group5 (amount=0.5)",
               "Group1 (amount=0.5)", "Group3 (amount=0.5)",
               "Group4 (amount=0.5)", "Group5 (amount=0.5)",
               #4th layer
               "Group6 (amount=0.75)", "Group1 (amount=0.75)", 
               "Group2 (amount=0.75)", "Group6 (amount=0.75)",
               "Group3 (amount=0.75)", "Group2 (amount=0.75)",
               "Group4 (amount=0.75)", "Group5 (amount=0.75)",
               "Group1 (amount=0.75)", "Group3 (amount=0.75)",
               "Group4 (amount=0.75)", "Group5 (amount=0.75)",
               #5th layer
               "Group6 (amount=1.0)", "Group1 (amount=1.0)", 
               "Group2 (amount=1.0)", "Group6 (amount=1.0)",
               "Group3 (amount=1.0)", "Group2 (amount=1.0)",
               "Group4 (amount=1.0)", "Group5 (amount=1.0)",
               "Group1 (amount=1.0)", "Group3 (amount=1.0)",
               "Group4 (amount=1.0)", "Group5 (amount=1.0)",
               #6th layer
               "Group6 (amount=2.0)", "Group1 (amount=2.0)", 
               "Group8 (amount=2.0)", "Group6 (amount=2.0)",
               "Group3 (amount=2.0)", "Group7 (amount=2.0)",
               "Group7 (amount=2.0)", "Group5 (amount=2.0)",
               "Group1 (amount=2.0)", "Group3 (amount=2.0)",
               "Group8 (amount=2.0)", "Group5 (amount=2.0)"), 
    target = c(# 1st layer
      "Group1 (amount=0.1)", "Group1 (amount=0.1)", 
      "Group1 (amount=0.1)", "Group1 (amount=0.1)",
      "Group1 (amount=0.1)", "Group1 (amount=0.1)",
      "Group1 (amount=0.1)", "Group1 (amount=0.1)",
      "Group1 (amount=0.1)", "Group1 (amount=0.1)",
      "Group2 (amount=0.1)", "Group1 (amount=0.1)",
      # 2nd layer
      "Group1 (amount=0.25)", "Group1 (amount=0.25)", 
      "Group1 (amount=0.25)", "Group1 (amount=0.25)",
      "Group1 (amount=0.25)", "Group1 (amount=0.25)",
      "Group1 (amount=0.25)", "Group1 (amount=0.25)",
      "Group1 (amount=0.25)", "Group1 (amount=0.25)",
      "Group2 (amount=0.25)", "Group1 (amount=0.25)",
      # 3rd layer
      "unmatched1 (amount=0.5)", "Group1 (amount=0.5)", 
      "Group2 (amount=0.5)", "unmatched4 (amount=0.5)",
      "Group3 (amount=0.5)", "Group2 (amount=0.5)",
      "Group4 (amount=0.5)", "Group5 (amount=0.5)",
      "Group1 (amount=0.5)", "Group3 (amount=0.5)",
      "Group4 (amount=0.5)", "Group5 (amount=0.5)",
      #4th layer
      "Group6 (amount=0.75)", "Group1 (amount=0.75)", 
      "Group2 (amount=0.75)", "Group6 (amount=0.75)",
      "Group3 (amount=0.75)", "Group2 (amount=0.75)",
      "Group4 (amount=0.75)", "Group5 (amount=0.75)",
      "Group1 (amount=0.75)", "Group3 (amount=0.75)",
      "Group4 (amount=0.75)", "Group5 (amount=0.75)",
      #5th layer
      "Group6 (amount=1.0)", "Group1 (amount=1.0)", 
      "Group2 (amount=1.0)", "Group6 (amount=1.0)",
      "Group3 (amount=1.0)", "Group2 (amount=1.0)",
      "Group4 (amount=1.0)", "Group5 (amount=1.0)",
      "Group1 (amount=1.0)", "Group3 (amount=1.0)",
      "Group4 (amount=1.0)", "Group5 (amount=1.0)",
      #6th layer
      "Group6 (amount=2.0)", "Group1 (amount=2.0)", 
      "Group8 (amount=2.0)", "Group6 (amount=2.0)",
      "Group3 (amount=2.0)", "Group7 (amount=2.0)",
      "Group7 (amount=2.0)", "Group5 (amount=2.0)",
      "Group1 (amount=2.0)", "Group3 (amount=2.0)",
      "Group8 (amount=2.0)", "Group5 (amount=2.0)",
      #7th layer
      "Group6 (amount=3.0)", "Group1 (amount=3.0)", 
      "Group8 (amount=3.0)", "Group6 (amount=3.0)",
      "Group3 (amount=3.0)", "Group7 (amount=3.0)",
      "Group7 (amount=3.0)", "Group5 (amount=3.0)",
      "Group1 (amount=3.0)", "Group3 (amount=3.0)",
      "Group8 (amount=3.0)", "Group5 (amount=3.0)"
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
                           colourScale = my_color,
                           sinksRight = FALSE)
p

## Changing subsidy threshold ----

# A connection data frame is a list of flows with intensity for each flow
links <-
  data.frame(
    source = c("Nippon Yusen (id 1)","Mitsubishi Kaiun (id 2)",
               "Osaka Syosen (id 3-9)", "Mitsui Senpaku (id 4-10)",
               "Nitto Syosen (id 5-13)", "Daido Kaiun (id 6-14)",
               "Kawasaki Kisen (id 7-15)", "Iino Kaiun (id 8-16)",
               "Yamashita Kisen (id 9-17)", "Shin-nihon Kisen (id 10-18)",
               "Nihon Yusosen (id 11-23)", "Nissan Kisen (id 12-24)",
               # 1st layer
               "Group1 (threshold=1mil)", "Group4 (threshold=1mil)", 
               "Group3 (threshold=1mil)", "Group1 (threshold=1mil)",
               "Group2 (threshold=1mil)", "Group3 (threshold=1mil)",
               "Group5 (threshold=1mil)", "Group6 (threshold=1mil)",
               "Group4 (threshold=1mil)", "Group2 (threshold=1mil)",
               "Group5 (threshold=1mil)", "Group6 (threshold=1mil)",
               # 2nd layer
               "Group1 (threshold=2mil)", "Group8 (threshold=2mil)", 
               "Group8 (threshold=2mil)", "Group1 (threshold=2mil)",
               "Group7 (threshold=2mil)", "Group8 (threshold=2mil)",
               "Group7 (threshold=2mil)", "Group1 (threshold=2mil)",
               "unmatched (threshold=2mil)", "Group7 (threshold=2mil)",
               "Group8 (threshold=2mil)", "Group7 (threshold=2mil)",
               # 3rd layer
               "Group1 (threshold=3mil)","Group8 (threshold=3mil)",
               "Group8 (threshold=3mil)","Group9 (threshold=3mil)",
               "Group9 (threshold=3mil)","Group8 (threshold=3mil)",
               "Group1 (threshold=3mil)", "Group1 (threshold=3mil)",
               "Group9 (threshold=3mil)", "Group1 (threshold=3mil)",
               "Group8 (threshold=3mil)", "Group9 (threshold=3mil)",
               #4th layer
               "Group10 (threshold=4mil)","Group11 (threshold=4mil)",
               "Group10 (threshold=4mil)","Group10 (threshold=4mil)",
               "Group11 (threshold=4mil)","Group11 (threshold=4mil)",
               "Group10 (threshold=4mil)", "Group11 (threshold=4mil)",
               "Group10 (threshold=4mil)", "Group11 (threshold=4mil)",
               "Group10 (threshold=4mil)", "unmatched12 (threshold=4mil)",
               #5th layer
               "Group12 (threshold=5mil)","Group12 (threshold=5mil)",
               "Group12 (threshold=5mil)","Group12 (threshold=5mil)",
               "Group11 (threshold=5mil)","unmatched5 (threshold=5mil)",
               "Group12 (threshold=5mil)", "Group12 (threshold=5mil)",
               "Group12 (threshold=5mil)", "Group12 (threshold=5mil)",
               "unmatched11 (threshold=5mil)", "unmatched12 (threshold=5mil)"
    ), 
    target = c(# 1st layer
      "Group1 (threshold=1mil)", "Group4 (threshold=1mil)", 
      "Group3 (threshold=1mil)", "Group1 (threshold=1mil)",
      "Group2 (threshold=1mil)", "Group3 (threshold=1mil)",
      "Group5 (threshold=1mil)", "Group6 (threshold=1mil)",
      "Group4 (threshold=1mil)", "Group2 (threshold=1mil)",
      "Group5 (threshold=1mil)", "Group6 (threshold=1mil)",
      # 2nd layer
      "Group1 (threshold=2mil)", "Group8 (threshold=2mil)", 
      "Group8 (threshold=2mil)", "Group1 (threshold=2mil)",
      "Group7 (threshold=2mil)", "Group8 (threshold=2mil)",
      "Group7 (threshold=2mil)", "Group1 (threshold=2mil)",
      "unmatched (threshold=2mil)", "Group7 (threshold=2mil)",
      "Group8 (threshold=2mil)", "Group7 (threshold=2mil)",
      # 3rd layer
      "Group1 (threshold=3mil)","Group8 (threshold=3mil)",
      "Group8 (threshold=3mil)","Group9 (threshold=3mil)",
      "Group9 (threshold=3mil)","Group8 (threshold=3mil)",
      "Group1 (threshold=3mil)", "Group1 (threshold=3mil)",
      "Group9 (threshold=3mil)", "Group1 (threshold=3mil)",
      "Group8 (threshold=3mil)", "Group9 (threshold=3mil)",
      #4th layer
      "Group10 (threshold=4mil)","Group11 (threshold=4mil)",
      "Group10 (threshold=4mil)","Group10 (threshold=4mil)",
      "Group11 (threshold=4mil)","Group11 (threshold=4mil)",
      "Group10 (threshold=4mil)", "Group11 (threshold=4mil)",
      "Group10 (threshold=4mil)", "Group11 (threshold=4mil)",
      "Group10 (threshold=4mil)", "unmatched12 (threshold=4mil)",
      #5th layer
      "Group12 (threshold=5mil)","Group12 (threshold=5mil)",
      "Group12 (threshold=5mil)","Group12 (threshold=5mil)",
      "Group11 (threshold=5mil)","unmatched5 (threshold=5mil)",
      "Group12 (threshold=5mil)", "Group12 (threshold=5mil)",
      "Group12 (threshold=5mil)", "Group12 (threshold=5mil)",
      "unmatched11 (threshold=5mil)", "unmatched12 (threshold=5mil)",
      #6th layer
      "Group13 (threshold=7.5mil)","Group13 (threshold=7.5mil)",
      "Group13 (threshold=7.5mil)","Group13 (threshold=7.5mil)",
      "Group13 (threshold=7.5mil)","Group13 (threshold=7.5mil)",
      "Group13 (threshold=7.5mil)", "Group13 (threshold=7.5mil)",
      "Group13 (threshold=7.5mil)", "Group13 (threshold=7.5mil)",
      "Group13 (threshold=7.5mil)", "Group13 (threshold=7.5mil)"
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
                           colourScale = my_color,
                           sinksRight = FALSE)
p

