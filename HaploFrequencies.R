#!/usr/bin/env Rscript

library(readxl)
library(writexl)
library(dplyr)
library(plyr)
library(tidyr)

library(htmltools)
library(maps)
library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(leaflet)
library(ggplot2)
library(leafpop)
library(leaflet.extras2)
library(tidygeocoder)
library(sf)
library(lattice)
library(reshape2)

library(RColorBrewer)
library(paletteer)

library(rcarbon)
library(stringr)


# HaploFreqencies.R
# 
# Description: Plots haplogroup frequencies for mtDNA or YDNA and creates a frequency table.
# 
# Procedure:
#     1. creates haplogroup frequency table by grouping samples from the AADR dataset
#     2. plots haplogroup frequency table in leaflet and creates pie charts
# 
# Input: AADR dataset, Ancient samples, mtDNA tree, YDNA tree
# Output: haplogroup frequency table 
# 
# Usage: ./HaploFrequencies.R mtDNA 500 Locality 2
# 
# Version: 1.00
# Date: 23.03.2024
# Author: Jakob Materna


find_basal_haplogroup <- function(row, column, df, basal_point) {
  # finds basal haplogroup from existing haplogroup
  haplogroup <- row[column]
  if (substring(haplogroup, 1, 3) %in% c("n/a", "..")) {
    return(NA)
  } 
  else {
    repeat {
      if (setequal(nchar(haplogroup), basal_point)) {
        break
      }
      if (is.na(df[df[1] == haplogroup][2])) {
        return(NA)
      }
      haplogroup <- df[df[1] == haplogroup][2]
    }
    return(haplogroup)
  }
}

create_age_groups <- function(row, step_size) {
  # creates age group based on input and groups samples in these groups
  age <- as.numeric(row[9])
  age_bot_BP <- round_any(age, step_size, ceiling)
  age_top_BP <- age_bot_BP - step_size
  if (age_bot_BP == (2000-step_size)) {
    age_bot <- as.character(2000-step_size)
    age_top <- "1950 AD"
  }
  else if (1950 - age > 0) {
    age_bot <-as.character(2000 - age_bot_BP)
    age_top <- paste(as.character(2000 - age_bot_BP + step_size), " AD", sep="")
  }
  else {
    age_bot <- as.character(abs(2000 - age_bot_BP))
    age_top <- paste(as.character(abs(2000 - age_bot_BP + step_size)), " BC", sep="")
  }
  name <- paste(row[15], " ", age_bot, "-", age_top, sep="")
  age_group_BP <- paste(age_bot_BP, "-", age_top_BP, sep="")
  age_group_ADBC <-  paste(age_bot, "-", age_top, sep="")
  return(c(name, age_group_BP, age_group_ADBC))
}
  
calculate_haplogroups <- function(AADR, dna, mode) {
  # calculates total number of occurences of YDNA or mrDNA haplogroups
  if (dna == "mtDNA") {
    keep <- "mtDNA Basal haplogroup"
    remove <- "YDNA Basal haplogroup"
  }
  else {
    keep <- "YDNA Basal haplogroup"
    remove <- "mtDNA Basal haplogroup"
  }
  AADR <- as.data.frame(AADR) %>% 
    drop_na(keep) %>% 
    select(!remove) %>% 
    pivot_wider(names_from=keep, 
                values_from=keep, 
                values_fn=length, 
                values_fill=0) %>% 
    mutate(Total=select(., -(1:7)) %>% rowSums(na.rm=TRUE))
  # if "country" is specified the data is further grouped
  if (mode == "Country") {
    AGG <- AADR
    AGG$`Locality` <- "various"
    AGG$`Lat.` <- 0
    AGG$`Long.` <- 0
    cols <- colnames(AGG)
    columns_to_sum <- colnames(AGG)[-(1:7)]
    AGG <- aggregate(AGG[, columns_to_sum], by=list(AGG$`Ancient pop name`,
                                                    AGG$`Locality`,                                                     AADR$`Political Entity`,
                                                    AGG$`Age group in BP format`,
                                                    AGG$`Age group in AD/BC format`,
                                                    AGG$`Lat.`,
                                                    AGG$`Long.`), FUN=sum)
    colnames(AGG) <- cols
    AGG$`Lat.` <- apply(AGG, 1, function(row) {subset(AADR, `Ancient pop name` == row[["Ancient pop name"]])[[1,"Lat."]]})
    AGG$`Long.` <- apply(AGG, 1, function(row) {subset(AADR, `Ancient pop name` == row[["Ancient pop name"]])[[1,"Long."]]})
    AADR <- AGG
  }
  return(AADR)
}

create_output_table <- function(step_size, dna, mode, basal_point) {
  # cleans and groups the data 
  # opens data files
  AADR <- read_excel("data/AADR.xlsx")
  Ancient <- read.csv("data/Ancient_samples.txt", sep="\t", header=FALSE, row.names=1, col.names=c("id", "name"))
  
  MtDNA <- read.csv("data/mt_tree.txt", sep="\t", quote = "")
  YDNA <- read.csv("data/y_tree.txt", sep="\t", quote = "")
  
  # reducing the dataframe to only contain relevant entries
  AADR <- AADR[AADR$`Genetic ID` %in% Ancient$name,]
  
  # finding basal haplogroup
  AADR$`mtDNA Basal haplogroup` <- apply(AADR, 1, find_basal_haplogroup, column="mtDNA haplogroup if >2x or published", df=MtDNA, basal_point=basal_point)
  AADR$`YDNA Basal haplogroup` <- apply(AADR, 1, find_basal_haplogroup, column="Y haplogroup (manual curation in ISOGG format)", df=YDNA, basal_point=basal_point)
  
  # grouping samples into populations
  age_group <- apply(AADR, 1, create_age_groups, step_size=step_size)
  AADR$`Ancient pop name` <- age_group[1,]
  AADR$`Age group in BP format` <- age_group[2,]
  AADR$`Age group in AD/BC format` <- age_group[3,]
  
  AADR <- select(AADR,
                 "Ancient pop name",
                 "Locality",
                 "Political Entity",
                 "Age group in BP format",
                 "Age group in AD/BC format",
                 "Lat.",
                 "Long.",
                 "mtDNA Basal haplogroup",
                 "YDNA Basal haplogroup")
  
  # calculating haplogroups
  result <- calculate_haplogroups(AADR, dna=dna, mode=mode)
  
  # removing unusable data
  result$`Lat.` <- as.numeric(result$`Lat.`)
  result$`Long.` <- as.numeric(result$`Long.`)
  result <- result %>% drop_na(`Lat.`, `Long.`)

  # writes result table 
  write_xlsx(result, "output.xlsx")
}



plot_data_column <- function(data) {
  # plots popup figures
  header <- data$`Ancient pop name`
  data <- t(as.data.frame(data[6:(length(data)-2)]))
  data <- cbind(rownames(data), data.frame(data, row.names=NULL))
  colnames(data) <- c("variable", "value")
  data$color <- colors_haplo
  data <- filter(data, value != 0)
  # creates pie charts 
  p <- ggplot(data, aes(x="", y=value, fill=variable)) +
    geom_bar(stat="identity", width=1) +
    geom_bar(data = subset(data, nrow(data) != 1),
             aes(x="", y=value, fill=variable),
             stat="identity", width=1, color="#FFFFFF") +
    scale_fill_manual(values=data$color) +
    coord_polar("y", start=0) +
    labs(x="", y="", title=header) +
    theme_void() +
    theme(plot.title=element_text(hjust=0.5, size=11),
          legend.position="bottom",
          legend.title=element_blank())
}


update_markers <- function(range) {
  # Hides/Shows data points when changing double slider.
  proxy_leaf <- leafletProxy("map")
  for (i in (1:nrow(points))) {
    age_limits <- as.numeric(str_split(points[[i, "Age group in BP format"]], pattern="-", simplify=TRUE))
    age_group <- points[[i, "Age group in AD/BC format"]]
    if (range[1] <= age_limits[1] & range[2] >= age_limits[2]) {
      proxy_leaf <- proxy_leaf %>% showGroup(age_group)
    }
    else {
      proxy_leaf <- proxy_leaf %>% hideGroup(age_group)
    }
  }
  return(proxy_leaf)
}

args <- commandArgs(trailingOnly=TRUE)
# args <- c("mtDNA", 500, "Country", 1)

if (length(args) == 4) {
  dna <- args[1]                         # mtDNA/YDNA
  step_size <- as.numeric(args[2])       # 100/500/1000
  mode <- args[3]                        # Locality/Country
  basal_point <- as.numeric(args[4])     # 1/2
  
  create_output_table(step_size, dna, mode, basal_point)
  
  data <- read_excel("output.xlsx")
} else if (length(args) == 1) {
  data <- read_excel(args[1])
}

# defining maps
maps <- data.frame(name=c("Stadia.Outdoors", "USGS.USImageryTopo", "Esri.WorldImagery"),
                   group=c("Cartographic", "Geopolitical", "Satellite"))

points <- st_as_sf(data, coords=c("Long.", "Lat."), crs=4326)
points <- points[order(as.numeric(str_split(points$`Age group in BP format`, pattern="-", simplify=TRUE)[,1]), decreasing = TRUE), ]
age_groups <- unique(points$`Age group in AD/BC format`)
age_groups_BP <- unique(points$`Age group in BP format`)

split_data <- split(points, points$`Age group in AD/BC format`) 
max_age <- as.numeric(str_split(age_groups_BP[1], pattern="-", simplify=TRUE)[,1])

# loading colors
colors_haplo <- paletteer_c("ggthemes::Temperature Diverging", ncol(data)-8)
colors_age <- data.frame(age_groups, color=paletteer_c("harrypotter::gryffindor", length(age_groups)))

# creating R shiny dashboard
ui <- dashboardPage(
  # creating user interface
  dashboardHeader(title="HaploFrequencies"),
  dashboardSidebar(disable=TRUE),
  dashboardBody(
    tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
    leafletOutput("map") %>% withSpinner(type=3, color.background="#ecf0f5", color="#3C8DBC"),
    absolutePanel(top=60, left=70,
                  sliderInput("range", "", min=0, max=max_age, step=step_size, value=range(c(0, max_age)),
                              ticks=T, timeFormat=T, post=" BP", width=400))
  )
)

server <- function(input, output, session) {
  # createing leaflet map
  leaf <- leaflet() %>%
    setView(lng=20, lat=30, zoom=3) %>%
    addLayersControl(baseGroups=maps["group"][,1],
                     options=layersControlOptions(collapsed=TRUE))
  # Ad different maps
  for (m in (1:nrow(maps))) {
    leaf <- leaf %>% addProviderTiles(maps[m,"name"], group=maps[m,"group"])
  }
  # Ad markers and popups with figures
  for (group in age_groups) {
    group_data <- filter(points, `Age group in AD/BC format` == group)
    graphs <- apply(group_data, 1, plot_data_column)
    color <- colors_age[colors_age$age_groups == group,][["color"]]
    leaf <- leaf %>% addCircleMarkers(data = group_data,
                                      radius=8,
                                      weight=1,
                                      opacity=1,
                                      color="black",
                                      fillColor=color,
                                      stroke = TRUE,
                                      fillOpacity = 0.8,
                                      group=group) %>%
      addPopupGraphs(graphs, group=group, width=200, height=170)
  }
  output$map <- renderLeaflet({leaf})
  observe(observeEvent(input$range, update_markers(input$range)))
}

# launching the shiny app
runApp(shinyApp(ui, server), launch.browser=TRUE)

