#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(googlesheets4)
library(ggplot2)
library(tidyverse)

# function to fix some inputs
fix_input_gs4 <- function(x) {
  unlist(
    lapply(x, function(x) {
      if (is.null(x[[1]])) {
        x[[1]] <- NA
      } else {
        x[[1]]
      }
    })
  )
}

link <- "https://docs.google.com/spreadsheets/d/1NM_22_EHo5-9a_lSGZNDSCKmuz6nS8jJNxzupJdMUns/edit?gid=0#gid=0"
googlesheets4::gs4_auth(cache = ".secrets")

# read in the growth curve sheet and the culture log, subset the culture log to only flasks in the growth curve sheet
df <- googlesheets4::read_sheet(link, sheet = "LightExp_GrowthCurves", na = "NA")

# clear residual start points, we will add new ones
df <- df |> subset(Hours != 0)

df2 <- googlesheets4::read_sheet(link, sheet = "Culture Log", na = "NA") |>
  subset(`Culture_ID#` %in% df$`Culture_ID#` & !is.na(`Calculated Start Cell density (cells/ml)`))

# create dummy rows for the start points
df <- df[NA,][1:nrow(df2),] |>
  dplyr::mutate(`Culture_ID#` = df2$`Culture_ID#`,
         `Measurement Date/Time` = df2$`Start_Date/Time`,
         `Start_Date/Time` = df2$`Start_Date/Time`,
         Hours = 0,
         CO2_PPM = df2$CO2_PPM,
         `Growth Temp. (C)` = df2$`Growth Temp. (C)`,
         `Light Intensity` = df2$`Light Intensity`,
         `Cell Density AquaPen-Blue (Cells/mL)` = df2$`Calculated Start Cell density (cells/ml)`,
         `Cell Density AquaPen-Red (Cells/mL)` = df2$`Calculated Start Cell density (cells/ml)`,
         `Cell Density Turner (Blue) (Cells/mL)` = df2$`Calculated Start Cell density (cells/ml)`,
         `Cell Density Turner (Red) (Cells/mL)` = df2$`Calculated Start Cell density (cells/ml)`
         ) |>
  rbind(df)

df$Strain <- stringr::str_split(df$`Culture_ID#`,pattern = "_") |>
  lapply(function(x){paste(x[1],x[2],sep = "_")}) |>
  fix_input_gs4()

df$`Nominal Hours` <- as.numeric(unlist(df$Hours))
df$Strain <- fix_input_gs4(df$Strain)
df$`Culture ID` <- df$`Culture_ID#`

y_choices <- c("Cell Density AquaPen-Blue (Cells/mL)",
               'Cell Density AquaPen-Red (Cells/mL)',
               'Cell Density Turner (Blue) (Cells/mL)',
               'Cell Density Turner (Red) (Cells/mL)'
                #"Fluor. Aquapen-Blue Mean",
                #"Fluor. Aquapen-Red Mean",
                #"Turner (Blue) Mean",
                #"Turner (Red) Mean"
                )

df$`Temp/CO2` <- paste0(df$`Growth Temp. (C)`,"C, ",df$CO2_PPM,"ppm")

#df$group <- paste(df$`Start_Date/Time`,df$`Temp/CO2`,df$`Light Intensity`,sep = "_")
df$group <- df$`Culture_ID#`

color_choices <- c("Growth Temp. (C)","Temp/CO2","Light Intensity","Start_Date/Time")

cutoffs <- c("Syn" = 3.3e6, "Ehux" = 5.16e5)

# Define UI for application that draws a histogram
ui <- shiny::fluidPage(

  # Application title
  shiny::titlePanel("Growth Data"),

  # Sidebar with a slider input for number of bins
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::selectInput("sp",
        "Species",
        choices = df$Strain,
        selected = "Syn_7002",
        multiple = T
      ),
      shiny::selectInput("incu",
                  "Temp/CO2",
                  choices = df$`Temp/CO2`,
                  selected = c("25.3C, 1135ppm","22C, 420ppm" ),
                  multiple = T
      ),
      shiny::selectInput("light",
                  "Light Condition",
                  choices = df$`Light Intensity`,
                  selected = c(10,100),
                  multiple = T
      ),
      shiny::selectInput("y",
        "Measurement (y-axis)",
        choices = y_choices,
        selected = 'Cell Density Turner (Blue) (Cells/mL)',
        multiple = F
      ),
      shiny::selectInput("flask",
                  "Highlight a Flask",
                  choices = df$`Culture_ID#`,
                  multiple = T
      ),
      shiny::selectInput("color",
                  "Color Option",
                  choices = color_choices,
                  multiple = F
      ),
      shiny::selectInput("cutoff",
                  "Cutoff",
                  choices = c("Syn","Ehux"),
                  multiple = F
      ),
      shiny::checkboxInput("logscale",
                   "Use log10 y-axis",
                   value = T
      ),
      shiny::checkboxInput("old",
                    "Exclude curves that don't have an inital innoculum calculation.",
                    value = T
      )
    ),

    # Show a plot of the generated distribution
    shiny::mainPanel(
      shiny::plotOutput("plot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$plot <- shiny::renderPlot({
    
    # cut out checkbox
    if (input$old) {
      # just get curves that have an inital value
      background <- df |> subset(`Culture_ID#` %in% df[which(df$Hours == 0),]$`Culture_ID#`)
    } else {
      background <- df 
    }
    
    #checkbox for logscale
    if(input$logscale){trans = "log10"} else {trans = 'identity'}
    
    #subset df for background points
    background <- background |>
      subset(Strain %in% input$sp & `Temp/CO2` %in% input$incu & `Light Intensity` %in% input$light) |>
      dplyr::select(`Nominal Hours`,
                    `Culture ID`,
                    group,
                    y = input$y,
                    color = input$color)
    
    # plotting ggplot object
    plot <- background |>
      ggplot2::ggplot(ggplot2::aes(`Nominal Hours`,
                 y = y,
                 group = group,
                 color = as.character(color)
      )) +
      #stat_smooth(aes(group = group),span = 0.6,se = F) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_hline(yintercept = as.numeric(cutoffs[input$cutoff]),linetype = 2) +
      ggplot2::labs(y = input$y, color = input$color) +
      ggplot2::scale_y_continuous(transform = trans) + 
      ggplot2::scale_x_continuous(n.breaks = 10)
    
    if (!is.null(input$flask)) {
      highlight <- background |> subset(`Culture ID` %in% input$flask)
      plot <- plot + ggplot2::geom_point(ggplot2::aes(`Nominal Hours`,
                                    y = y,
                                    group = as.character(color)
      ),data = highlight,color = "black") +
        stat_smoothstat_smooth(aes(group = group),span = 0.6,se = F,data = highlight,color = "black") 
      ggplot2::geom_line(aes(`Nominal Hours`,
                     y = y,
                     group = as.character(color)),
                 data = highlight,color = "black")
    }
    
    plot
  }, res = 100)
}

# Run the application
shiny::shinyApp(ui = ui, server = server)
