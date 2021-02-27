# shiny app for interactive inspection of data

library(shiny)
library(extrafont)
source("cw_functions.R")

interactive_cw <- function(ml = data) {
    theme_common <- 
        theme_bw() +
        theme(text = element_text(size = 20, family = "Arial"),
              panel.grid = element_blank(),
              legend.position = "bottom")
    
    Genotypes <- ml$info$Genotype %>% unique()
    Samples <- ml$info$Pyrat_id %>% unique() %>% sort()
    
    survivaldata <- survival_data(ml)$entries
    totalentries <- entries_data(ml)$total
    subtypes <- entries_data(ml)$subtypes
    # cw <- ml$cw
    
    ui <- fluidPage(
        
        titlePanel("Shiny Cognition Wall"),
        
        sidebarLayout(
            sidebarPanel(
                checkboxGroupInput("genotypes2show", "Genotypes:",
                                   choices = Genotypes, selected = Genotypes),
                checkboxGroupInput("samples2show", "Samples:",
                                   choices = Samples, 
                                   # inline = T,
                                   selected = Samples),
                numericInput("maxvalue", "Maximum value for X-axis:", min = 500, max = 2000, value = 1200),
                numericInput("ticks", "Ticks per:", min = 50, max = 1000, value = 200) #,
                # actionButton("applysettings","Apply!")
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Collective Evaluation",
                             fluidRow(plotOutput("survivalplot")),
                             br(),
                             fluidRow(
                                 column(6, plotOutput("totalplot")),
                                 column(6, plotOutput("subtypesplot"))
                             )
                    ),
                    tabPanel("Individual Inspection",
                             plotOutput("accuracyplot")
                    )
                )
            )
        )
    )
    
    server <- function(input, output) {
        
        output$survivalplot <- renderPlot({
            data <- survivaldata %>%
                filter(Genotype %in% input$genotypes2show) %>%
                filter(Pyrat_id %in% input$samples2show) %>% 
                group_by(Phase,Genotype) %>%
                mutate(Fraction = (1:n())/n()) %>% ungroup()
            
            discrimination <- survdiff(Surv(Entries_stat,Status) ~ Genotype,
                                       data = data %>%
                                           filter(Phase == "Discrimination"))
            reversal <- survdiff(Surv(Entries_stat,Status) ~ Genotype,
                                 data = data %>%
                                     filter(Phase == "Reversal"))
            
            surv_stats <- list(
                discrimination = list(output = discrimination,
                                      details = tibble(
                                          Phase = "Discrimination",
                                          chi_sq = round(discrimination$chisq,4),
                                          p_value = round(pchisq(discrimination$chisq, df = length(Genotypes)-1, lower.tail = F),4))
                ),
                reversal = list(output = reversal,
                                details = tibble(
                                    Phase = "Reversal",
                                    chi_sq = round(reversal$chisq,4),
                                    p_value = round(pchisq(reversal$chisq, df = length(Genotypes)-1, lower.tail = F),4))
                )
            )
            
            annot <-
                bind_rows(surv_stats$discrimination$details, surv_stats$reversal$details) %>%
                mutate(
                    Genotype = as.character(Genotypes[1]),
                    p_value = paste("p =", p_value)
                )
            
            data %>%
                ggplot(aes(Entries_stat, Fraction*100, color = Genotype)) +
                geom_step(size = 1) +
                facet_grid(.~ Phase) +
                labs(x = NULL, y = "Proportion of mice finished (%)",
                     color = NULL) +
                theme_common +
                geom_text(data = annot, color = "black", hjust = 1, vjust = -1, size = 6,
                          mapping = aes(x = input$maxvalue, y = -Inf, label = p_value)) +
                scale_x_continuous(limits = c(0, input$maxvalue), breaks = seq(0, input$maxvalue,input$ticks))
                
        })
        
        output$totalplot <- renderPlot({
            data <- totalentries %>% 
                filter(Genotype %in% input$genotypes2show) %>%
                filter(Pyrat_id %in% input$samples2show)
            
            data %>%
                ggplot(aes(Phase, tEntries, fill = Genotype)) +
                labs(x = "", y = "Total entries",
                     fill = NULL) +
                theme_common +
                stat_summary(geom = "bar", fun = mean, position = position_dodge(.9)) +
                stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0, size = 1)
        })
        
        output$subtypesplot <- renderPlot({
            data <- subtypes %>% 
                filter(Genotype %in% input$genotypes2show) %>%
                filter(Pyrat_id %in% input$samples2show)
            
            data %>%
                ggplot(aes(Entry_type, Entries, fill = Genotype)) +
                labs(x = "Entry type", y = "Entries",
                     fill = NULL) +
                theme_common +
                stat_summary(geom = "bar", fun = mean, position = position_dodge(.9)) +
                stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0, size = 1)
        })
        
        output$accuracyplot <- renderPlot({
            data <- ml$cw %>% 
                filter(Genotype %in% input$genotypes2show) %>%
                filter(Pyrat_id %in% input$samples2show)
            
            data %>%
                ggplot(aes(Entry_id, Accuracy, group = Pyrat_id)) +
                geom_line() +
                geom_hline(yintercept = .50, linetype = "dashed") +
                theme_common +
                labs(x = "Entry ID") +
                scale_y_continuous(breaks = seq(0,1,.50)) +
                geom_point(data = ml$cw %>% 
                               filter(!is.na(Reward) & Genotype %in% input$genotypes2show & Pyrat_id %in% input$samples2show),
                           aes(Entry_id, Accuracy), color = "purple", size = 0.2) +
                geom_point(data = ml$crit80 %>% 
                               filter(Criterium == "Reached" & Genotype %in% input$genotypes2show & Pyrat_id %in% input$samples2show),
                           aes(Entry_id, Accuracy), color = "red", size = 2) +
                facet_grid(Pyrat_id ~ Phase)
        }, height = 2000)
    }
    
    shinyApp(ui, server)
}

ml <- readRDS("data/adrian_fmr1_p80.RDS")
interactive_cw(ml = ml)
