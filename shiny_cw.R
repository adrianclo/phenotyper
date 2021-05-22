# shiny app for interactive inspection of cw data

library(shiny)
library(extrafont)
library(ggrepel)
# library(plotly)

interactive_cw <- function(ml = data, Genotypes = NULL) {
    theme_common <- 
        theme_bw() +
        theme(text = element_text(size = 20, family = "Arial"),
              panel.grid = element_blank(),
              legend.position = "bottom")
    
    if(is.null(Genotypes)) { 
        Genotypes <- ml$info$Genotype %>% unique() %>% factor() 
    }
    
    Samples <- ml$info$Pyrat_id %>% unique() %>% sort()
    
    # this dataset is used for two graphs
    # survivaldata <- survival_data(ml, factor_levels = Genotypes)$entries
    
    totalentries <- entries_data(ml, factor_levels = Genotypes)$total
    subtypes <- entries_data(ml, factor_levels = Genotypes)$subtypes
    
    ui <- fluidPage(
        
        # https://stackoverflow.com/questions/29738975/how-to-align-a-group-of-checkboxgroupinput-in-r-shiny
        tags$head(
            tags$style(
                HTML(".checkbox-inline { 
                    margin-left: 0px;
                    margin-right: 10px;
                    }
                      .checkbox-inline+.checkbox-inline {
                    margin-left: 0px;
                    margin-right: 10px;
                    }"
                )
            ) 
        ),
        
        titlePanel("Shiny Cognition Wall"),
        
        sidebarLayout(
            sidebarPanel(
                checkboxGroupInput("genotypes2show", "Genotypes:",
                                   choices = levels(Genotypes), 
                                   selected = Genotypes),
                checkboxGroupInput("samples2show", "Samples:",
                                   choices = Samples, 
                                   inline = T,
                                   selected = Samples),
                sliderInput("thresh", "Task criterium:", 
                            min = 40, max = 100, step = 5,
                            value = 80),
                numericInput("maxvalue", "Maximum value for X-axis:", 
                             min = 500, max = 2000, value = 2000),
                numericInput("ticks", "Ticks per:", 
                             min = 50, max = 1000, value = 200),
                checkboxInput("minmax", "Show min-max", value = T),
                submitButton("Apply!")
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Collective evaluation",
                             fluidRow(h3("Group sizes"),
                                      tableOutput("sample_n")),
                             fluidRow(h3("Survival data"),
                                      plotOutput("survivalplot")),
                             br(),
                             fluidRow(plotOutput("survivalbar")),
                             br(),
                             fluidRow(h3("Entry data"),
                                      column(6, plotOutput("totalplot")),
                                      column(6, plotOutput("subtypesplot")))
                    ),
                    tabPanel("Sample accuracy evaluation",
                             plotOutput("accuracyplot"),
                    ),
                    tabPanel("Sample information",
                             tableOutput("sample_info"),
                    )
                )
            )
        )
    )
    
    server <- function(input, output) {
        # reactive survival data ----
        survivaldata2 <- reactive({
            ml %>%
                new_threshold(value = input$thresh / 100) %>%
                survival_data() %>%
                pluck("entries")
        })
        
        # n ----
        output$sample_n <- renderTable({
            data <- ml$info %>%
                filter(Genotype %in% input$genotypes2show) %>%
                filter(Pyrat_id %in% input$samples2show) %>% 
                count(Genotype)
            
            genotypes <- data$Genotype %>% unique() %>% as.character()
            genotypes <- levels(Genotypes)[which(levels(Genotypes) %in% genotypes)]
            data <- data %>% 
                mutate(Genotype = factor(Genotype, levels = genotypes)) %>% 
                arrange(Genotype)
            data
        })
        
        # survival plot ----
        output$survivalplot <- renderPlot({
            data <- 
                # survivaldata %>%
                survivaldata2() %>% 
                filter(Genotype %in% input$genotypes2show) %>%
                filter(Pyrat_id %in% input$samples2show) %>%
                group_by(Phase,Genotype) %>%
                mutate(Fraction = (1:n())/n()) %>% ungroup()
            
            genotypes <- data$Genotype %>% unique() %>% as.character()
            genotypes <- levels(Genotypes)[which(levels(Genotypes) %in% genotypes)]
            data <- data %>% 
                mutate(Genotype = factor(Genotype, levels = genotypes))
            
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
                                          chi_sq = round(discrimination$chisq, 4),
                                          p_value = round(pchisq(discrimination$chisq, 
                                                                 df = length(genotypes)-1, 
                                                                 lower.tail = F), 4))
                ),
                reversal = list(output = reversal,
                                details = tibble(
                                    Phase = "Reversal",
                                    chi_sq = round(reversal$chisq, 4),
                                    p_value = round(pchisq(reversal$chisq, 
                                                           df = length(genotypes)-1, 
                                                           lower.tail = F), 4))
                )
            )
            
            annot <-
                bind_rows(surv_stats$discrimination$details, 
                          surv_stats$reversal$details) %>%
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
                scale_x_continuous(limits = c(0, input$maxvalue), 
                                   breaks = seq(0, input$maxvalue,input$ticks))
            
        })
        
        # survival bar ----
        output$survivalbar <- renderPlot({
            data <- 
                # survivaldata %>%
                survivaldata2() %>% 
                filter(Genotype %in% input$genotypes2show) %>%
                filter(Pyrat_id %in% input$samples2show) %>%
                group_by(Phase,Genotype) %>%
                mutate(Fraction = (1:n())/n()) %>% ungroup()
            
            # test line
            # data <- filter(survivaldata, !Genotype %in% c("WT_CONTROL (previous)", "KO_CONTROL (previous)"))
            
            genotypes <- data$Genotype %>% unique() %>% as.character()
            genotypes <- levels(Genotypes)[which(levels(Genotypes) %in% genotypes)]
            data <- data %>% 
                mutate(Genotype = factor(Genotype, levels = genotypes))
            
            minmaxdata <- data %>% 
                group_by(Genotype, Phase) %>% 
                filter(Entries_stat == min(Entries_stat) | 
                           Entries_stat == max(Entries_stat)) %>% 
                ungroup()
            
            g1 <- data %>% 
                ggplot(aes(Genotype, Entries_stat, fill = Genotype)) +
                theme_common +
                labs(y = "Entries to 80% criterium", x = NULL) +
                theme(legend.position = "none", 
                      axis.text.x = element_text(angle = 45, 
                                                 hjust = 1)) +
                facet_grid(.~Phase) +
                stat_summary(geom = "bar", fun = mean, position = position_dodge(.9)) +
                stat_summary(geom = "errorbar", fun.data = mean_se, color = "black", 
                             position = position_dodge(.9), width = 0, size = 1) +
                geom_point(position = position_dodge(.9))
            
            if(input$minmax) { 
                g1 +
                    ggrepel::geom_text_repel(data = minmaxdata,
                                             label = minmaxdata$Pyrat_id, 
                                             position = position_dodge(.9))
            } else { g1 }
        })
        
        # total entries ----
        output$totalplot <- renderPlot({
            data <- totalentries %>% 
                filter(Genotype %in% input$genotypes2show) %>%
                filter(Pyrat_id %in% input$samples2show)
            
            genotypes <- data$Genotype %>% unique() %>% as.character()
            genotypes <- levels(Genotypes)[which(levels(Genotypes) %in% genotypes)]
            data <- data %>% 
                mutate(Genotype = factor(Genotype, levels = genotypes))
            
            minmaxdata <- data %>% 
                group_by(Genotype, Phase) %>% 
                filter(tEntries == min(tEntries) | tEntries == max(tEntries)) %>% 
                ungroup()
            
            g2 <- data %>%
                ggplot(aes(Phase, tEntries, fill = Genotype)) +
                labs(x = "", y = "Total entries",
                     fill = NULL) +
                theme_common +
                theme(legend.position = "none") +
                stat_summary(geom = "bar", fun = mean, position = position_dodge(.9)) +
                stat_summary(geom = "errorbar", fun.data = mean_se, 
                             position = position_dodge(.9), width = 0, size = 1) +
                geom_point(position = position_dodge(.9)) 
            
            if(input$minmax) { 
                g2 +
                    ggrepel::geom_text_repel(data = minmaxdata,
                                             label = minmaxdata$Pyrat_id, 
                                             position = position_dodge(.9))
            } else { g2 }
        })
        
        # subtype entries ----
        output$subtypesplot <- renderPlot({
            data <- subtypes %>% 
                filter(Genotype %in% input$genotypes2show) %>%
                filter(Pyrat_id %in% input$samples2show)
            
            genotypes <- data$Genotype %>% unique() %>% as.character()
            genotypes <- levels(Genotypes)[which(levels(Genotypes) %in% genotypes)]
            data <- data %>% 
                mutate(Genotype = factor(Genotype, levels = genotypes))
            
            minmaxdata <- data %>% 
                group_by(Genotype, Entry_type) %>% 
                filter(Entries == min(Entries) | Entries == max(Entries)) %>% 
                ungroup()
            
            g3 <- data %>%
                ggplot(aes(Entry_type, Entries, fill = Genotype)) +
                labs(x = "Entry type", y = "Entries",
                     fill = NULL) +
                theme_common +
                theme(legend.position = "none") +
                stat_summary(geom = "bar", fun = mean, position = position_dodge(.9)) +
                stat_summary(geom = "errorbar", fun.data = mean_se, 
                             position = position_dodge(.9), width = 0, size = 1) +
                geom_point(position = position_dodge(.9))
            
            if(input$minmax) { 
                g3 +
                    ggrepel::geom_text_repel(data = minmaxdata,
                                             label = minmaxdata$Pyrat_id, 
                                             position = position_dodge(.9))
            } else { g3 }
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
                               filter(!is.na(Reward) & 
                                          Genotype %in% input$genotypes2show & 
                                          Pyrat_id %in% input$samples2show),
                           aes(Entry_id, Accuracy), color = "purple", size = .2) +
                geom_point(data = ml$crit80 %>% 
                               filter(Criterium == "Reached" & 
                                          Genotype %in% input$genotypes2show & 
                                          Pyrat_id %in% input$samples2show),
                           aes(Entry_id, Accuracy), color = "red", size = 2) +
                facet_grid(Pyrat_id ~ Phase)
        }, height = 2000)
        
        # sample info ----
        output$sample_info <- renderTable({
            data <- ml$info %>%
                select(-QC) %>% 
                mutate(Pyrat_id = as.character(Pyrat_id),
                       Genotype = factor(Genotype, levels = levels(Genotypes))) %>% 
                arrange(Genotype)
            data
        })
    }
    
    shinyApp(ui, server)
}

# interactive_cw(ml = ml)