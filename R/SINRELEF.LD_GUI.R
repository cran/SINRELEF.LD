SINRELEF.LD_GUI <- function(){

  # Define UI
  ui <- fluidPage(
    titlePanel(div("SINRELEF", style = "color: #304674" )),
    sidebarLayout(
      sidebarPanel(
        tags$style(".well {background-color:#A7C7E7;}"),
        # Radio with 2 options
        radioButtons("selector", "Model:",
                     choices = c("Linear", "Graded"),
                     selected = "Linear"),

        # Input panel

        fileInput("L", "Lambdas vector (L)"),

        conditionalPanel(
          condition = "input.selector == 'Linear'",
          fileInput("PSI", "Residual variances vector (PSI)")
        ),

        conditionalPanel(
          condition = "input.selector == 'Graded'",
          fileInput("THRES", "Thresholds (THRES)")
        ),

        conditionalPanel(
          condition = "input.selector == 'Graded'",
          textInput("ncat", "Categories (ncat)")
        ),

        fileInput("doublet_list", "Doublets vector (doublet_list)"),

        fileInput("cor_doublet", "Doublets correlation (cor_doublet)"),

        textInput("N", "Sample size (N)"),

        textInput("CI", "Confidence Interval (CI)", 90),

        actionButton("compute_button", "Compute")


      ),
      mainPanel(

        #div(id = "spinner", class = "shiny-spinner", style = "display:none;"),

        withSpinner(verbatimTextOutput("output_results", placeholder = TRUE)),

        useShinyjs()

        #useWaiter(),
        #verbatimTextOutput("output_results")

        #verbatimTextOutput("output_results") %>% withSpinner()

        #hidden(div(id = 'spinner', withSpinner(verbatimTextOutput(outputId = "output_results"))))
      )
    )
  )

  # Define server
  server <- function(input, output) {

    output$output_results <- renderText({
      "Welcome to SINRELEF GUI version!\n
Please, configure the analysis using the side menu.\n
Keep in mind that graded model can be potentially very slow depending on data size.\n
If you have any doubts, please revise SINRELEF.LD documentation"
    })

    observeEvent(input$compute_button, {

      runjs("window.scrollTo(0, 0);")

      #waiter_show(html = waiting_screen, color = "black")


      if (is.null(input$N) || input$N == "") {
        showModal(modalDialog(
          title = "Warning",
          "Please, specify the sample size (N) on the appropiate field.",
          easyClose = TRUE
        ))
      }
      else {

        if (is.null(input$L) || is.null(input$doublet_list) || is.null(input$cor_doublet)){
          showModal(modalDialog(
            title = "Warning",
            "Some requiered fields are empty. Please, select a data file for each field.",
            easyClose = TRUE
          ))
        }
        else {

          # Sanity checks, change all sep to \t
          buff <- readLines(input$L$datapath)
          buff <- gsub("[;, ]", "\t",buff)
          L <- as.matrix(read.table(text = buff))

          #L <- as.matrix(read.table(input$L$datapath))
          buff <- readLines(input$doublet_list$datapath)
          buff <- gsub("[;, ]", "\t",buff)
          doublet_list <- as.matrix(read.table(text = buff))

          #doublet_list <- as.matrix(read.table(input$doublet_list$datapath))

          buff <- readLines(input$cor_doublet$datapath)
          buff <- gsub("[;, ]", "\t",buff)
          cor_doublet <- as.matrix(read.table(text = buff))

          #cor_doublet <- as.matrix(read.table(input$cor_doublet$datapath))
          N <- as.numeric(input$N)
          CI <- as.numeric(input$CI)

          if (input$selector == "Linear"){

            if (is.null(input$PSI)){
              showModal(modalDialog(
                title = "Warning",
                "PSI field is empty. Please, select a data file for PSI field.",
                easyClose = TRUE
              ))
            }
            else {

              buff <- readLines(input$PSI$datapath)
              buff <- gsub("[;, ]", "\t",buff)
              PSI <- as.matrix(read.table(text = buff))

              #PSI <- as.matrix(read.table(input$PSI$datapath))

              # Waiting bar
              output$output_results <- renderPrint({
                SINRELEF.LD(L=L, PSI=PSI, model = 'linear', doublet_list=doublet_list,
                            cor_doublet=cor_doublet, N=N, CI = CI, display = TRUE)
              })
              #waiter_hide()
            }


          }
          else {

            if (is.null(input$THRES)){
              showModal(modalDialog(
                title = "Warning",
                "THRES field is empty. Please, select a data file for THRES field.",
                easyClose = TRUE
              ))
            }

            else {
              if (is.null(input$ncat) || input$ncat == "") {
                showModal(modalDialog(
                  title = "Warning",
                  "Please, specify the number of categories (ncat) on the appropiate field.",
                  easyClose = TRUE
                ))
              }
              else {

                buff <- readLines(input$THRES$datapath)
                buff <- gsub("[;, ]", "\t",buff)
                THRES <- as.matrix(read.table(text = buff))

                #THRES <- as.matrix(read.table(input$THRES$datapath))
                ncat <- as.numeric(input$ncat)

                # Waiting bar

                output$output_results <- renderPrint({
                result <- SINRELEF.LD(L=L, model = 'graded', THRES = THRES, ncat = ncat, doublet_list=doublet_list,
                          cor_doublet=cor_doublet, N=N, CI = CI, display = TRUE)
                })
                #waiter_hide()
              }
            }

          }
        }
      }

    })
  }


  # Crea la aplicación Shiny
  shinyApp(ui, server)

}
