## Distribuição normal bivariada: superfície 3D com as marginais
## projetadas nas paredes e as curvas de nível projetadas no chão.
## Disciplina BIE5781 - IB - USP
##
## Para rodar: shiny::runApp("normal_bivariada")
## Pacotes necessários: shiny, plotly

library(shiny)
library(plotly)

## Limites dos eixos x e y: fixos, para que se veja a distribuição mudar
## de posição quando as médias mudam
lim <- c(-7, 7)
grade <- seq(lim[1], lim[2], length.out = 81)

## Densidade da normal bivariada numa grade de valores
## (matriz: linhas = valores de x2, colunas = valores de x1, como o plotly espera)
dens_bivariada <- function(x, y, mu1, mu2, s1, s2, rho) {
    z1 <- outer(y, x, function(y, x) (x - mu1) / s1)
    z2 <- outer(y, x, function(y, x) (y - mu2) / s2)
    q <- (z1^2 - 2 * rho * z1 * z2 + z2^2) / (1 - rho^2)
    exp(-q / 2) / (2 * pi * s1 * s2 * sqrt(1 - rho^2))
}

## Números com vírgula decimal, para o LaTeX (a vírgula entre chaves evita
## o espaço que o LaTeX coloca depois de vírgulas)
fmt <- function(x) sub(",", "{,}", format(round(x, 2), nsmall = 2, decimal.mark = ","),
                       fixed = TRUE)

## Valores iniciais dos parâmetros
inicio <- list(mu1 = 0, mu2 = 0, s1 = 1, s2 = 1, rho = 0)

## Animação dos controles: um passo a cada 500 ms, em ciclo
anima <- animationOptions(interval = 500, loop = TRUE)

ui <- fluidPage(
    withMathJax(),
    titlePanel("Distribuição normal bivariada"),
    sidebarLayout(
        sidebarPanel(
            width = 3,
            h4("Médias"),
            sliderInput("mu1", "\\(\\mu_1\\)", min = -2, max = 2, value = inicio$mu1,
                        step = 0.25, animate = anima),
            sliderInput("mu2", "\\(\\mu_2\\)", min = -2, max = 2, value = inicio$mu2,
                        step = 0.25, animate = anima),
            h4("Desvios-padrão"),
            sliderInput("s1", "\\(\\sigma_1\\)", min = 0.5, max = 2, value = inicio$s1,
                        step = 0.1, animate = anima),
            sliderInput("s2", "\\(\\sigma_2\\)", min = 0.5, max = 2, value = inicio$s2,
                        step = 0.1, animate = anima),
            h4("Correlação"),
            sliderInput("rho", "\\(\\rho\\)", min = -0.95, max = 0.95, value = inicio$rho,
                        step = 0.05, animate = anima),
            hr(),
            checkboxInput("z.fixo", "Congelar a escala do eixo z", value = FALSE),
            checkboxInput("marginais", "Mostrar as marginais nas paredes", value = TRUE),
            actionButton("reinicia", "Voltar aos valores iniciais"),
            br(), br(),
            helpText("Arraste o gráfico para girá-lo. O botão ▶ de cada",
                     "controle anima o parâmetro.")
        ),
        mainPanel(
            width = 9,
            plotlyOutput("grafico", height = "620px"),
            fluidRow(
                column(6, uiOutput("modelo")),
                column(6, uiOutput("covariancia"))
            )
        )
    )
)

server <- function(input, output, session) {

    observeEvent(input$reinicia, {
        for (p in names(inicio)) updateSliderInput(session, p, value = inicio[[p]])
    })

    ## Densidade conjunta e marginais com os valores atuais dos parâmetros
    dens <- reactive({
        z <- dens_bivariada(grade, grade, input$mu1, input$mu2,
                            input$s1, input$s2, input$rho)
        list(z = z,
             marg1 = dnorm(grade, input$mu1, input$s1),
             marg2 = dnorm(grade, input$mu2, input$s2))
    })

    ## Limite superior do eixo z: automático, ou congelado no valor do
    ## momento em que a opção foi marcada
    z.auto <- reactive({
        d <- dens()
        max(c(d$z, if (input$marginais) c(d$marg1, d$marg2))) * 1.05
    })
    z.congelado <- reactiveVal(NULL)
    observeEvent(input$z.fixo, {
        z.congelado(if (input$z.fixo) isolate(z.auto()) else NULL)
    })
    z.max <- reactive(if (is.null(z.congelado())) z.auto() else z.congelado())

    ## O gráfico é desenhado uma vez, com os valores iniciais. Depois, só
    ## os dados são atualizados (observer abaixo), o que mantém o ângulo
    ## de visão escolhido com o mouse e deixa a animação mais suave
    output$grafico <- renderPlotly({
        d <- isolate(dens())
        plot_ly() |>
            add_surface(x = grade, y = grade, z = d$z,
                        colorscale = "Viridis", showscale = FALSE,
                        name = "Densidade conjunta",
                        hovertemplate = paste0("x<sub>1</sub> = %{x:.2f}<br>",
                                               "x<sub>2</sub> = %{y:.2f}<br>",
                                               "densidade = %{z:.3f}<extra></extra>"),
                        ## Curvas de nível projetadas no chão
                        contours = list(z = list(show = TRUE, usecolormap = TRUE,
                                                 project = list(z = TRUE),
                                                 start = 0, end = max(d$z),
                                                 size = max(d$z) / 8))) |>
            ## Marginal de X1 na parede do fundo (x2 no limite superior)
            add_trace(type = "scatter3d", mode = "lines",
                      x = grade, y = rep(lim[2], length(grade)), z = d$marg1,
                      line = list(color = "#dc322f", width = 6),
                      name = "Marginal de X<sub>1</sub>",
                      hovertemplate = "x<sub>1</sub> = %{x:.2f}<br>densidade = %{z:.3f}<extra></extra>") |>
            ## Marginal de X2 na parede lateral (x1 no limite inferior)
            add_trace(type = "scatter3d", mode = "lines",
                      x = rep(lim[1], length(grade)), y = grade, z = d$marg2,
                      line = list(color = "#268bd2", width = 6),
                      name = "Marginal de X<sub>2</sub>",
                      hovertemplate = "x<sub>2</sub> = %{y:.2f}<br>densidade = %{z:.3f}<extra></extra>") |>
            layout(
                legend = list(x = 0.02, y = 0.98),
                margin = list(l = 0, r = 0, t = 10, b = 0),
                scene = list(
                    xaxis = list(title = "X<sub>1</sub>", range = lim),
                    yaxis = list(title = "X<sub>2</sub>", range = lim),
                    zaxis = list(title = "Densidade", range = c(0, isolate(z.max()))),
                    aspectmode = "manual", aspectratio = list(x = 1, y = 1, z = 0.7),
                    camera = list(eye = list(x = 1.6, y = -1.6, z = 1.1))
                )
            ) |>
            config(displaylogo = FALSE)
    })

    ## Atualiza os dados do gráfico quando os parâmetros mudam
    proxy <- plotlyProxy("grafico", session)
    observe({
        d <- dens()
        mostra <- input$marginais
        proxy |>
            ## traço 0: superfície; traços 1 e 2: marginais
            plotlyProxyInvoke("restyle", list(z = list(d$z),
                                              contours.z.end = max(d$z),
                                              contours.z.size = max(d$z) / 8), list(0)) |>
            plotlyProxyInvoke("restyle", list(z = list(d$marg1), visible = mostra), list(1)) |>
            plotlyProxyInvoke("restyle", list(z = list(d$marg2), visible = mostra), list(2)) |>
            plotlyProxyInvoke("relayout", list(scene.zaxis.range = c(0, z.max()),
                                               showlegend = mostra))
    })

    ## Especificação da distribuição com os valores atuais
    output$modelo <- renderUI({
        withMathJax(sprintf(
            "$$\\begin{pmatrix} X_1 \\\\ X_2 \\end{pmatrix} \\sim
              \\text{MVNormal}\\left( \\begin{pmatrix} %s \\\\ %s \\end{pmatrix}, \\, \\mathbf{S} \\right)$$",
            fmt(input$mu1), fmt(input$mu2)))
    })

    ## Matriz de covariâncias, por extenso e com os valores atuais
    output$covariancia <- renderUI({
        cov12 <- input$rho * input$s1 * input$s2
        withMathJax(sprintf(
            "$$\\mathbf{S} = \\begin{pmatrix} \\sigma_1^2 & \\rho \\, \\sigma_1 \\sigma_2 \\\\
              \\rho \\, \\sigma_1 \\sigma_2 & \\sigma_2^2 \\end{pmatrix} =
              \\begin{pmatrix} %s & %s \\\\ %s & %s \\end{pmatrix}$$",
            fmt(input$s1^2), fmt(cov12), fmt(cov12), fmt(input$s2^2)))
    })
}

shinyApp(ui, server)
