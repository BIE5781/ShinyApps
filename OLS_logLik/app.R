## Ajuste de uma reta por mínimos quadrados, "na mão"
## Dados: medidas do meridiano de Paris (Stigler 1999, tab. 17.1)
## Coloque este arquivo e o "Stigler_99_tab_17_1.csv" na mesma pasta
## e rode com shiny::runApp()

library(shiny)
library(dplyr)
library(stringr)
library(ggplot2)

## ---- Dados ---------------------------------------------------------------
## Conversão: 1 módulo = 2 toesas; 1 toesa = 1,949 m
m.por.toesa <- 1.949
toesas.por.modulo <- 2
km.por.modulo <- toesas.por.modulo * m.por.toesa / 1000

tab1 <-
    read.csv("Stigler_99_tab_17_1.csv") |>
    mutate(
        deg = as.numeric(str_extract(Midpoint, "^\\d+")),
        min = as.numeric(str_extract(Midpoint, "(?<=°|◦)\\s*\\d+")),
        sec = as.numeric(str_extract(Midpoint, "(?<=′|'|´)\\s*\\d+")),
        decimal_deg = deg + (min / 60) + (sec / 3600),
        Radians = decimal_deg * (pi / 180)) |>
    select(-deg, -min, -sec, -decimal_deg) |>
    mutate(y = Modules * km.por.modulo / Degrees,   # km por grau
           x = (sin(Radians))^2) |>
    select(X, x, y)

## Fulcro: ponto (média de x, média de y)
x.m <- mean(tab1$x)
y.m <- mean(tab1$y)

## Soma dos quadrados dos resíduos de uma reta que passa pelo fulcro
sqr <- function(b) {
    pred <- y.m + b * (tab1$x - x.m)
    sum((tab1$y - pred)^2)
}

## Limites fixos dos eixos, para o gráfico não "pular" ao mudar a inclinação
x.lim <- range(tab1$x) + c(-0.02, 0.02)
y.lim <- range(tab1$y) + c(-0.25, 0.25)

## Verossimilhança do modelo Normal para uma reta que passa pelo fulcro.
## O desvio-padrão fica fixo na sua estimativa de máxima verossimilhança
## (raiz da soma dos quadrados mínima dividida por n).
n <- nrow(tab1)
sigma.mv <- sqrt(sum(resid(lm(y ~ x, data = tab1))^2) / n)

veross <- function(b) {
    pred <- y.m + b * (tab1$x - x.m)
    prod(dnorm(tab1$y, mean = pred, sd = sigma.mv))
}

## Curvas da soma dos quadrados e da verossimilhança para os painéis opcionais
b.grid <- seq(-2, 6, by = 0.01)
curva.sqr <- data.frame(b = b.grid, sqr = sapply(b.grid, sqr))
curva.veross <- data.frame(b = b.grid, L = sapply(b.grid, veross))

fmt <- function(v, dig = 1)
    format(round(v, dig), big.mark = ".", decimal.mark = ",", nsmall = dig)

## Formata com algarismos significativos (para valores muito pequenos)
fmt.sig <- function(v, sig = 4)
    format(signif(v, sig), decimal.mark = ",", big.mark = ".", scientific = FALSE)

## ---- Interface ------------------------------------------------------------
ui <- fluidPage(
    titlePanel("Qual reta descreve melhor as medidas do meridiano?"),
    sidebarLayout(
        sidebarPanel(
            sliderInput("b", "Inclinação da reta (b, em km)",
                        min = -2, max = 6, value = 0, step = 0.02,
                        animate = animationOptions(interval = 100)),
            checkboxInput("residuos", "Mostrar os resíduos", value = TRUE),
            checkboxInput("curva", "Mostrar a soma dos quadrados para cada inclinação",
                          value = FALSE),
            checkboxInput("veross", "Mostrar a verossimilhança para cada inclinação",
                          value = FALSE),
            hr(),
            h4("Soma dos quadrados dos resíduos (km²)"),
            div(style = "font-size: 2.2em; font-weight: bold;",
                textOutput("sqr_valor")),
            br(),
            uiOutput("equacao"),
            hr(),
            h4("Previsões do modelo"),
            uiOutput("metro")
        ),
        mainPanel(
            plotOutput("grafico", height = "450px"),
            uiOutput("paineis")
        )
    )
)

## ---- Servidor -------------------------------------------------------------
server <- function(input, output, session) {

    ajuste <- reactive({
        b <- input$b
        a <- y.m - b * x.m
        dados <- mutate(tab1, pred = y.m + b * (x - x.m))
        list(dados = dados, a = a, b = b, sqr = sqr(b))
    })

    output$sqr_valor <- renderText(fmt.sig(ajuste()$sqr))

    output$metro <- renderUI({
        aj <- ajuste()
        dist.km  <- 90 * (aj$a + aj$b / 2)      # equador ao polo, em km
        dist.m   <- dist.km * 1000               # em metros
        metro    <- dist.m / 1e7                 # metro estimado
        dif      <- metro - 1
        tagList(
            p("Distância do equador ao polo norte:", br(),
              strong(paste(fmt(dist.km, 2), "km")), br(),
              strong(paste(fmt(dist.m, 0), "m"))),
            p("Metro estimado (distância / 10", tags$sup(7), "):", br(),
              strong(paste(fmt(metro, 5), "m"))),
            p("Diferença em relação a 1 m:", br(),
              strong(paste(fmt(dif, 5), "m", sprintf("(%s mm)", fmt(dif * 1000, 2)))))
        )
    })

    output$equacao <- renderUI({
        aj <- ajuste()
        withMathJax(sprintf("$$y = %s + %s \\, x$$",
                            fmt(aj$a, 3), fmt(aj$b, 2)))
    })

    output$grafico <- renderPlot({
        aj <- ajuste()
        reta <- data.frame(x = x.lim, y = y.m + aj$b * (x.lim - x.m))

        p <- ggplot(aj$dados, aes(x, y))
        if (input$residuos)
            p <- p + geom_segment(aes(xend = x, yend = pred),
                                  color = "red", linetype = "dashed", linewidth = 0.8)
        p +
            geom_line(data = reta, color = "blue", linewidth = 1.2) +
            geom_point(size = 4) +
            annotate("point", x = x.m, y = y.m, shape = 4, size = 6, stroke = 2,
                     color = "blue") +
            annotate("text", x = x.m, y = y.m, label = "fulcro",
                     vjust = -1.2, color = "blue", size = 5) +
            coord_cartesian(xlim = x.lim, ylim = y.lim, expand = FALSE) +
            xlab(expression(sen^2 ~ "(latitude)")) +
            ylab("Comprimento de um grau (km)") +
            theme_bw(base_size = 18)
    })

    ## Painéis opcionais: lado a lado se os dois estiverem ligados
    output$paineis <- renderUI({
        graficos <- list()
        if (input$curva)
            graficos <- c(graficos, list(plotOutput("grafico_sqr", height = "320px")))
        if (input$veross)
            graficos <- c(graficos, list(plotOutput("grafico_veross", height = "320px")))
        if (length(graficos) == 0) return(NULL)
        do.call(fluidRow, lapply(graficos, column, width = 12 / length(graficos)))
    })

    output$grafico_veross <- renderPlot({
        aj <- ajuste()
        L.atual  <- veross(aj$b)
        expoente <- floor(log10(L.atual))
        mantissa <- format(round(L.atual / 10^expoente, 2),
                           decimal.mark = ",", nsmall = 2)
        ggplot(curva.veross, aes(b, L)) +
            geom_line(linewidth = 1) +
            annotate("point", x = aj$b, y = L.atual, color = "red", size = 4) +
            ggtitle(bquote(L(b) == .(mantissa) %*% 10^.(expoente))) +
            labs(caption = paste0("σ fixo na estimativa de MV: ",
                                  fmt(sigma.mv * 1000, 1), " m")) +
            xlab("Inclinação (b, em km)") +
            ylab("Verossimilhança") +
            theme_bw(base_size = 16)
    })

    output$grafico_sqr <- renderPlot({
        aj <- ajuste()
        ggplot(curva.sqr, aes(b, sqr)) +
            geom_line(linewidth = 1) +
            annotate("point", x = aj$b, y = aj$sqr, color = "red", size = 4) +
            xlab("Inclinação (b, em km)") +
            ylab("Soma dos quadrados\ndos resíduos (km²)") +
            theme_bw(base_size = 16)
    })
}

shinyApp(ui, server)
