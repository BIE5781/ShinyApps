## Ajuste Bayesiano explícito da inclinação de uma reta que passa pelo fulcro (priors conjugadas)
## Dados: medidas do meridiano de Paris (Stigler 1999, tab. 17.1)
## Coloque este arquivo e o "Stigler_99_tab_17_1.csv" na mesma pasta
## e rode com shiny::runApp()

library(shiny)
library(dplyr)
library(stringr)
library(ggplot2)

## ---- Dados ---------------------------------------------------------------
tab1 <-
    read.csv("Stigler_99_tab_17_1.csv") |>
    mutate(
        deg = as.numeric(str_extract(Midpoint, "^\\d+")),
        min = as.numeric(str_extract(Midpoint, "(?<=°|◦)\\s*\\d+")),
        sec = as.numeric(str_extract(Midpoint, "(?<=′|'|´)\\s*\\d+")),
        decimal_deg = deg + (min / 60) + (sec / 3600),
        Radians = decimal_deg * (pi / 180)) |>
    select(-deg, -min, -sec, -decimal_deg) |>
    mutate(y = Modules/Degrees, x = (sin(Radians))^2) |>
    select(X, x, y)

## Reta que passa pelo fulcro: y = y.m + b (x - x.m); só a inclinação b é desconhecida
x.m <- mean(tab1$x)
y.m <- mean(tab1$y)
n <- nrow(tab1)
z <- tab1$x - x.m
d <- tab1$y - y.m
Sxx <- sum(z^2)
b.mv <- sum(z * d) / Sxx                              # estimativa de máxima verossimilhança
rss <- sum((d - b.mv * z)^2)
sigma.mv <- sqrt(rss / (n - 1))
sigma.lim <- c(0.1, 20)

x.lim <- range(tab1$x) + c(-0.02, 0.02)
y.lim <- range(tab1$y) + c(-60, 60)

fmt <- function(v, dig = 1)
    format(round(v, dig), big.mark = ".", decimal.mark = ",", nsmall = dig)

## Conversão para metros: 1 módulo = 2 toesas; 1 toesa = 1,949 m
m.por.modulo <- 2 * 1.949

dinvgamma <- function(x, a, b) exp(a * log(b) - lgamma(a) - (a + 1) * log(x) - b / x)
## t de locação-escala (df = Inf dá a normal)
dt_ls <- function(x, df, m, s) dt((x - m) / s, df) / s

## ---- Interface ------------------------------------------------------------
ui <- fluidPage(
    titlePanel("Inferência Bayesiana para a inclinação da reta do meridiano"),
    sidebarLayout(
        sidebarPanel(
            h4("Prior da inclinação: b ~ N(mu, σ)"),
            sliderInput("m_b", "Média de b", min = -1000, max = 1500, value = 0, step = 10),
            sliderInput("s_b", "Desvio padrão de b", min = 1, max = 2000, value = 500, step = 1),
            hr(),
            checkboxInput("var_desc", "Variância desconhecida (prior para σ²)", FALSE),
            conditionalPanel(
                "!input.var_desc",
                sliderInput("sigma", "Desvio padrão conhecido (σ)",
                            min = sigma.lim[1], max = sigma.lim[2],
                            value = round(sigma.mv, 1), step = 0.1)),
            conditionalPanel(
                "input.var_desc",
                sliderInput("s0", "Palpite para σ (módulos)",
                            min = sigma.lim[1], max = sigma.lim[2],
                            value = round(sigma.mv, 1), step = 0.1),
                sliderInput("nu0", "Confiança no palpite (ν0, obs. equivalentes)",
                            min = 0.5, max = 50, value = 1, step = 0.5),
                helpText("σ² ~ Inv-Gama(ν0/2, ν0·σ0²/2); b | σ² ~ N(m, σ²·s²/σ̂²), ",
                         "com σ̂ o desvio padrão residual do ajuste de mínimos quadrados.")),
            checkboxInput("sorteios", "Mostrar retas sorteadas da posterior", TRUE),
            hr(),
            uiOutput("resumo"),
            hr(),
            h4("Previsões do modelo"),
            uiOutput("metro")
        ),
        mainPanel(
            plotOutput("grafico", height = "400px"),
            plotOutput("grafico_dens", height = "350px")
        )
    )
)

## ---- Servidor -------------------------------------------------------------
server <- function(input, output, session) {

    ## Posterior de b e, se for o caso, de σ². Marginais de b: list(df, m, s)
    post <- reactive({
        req(input$s_b > 0)
        m0 <- input$m_b
        if (!input$var_desc) {
            req(input$sigma > 0)
            prec <- 1 / input$s_b^2 + Sxx / input$sigma^2
            list(
                pri = list(df = Inf, m = m0, s = input$s_b),
                ver = list(df = Inf, m = b.mv, s = input$sigma / sqrt(Sxx)),
                pos = list(df = Inf, m = (m0 / input$s_b^2 + sum(z * d) / input$sigma^2) / prec,
                           s = sqrt(1 / prec)),
                var = NULL)
        } else {
            req(input$s0 > 0)
            k0 <- sigma.mv^2 / input$s_b^2
            a0 <- input$nu0 / 2
            b0 <- input$nu0 * input$s0^2 / 2
            kn <- k0 + Sxx
            mn <- (k0 * m0 + sum(z * d)) / kn
            an <- a0 + n / 2
            bn <- b0 + 0.5 * (sum(d^2) + k0 * m0^2 - kn * mn^2)
            list(
                pri = list(df = 2 * a0, m = m0, s = sqrt(b0 / (a0 * k0))),
                ver = list(df = n - 1, m = b.mv, s = sqrt(rss / (n - 1) / Sxx)),
                pos = list(df = 2 * an, m = mn, s = sqrt(bn / (an * kn))),
                var = list(pri = c(a0, b0), ver = c((n - 1) / 2, rss / 2), pos = c(an, bn)))
        }
    })

    output$resumo <- renderUI({
        pos <- post()$pos
        q <- qt(0.975, pos$df) * pos$s
        tagList(
            h4("Posterior de b (média e IC 95%)"),
            p(strong(sprintf("%s [%s; %s]", fmt(pos$m), fmt(pos$m - q), fmt(pos$m + q)))),
            p("Estimativa de mínimos quadrados:", strong(fmt(b.mv)))
        )
    })

    output$metro <- renderUI({
        pos <- post()$pos
        ## distância equador-polo = 90 * (a + b/2), com a = y.m - b x.m
        k <- 90 * (0.5 - x.m)
        mu <- 90 * y.m + k * pos$m
        q <- qt(0.975, pos$df) * abs(k) * pos$s
        conv <- m.por.modulo / 1e7
        tagList(
            p("Distância do equador ao polo:", br(),
              strong(paste(fmt(mu, 0), "módulos"))),
            p("Metro estimado:", br(),
              strong(sprintf("%s m (± %s)", fmt(mu * conv, 5), fmt(q * conv, 5)))),
            p("Diferença em relação a 1 m:", br(),
              strong(paste(fmt((mu * conv - 1) * 1000, 2), "mm")))
        )
    })

    output$grafico <- renderPlot({
        pos <- post()$pos
        xs <- seq(x.lim[1], x.lim[2], length.out = 100)
        q <- qt(0.975, pos$df) * pos$s
        faixa <- data.frame(x = xs, fit = y.m + pos$m * (xs - x.m),
                            q = q * abs(xs - x.m))

        g <- ggplot()
        if (input$sorteios) {
            k <- 50
            b <- pos$m + pos$s * (if (is.finite(pos$df)) rt(k, pos$df) else rnorm(k))
            ret <- data.frame(id = rep(1:k, each = 2), x = rep(x.lim, k),
                              y = as.vector(sapply(b, \(bi) y.m + bi * (x.lim - x.m))))
            g <- g + geom_line(data = ret, aes(x, y, group = id), alpha = 0.25, color = "blue")
        }
        g +
            geom_ribbon(data = faixa, aes(x, ymin = fit - q, ymax = fit + q),
                        fill = "blue", alpha = 0.2) +
            geom_line(data = faixa, aes(x, fit), color = "blue", linewidth = 1.2) +
            geom_point(data = tab1, aes(x, y), size = 4) +
            annotate("point", x = x.m, y = y.m, shape = 4, size = 6, stroke = 2, color = "blue") +
            coord_cartesian(xlim = x.lim, ylim = y.lim, expand = FALSE) +
            xlab(expression(sen^2 ~ "(latitude)")) +
            ylab("Comprimento de um grau (km)") +
            theme_bw(base_size = 18)
    })

    ## Prior, verossimilhança e posterior
    output$grafico_dens <- renderPlot({
        p <- post()
        nomes <- c("Prior", "Verossimilhança", "Posterior")
        rng <- range(p$ver$m + c(-5, 5) * p$ver$s, p$pos$m + c(-5, 5) * p$pos$s)
        vb <- seq(rng[1], rng[2], length.out = 500)
        db <- do.call(rbind, lapply(1:3, \(i) {
            m <- list(p$pri, p$ver, p$pos)[[i]]
            data.frame(v = vb, dens = dt_ls(vb, m$df, m$m, m$s), dist = nomes[i],
                       par = "b (inclinação)")
        }))
        dens <- db

        if (!is.null(p$var)) {
            top <- 4 * sigma.mv^2
            vs <- seq(top / 500, top, length.out = 500)
            dv <- do.call(rbind, lapply(1:3, \(i) {
                ab <- p$var[[c("pri", "ver", "pos")[i]]]
                data.frame(v = vs, dens = dinvgamma(vs, ab[1], ab[2]), dist = nomes[i],
                           par = "σ² (variância)")
            }))
            dens <- rbind(db, dv)
        }
        dens$dist <- factor(dens$dist, nomes)
        ggplot(dens, aes(v, dens, color = dist, fill = dist)) +
            geom_area(alpha = 0.15, position = "identity") +
            geom_line(linewidth = 1.1) +
            facet_wrap(~par, scales = "free") +
            scale_color_manual(values = c("grey40", "darkgreen", "blue")) +
            scale_fill_manual(values = c("grey40", "darkgreen", "blue")) +
            labs(x = NULL, y = "Densidade", color = NULL, fill = NULL) +
            theme_bw(base_size = 16) + theme(legend.position = "bottom")
    })
}

shinyApp(ui, server)
