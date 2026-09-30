# Aplicativos Shiny | BIE5781

Este repositório reúne aplicativos interativos em R para apoiar as atividades da disciplina **BIE5781 - Introdução à modelagem estatística na pesquisa em biologia**. A proposta é explorar ideias de modelagem estatística manipulando parâmetros e observando, em tempo real, como eles afetam ajustes, resíduos e previsões.

## Aplicativo disponível

### OLS: ajuste de uma reta por mínimos quadrados

O aplicativo em [`OLS`](OLS/) usa dados históricos das medidas do meridiano de Paris, reproduzidos da tabela 17.1 de Stigler (1999). A interface permite variar a inclinação de uma reta que passa pelo ponto médio dos dados e observar:

- como a reta se ajusta às observações;
- os resíduos de cada observação;
- a soma dos quadrados dos resíduos para cada inclinação;
- a inclinação que minimiza essa soma e as previsões associadas ao modelo.

## Executar localmente

Você precisa ter o [R](https://cran.r-project.org/) instalado. O aplicativo usa os pacotes `shiny`, `dplyr`, `stringr` e `ggplot2`.

1. Abra um terminal na pasta deste repositório. No RStudio, você também pode abrir esta pasta como projeto ou definir o diretório de trabalho para ela.
2. Instale os pacotes necessários, se ainda não estiverem instalados:

   ```r
   install.packages(c("shiny", "dplyr", "stringr", "ggplot2"))
   ```

3. Inicie o aplicativo a partir da raiz do repositório:

   ```r
   shiny::runApp("OLS")
   ```

   Ou, no terminal, execute:

   ```sh
   Rscript -e 'shiny::runApp("OLS")'
   ```

O R abrirá o aplicativo no navegador. Para encerrá-lo, interrompa a execução no terminal ou no console do RStudio (normalmente, `Esc` ou o botão de parar).

> Mantenha [`Stigler_99_tab_17_1.csv`](OLS/Stigler_99_tab_17_1.csv) na pasta `OLS`, junto de `app.R`: o aplicativo carrega esse arquivo ao iniciar.

## Executar a versão publicada no GitHub

O arquivo [`run_apps.R`](run_apps.R) chama `shiny::runGitHub()` para buscar o aplicativo no repositório GitHub `BIE5781/ShinyApps`. Esse modo requer conexão com a internet e executa a versão publicada no GitHub; para testar os arquivos da sua cópia local, use `shiny::runApp("OLS")` conforme as instruções acima.

## Estrutura

```text
.
├── run_apps.R
└── OLS/
    ├── app.R
    └── Stigler_99_tab_17_1.csv
```

Novos aplicativos da disciplina podem ser organizados em subpastas próprias, cada uma com seu arquivo `app.R` e os dados de que precisar.