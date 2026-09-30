# Run local
if(!require(shiny)){install.packages("shiny"); library(shiny)}
runGitHub("ShinyApps", "BIE5781", subdir = "OLS")
