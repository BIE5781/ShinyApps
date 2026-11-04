# Run local
if(!require(shiny)){install.packages("shiny"); library(shiny)}
runGitHub("ShinyApps", "BIE5781", subdir = "OLS")
runGitHub("ShinyApps", "BIE5781", subdir = "Bayes")

# Run from local files (set the working directory to this folder)
runApp("OLS")
runApp("Bayes")
