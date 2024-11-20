

# model selection ---------------------------------------------------------
# render html
rmarkdown::render(
  input = "modelSelection/modelSelection.Rmd",
  output_file = "modelSelection.html"
)

# extract code
knitr::purl(input = "modelSelection/modelSelection.Rmd", output = "modelSelection/modelSelection.R", documentation = 1)
