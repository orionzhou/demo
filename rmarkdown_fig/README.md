# Generate PDF Report using R Markdown

* Install R Markdown as an R library
In [RStudio](https://www.rstudio.com) terminal panel type:

    install.packages("rmarkdown")

then load it:

    library(rmarkdown)

* Download `demo.Rmd`, `rmd-latex-ms.tex` (a template file) and `maize.png` (the figure(s) you want to insert into PDF);
* Open the `demo.Rmd` in RStudio;
* (Optional) Make modifications (insert additional figures, etc.);
* Select `Knit` -> `Knit to PDF` in the menu;

<img src="rmd_fig.png" width="400">

* You should then get the output [demo.pdf](demo.pdf)

Check [this tutorial](https://rmarkdown.rstudio.com/lesson-1.html) for more R Markdown tutorial.
