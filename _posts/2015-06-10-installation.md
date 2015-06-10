---
layout: post
title: Installation
author: "Tobias Madsen"
tags: [tutorial]
---



<!--
Build post
library(knitr)
render_jekyll()
knit(input = "../gh-pages-dgRaph/knitr//installation.Rmd", output = "../gh-pages-dgRaph/_posts/2015-06-10-installation.md")
-->

First step install the `devtools` package:


{% highlight r %}
install.packages("devtools")
{% endhighlight %}

With `devtools` in place. Use `devtools` to install `dgRaph` directly from github:


{% highlight r %}
devtools::install_github("TobiasMadsen/dgRaph")
{% endhighlight %}

If compilation fails please raise an issue on [github](https://github.com/TobiasMadsen/dgRaph/issues).

