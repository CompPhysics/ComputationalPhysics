#!/bin/sh
# Create HTML version of quizzes
name=cpp

doconce format html $name --html_style=bootstrap --html_code_style=inherit --quiz_horizontal_rule=on
echo "load $name.html into a browser, e.g., google-chrome $name.html"

# If !split is used to split the document into multiple pages, run
#doconce split_html quiz.html

# PDF via latex including full solutions
#doconce format pdflatex $name --without_answers
#doconce ptex2tex $name envir=minted
#pdflatex -shell-escape $name
