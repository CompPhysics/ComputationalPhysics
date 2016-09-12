#!/bin/sh
set -x

function system {
  "$@"
  if [ $? -ne 0 ]; then
    echo "make.sh: unsuccessful command $@"
    echo "abort!"
    exit 1
  fi
}

if [ $# -eq 0 ]; then
echo 'bash make.sh slides1|slides2'
exit 1
fi

name=$1
rm -f *.tar.gz

opt="--encoding=utf-8"
opt=

rm -f *.aux



# Plain HTML documents


html=${name}
system doconce format html $name --pygments_html_style=default --html_style=bloodish --html_links_in_new_window --html_output=$html $opt
system doconce split_html $html.html --method=space10


# IPython notebook
system doconce format ipynb $name $opt

# Ordinary plain LaTeX document
system doconce format pdflatex $name --print_latex_style=trac --latex_admon=paragraph $opt
system doconce ptex2tex $name envir=print
# Add special packages
doconce subst "% Add user's preamble" "\g<1>\n\\usepackage{simplewick}" $name.tex
doconce replace 'section{' 'section*{' $name.tex
pdflatex -shell-escape $name
pdflatex -shell-escape $name
mv -f $name.pdf ${name}.pdf
cp $name.tex ${name}.tex



