The user manual is converted into a pdf called converterManual.pdf by typing

pandoc user_manual.md -o converterManual.pdf --pdf-engine=xelatex   -V monofont="DejaVu Sans Mono" -V mainfont="TeX Gyre Pagella" -H <(echo '\usepackage{fvextra}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,breakanywhere,commandchars=\\\{\}}')
