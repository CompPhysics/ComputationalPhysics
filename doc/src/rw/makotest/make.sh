name=rw
# Textbook with Python code
doconce format html $name --html_style=bootstrap CODE=Python HEAD_TYPE=0 --html_output=${name}-text-py
# Textbook with C++ code
doconce format html $name --html_style=bootstrap CODE=C++ HEAD_TYPE=0 --html_output=${name}-text-cpp
# Slides with Python code
doconce format html $name --html_style=bootstrap CODE=Python HEAD_TYPE=1 --html_output=${name}-slides-py


