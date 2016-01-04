import os
src_path = 'src'

# Mako function for flexible subsection headings
def ssh(heading, onoff):
    return '===== ' + heading.strip() + ' =====' if (HEAD_TYPE == 'slides' or (HEAD_TYPE == 'textbook' and onoff == 1)) else ''

# Here is how to edit a slides file with lots of subsection headings:
# doconce subst -m '^===== (.+?) =====' '${ssh("\g<1>", 1)}' *.do.txt
# Then edit 1 to 0 in ssh calls where you don't want the heading in a
# textbook-like version of the document.
# Run doconce format html mydoc HEAD_TYPE=slides
# to get all headings, else set HEAD_TYPE=textbook to rule out those with
# the 0 argument

# Mako function for flexible inclusion of code in different languages.
# Programs are located in src/py, src/cpp, src/f.
# Set CODE=Python/Fortran/C++ on the command line to select language.

def code(code='', filename='', language=None, from_regex=None, to_regex=None):
    # code can be a filename or computer code
    if language is None:
        language = CODE  # Use global language if not specified
    if language == 'Python':
        if filename:
            filename += '.py'
            # Include from file
            text = '@@@CODE src/py/%s' % filename
            if from_regex is not None and to_regex is not None:
                # Include just a portion of the file
                text += ' fromto: %s@%s' % (from_regex, to_regex)
            elif from_regex is not None and to_regex is None:
                text += ' fromto: %s@' % (from_regex)
        else:
            # The code argumnet holds the actual computer code,
            # assume it's just a code snippet (not complete program)
            text = '!bc pycod\n%s\n!ec' % code.strip()
    elif language == 'Fortran':
        if filename:
            for ext in '.f', '.f90':
                if os.path.isfile(os.path.join(src_path, 'f', filename + ext)):
                    filename += ext
                    break
            # Include from file
            text = '@@@CODE src/f/%s' % filename
            if from_regex is not None and to_regex is not None:
                # Include just a portion of the file
                text += ' fromto: %s@%s' % (from_regex, to_regex)
            elif from_regex is not None and to_regex is None:
                text += ' fromto: %s@' % (from_regex)
        else:
            # The code argumnet holds the actual computer code,
            # assume it's just a code snippet (not complete program)
            text = '!bc fcod\n%s\n!ec' % code.strip()
    elif language == 'C++':
        if filename:
            for ext in '.cpp', '.c++', '.cxx':
                if os.path.isfile(os.path.join(src_path, 'cpp', filename + ext)):

                    filename += ext
                    break
            # Include from file
            text = '@@@CODE src/cpp/%s' % filename
            if from_regex is not None and to_regex is not None:
                # Include just a portion of the file
                text += ' fromto: %s@%s' % (from_regex, to_regex)
            elif from_regex is not None and to_regex is None:
                text += ' fromto: %s@' % (from_regex)
        else:
            # The code argumnet holds the actual computer code,
            # assume it's just a code snippet (not complete program)
            text = '!bc cppcod\n%s\n!ec' % code.strip()
    else:
        print 'language=%s is illegal' % language
    return text
