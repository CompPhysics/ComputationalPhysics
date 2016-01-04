# Test code.py
import code

def test_ssh():
    heading  = 'Demo of code in different languages; Python'
    formatted_heading = '===== ' + heading + ' ====='
    code.HEAD_TYPE = 'textbook'
    assert code.ssh(heading, 1) == formatted_heading
    assert code.ssh(heading, 0) == ''
    code.HEAD_TYPE = 'slides'
    assert code.ssh(heading, 1) == formatted_heading
    assert code.ssh(heading, 0) == formatted_heading

def test_code():
    code.src_path = '..'  # relative path for src subdir
    code.CODE = 'Python'
    computed = code.code(filename='apb')
    expected = '@@@CODE src/py/apb.py'
    assert computed == expected
    computed = code.code(filename='apb', from_regex='a =', to_regex='a + b')
    expected = '@@@CODE src/py/apb.py fromto: a =@a + b'
    assert computed == expected
    code.CODE = 'C++'
    computed = code.code(filename='apb')
    expected = '@@@CODE src/cpp/apb.cpp'
    assert computed == expected
    computed = code.code(filename='apb', from_regex='a =', to_regex='a + b')
    expected = '@@@CODE src/cpp/apb.cpp fromto: a =@a + b'
    assert computed == expected
    computed = code.code(filename='demo', language='C++')
    text = """
#include <iostream>
using namespace std;

int main()
{
  a = 1;
  b = 2;
  cout << a + b;
  return 0;
}
"""
    computed = code.code(language='C++', code=text)
    expected = '!bc cppcod\n%s\n!ec' % text.strip()
    assert computed == expected

test_ssh()
test_code()
