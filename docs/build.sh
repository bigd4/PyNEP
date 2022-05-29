rm -rf build/* html
# ls source/*.rst | grep -v "index.rst" | xargs rm -rf
sphinx-apidoc -f -o source/ ../pynep/
make html
cp -r build/html .
