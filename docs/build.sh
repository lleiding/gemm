#!/bin/bash

echo "Preprocessing..."

# include a date stamp with the latest commit
sed -i -e "s/\*Last updated:.*/\*Last updated: $(git log --format="%cd (commit %h)" --date=short -1)\*  /" src/index.md 

echo "Building documentation..."

julia make.jl


# Disabling the pretty-url feature of `makedocs` doesn't work, so we have to
# revert it manually

echo "Postprocessing..."

sed -i -e "s/processes\//processes\/index.html/g" build/index.html build/search_index.js
sed -i -e "s/io\//io\/index.html/g" build/index.html build/search_index.js
sed -i -e "s/extensions\//extensions\/index.html/g" build/index.html build/search_index.js
sed -i -e "s/framework\//framework\/index.html/g" build/index.html build/search_index.js
sed -i -e "s/search\//search\/index.html/g" build/index.html build/search_index.js

sed -i -e "s/href=\"..\/\"/href=\"..\/index.html\"/g" build/*/index.html
sed -i -e "s/..\/processes\//..\/processes\/index.html/g" build/*/index.html
sed -i -e "s/..\/io\//..\/io\/index.html/g" build/*/index.html
sed -i -e "s/..\/extensions\//..\/extensions\/index.html/g" build/*/index.html
sed -i -e "s/..\/framework\//..\/framework\/index.html/g" build/*/index.html
sed -i -e "s/..\/search\//..\/search\/index.html/g" build/*/index.html

echo "Done."
