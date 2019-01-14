#!/bin/bash

julia make.jl


# Disabling the pretty-url feature of `makedocs` doesn't work, so we have to
# revert it manually

sed -i -e "s/href=\"aux\/\"/href=\"aux\/index.html\"/g" build/index.html
sed -i -e "s/href=\"io\/\"/href=\"io\/index.html\"/g" build/index.html
sed -i -e "s/href=\"model\/\"/href=\"model\/index.html\"/g" build/index.html
sed -i -e "s/href=\"search\/\"/href=\"search\/index.html\"/g" build/index.html

sed -i -e "s/href=\"..\/\"/href=\"..\/index.html\"/g" build/*/index.html
sed -i -e "s/href=\"..\/aux\/\"/href=\"..\/aux\/index.html\"/g" build/*/index.html
sed -i -e "s/href=\"..\/io\/\"/href=\"..\/io\/index.html\"/g" build/*/index.html
sed -i -e "s/href=\"..\/model\/\"/href=\"..\/model\/index.html\"/g" build/*/index.html
sed -i -e "s/href=\"..\/search\/\"/href=\"..\/search\/index.html\"/g" build/*/index.html
