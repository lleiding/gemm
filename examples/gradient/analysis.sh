tail -n 1 2019-0*_sg*/*tsv | grep -E '^1000' | cut -f 103 | sort | uniq -c | grep -vw 1 > twothree
