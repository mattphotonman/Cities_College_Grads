#!/bin/bash
#Finds the line numbers of the cities from city_names_from_college_percent.txt in city_population_wiki_2010_census.txt.
#For each city in city_names_from_college_percent.txt writes this city name, and then the candidate lines and line numbers that could be this city from city_population_wiki_2010_census.txt.

n=`cat city_names_from_college_percent.txt | wc -l`
touch candidate_line_numbers.txt
rm candidate_line_numbers.txt
for ((linenum=1;linenum<=n;linenum++)); do
    echo "$linenum `print_line $linenum city_names_from_college_percent.txt`" >> candidate_line_numbers.txt
    cp city_population_wiki_2010_census.txt tmp
    for word in $( print_line $linenum city_names_from_college_percent.txt ); do
	grep $word tmp > tmp2
	mv tmp2 tmp
    done
    cat tmp >> candidate_line_numbers.txt
    echo >> candidate_line_numbers.txt
done
rm tmp
