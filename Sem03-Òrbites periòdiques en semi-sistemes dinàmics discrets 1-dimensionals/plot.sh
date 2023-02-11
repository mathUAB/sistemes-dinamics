#!/bin/zsh
gcc -Wall -pedantic -o sem3 sem3.c -lm
./sem3 $1
returnvalue=$?
if [[ $returnvalue -eq 0 && $2 -eq 1 ]]; then
  echo "Computation ended successfully. Ploting..."
  gnuplot -p -e "
  set encoding utf8; 
  set terminal qt size 1750, 850;
  set xlabel 'µ'; 
  set ylabel 'x'; 
  set key left top; 
  plot 'data.dat' title 'x_0 = $1' pointtype 7 pointsize 0.001 linecolor 6"
  # plot [4:5.828427124] 'data.dat' title 'x_0 = $1' pointtype 7 pointsize 0.001 linecolor 2"
elif [[ $returnvalue -ne 0 ]];then
  echo "Error!"
fi