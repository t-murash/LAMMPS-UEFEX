#!/bin/bash
INI=0
END=100
for((i=$INI;i<=$END;i++)); do
    sed -i "1c step=$i" anim.py
    python anim.py
done

convert -delay 5 -loop 0 figure.*.png movie.gif
