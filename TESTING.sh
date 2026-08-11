#!/bin/bash
time_limit=1800
maxit=100
maxit_tabu=100

#for layers in "50"
#do

   for k in {1..3}
   do

#      for f in instance/${layers}-layers/*
#      do

	f="instance/incgraph_6_0.15_15_20_2.00_2.txt"

         name="${f%.*}"
         name="${name#*/}"
         name="${name#*/}"

         echo "instance:$name - k_$k"

      algorithm=grasp3	
	   
      alpha=0;

      ls="best"

      algorithm=grasp3	
      ./Release/C-IGDP $f $algorithm $k $alpha $ls $maxit $name $name.txt $name\_$k.txt $time_limit

#     done
   done
#done
