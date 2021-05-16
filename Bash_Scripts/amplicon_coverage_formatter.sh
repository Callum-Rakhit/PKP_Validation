##### Get the coverage per primer pair metrics #####

# Copy, rename and split the amplicon coverage metrics file
function aggregate_amplion_coverage {
  for i in $(ls $1); do for x in $(ls $1/$i); 
  do cp $1/$i/$x/metrics/*.counts.PRIMER $2 &&
  for y in $1/$i/$x/metrics/*.counts.PRIMER;
  do z=$(basename $y) && a=$(echo $z | cut -c1-8) &&
  mv $2/$a* $2/${i}_${x}_amplicon_coverage_temp &&
  awk '{print $1, $3}' $2/${i}_${x}_amplicon_coverage_temp | sed '1,6d' | head -n -1 > $2/${i}_${x}_amplicon_coverage;
  done;
  done;
  done
}

# Final function to call all previous functions by iterating for the sample folders
function make_plots {
  aggregate_coverage_files $1 $2
  aggregate_amplion_coverage  $1 $2
  aggregate_VAF $1 $2
  # rm *temp* && rm *verbose* && rm *metrics*
}

make_plots ./analysis/20191003_190447/ ./metrics_extraction_for_validation/