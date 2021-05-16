# Get information about total reads and % ROI covered, plus headers
for i in $(ls /home/callum/PKP_validation/PKP_metrics_files/*.metrics.exoncoverage); 
do grep -A 2 "# COVERAGE BINS" $i | tail -1 >> exoncoverage.summary;
done

# Remove the last column (extra addition for the subsamples)
# sed -i -r 's/\S+//19' exoncoverage.summaryplus 

for i in $(ls /home/callum/PKP_validation/PKP_metrics_files/*.metrics.READS); 
do grep -A 2 "# BASIC FILTERING METRICS" $i | tail -1 >> READS.summary;
done

# Convert names to lower case
find . -name '*.*' -exec sh -c 'a=$(echo "$0" | sed -r "s/([^.]*)\$/\L\1/"); [ "$a" != "$0" ] && mv "$0" "$a"' {} \;

# Get the duplication rates from the .json
for i in $(ls /home/callum/PKP_validation/PKP_metrics_files/*.json); 
# do grep -A 1 "READS_OUTPUT_FRAC" $i | awk '{print $2}' | cut -d \, -f 1 | awk -F '"TAS"' '{print $1}' | awk 'NF' >> duplication.summary; # for UMI and PCR duplicates
do cat $i | grep '"READS_OUTPUT_FRAC":' | cut -c 41-46 | cut -d \, -f 1 >> duplication.summary; # Just PCR
done

# Get the duplication rates from the .json
for i in $(ls /home/callum/PKP_validation/PKP_metrics_files/*.json); 
# do grep -A 1 "READS_OUTPUT_FRAC" $i | awk '{print $2}' | cut -d \, -f 1 | awk -F '"TAS"' '{print $1}' | awk 'NF' >> duplication.summary; # for UMI and PCR duplicates
do cat $i | grep '"READS_OUTPUT_FRAC":' | cut -c 41-46 | cut -d \, -f 1 >> duplication.summary; # Just PCR
done