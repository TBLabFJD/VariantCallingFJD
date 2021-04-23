
printf '\n\nFiltering by population frequencies...\n'

perl $FILTER_VEP \
-i $VCF_FILTERED -o $VCF_FILTERED_2 \
--filter "(DP > 10) and (MAX_AF < 0.01 or not MAX_AF) and (CANONICAL is YES) and (Consequence is 3_prime_UTR_variant or Consequence is 5_prime_UTR_variant or Consequence is intron_variant or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is splice_region_variant or Consequence is synonymous_variant or Consequence is missense_variant or Consequence is inframe_deletion or Consequence is inframe_insertion or Consequence is stop_gained or Consequence is frameshift_variant or Consequence is coding_sequence_variant)" \
--force_overwrite 


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\n VEP FREQUENCY FILTERING 1 for '${sample}' DONE\n' 
	rm $VCF_FILTERED

else
	printf "\nERROR: PROBLEMS WITH VEP FREQUENCY FILTERING 1"
	exit 1
fi


perl $FILTER_VEP \
-i $VCF_FILTERED_2 -o $VCF_FILTERED_3 \
--filter "MAX_AF < 0.01 or not MAX_AF" \
--force_overwrite 


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP FREQUENCY FILTERING 2 for '${sample}' DONE\n' 
	rm $VCF_FILTERED_2

else
	printf "\nERROR: PROBLEMS WITH VEP FREQUENCY FILTERING 2"
	exit 1
fi

