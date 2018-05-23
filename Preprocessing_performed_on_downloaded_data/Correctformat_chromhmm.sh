gunzip *.gz

for file in *.bed
do
awk '{gsub(/^chr/,""); print}' $file | bgzip > ${file}.gz
rm $file
done

for file in *.bed.gz
do
tabix -p bed $file
done
