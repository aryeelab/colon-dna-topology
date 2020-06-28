cd ../terra

# Get full file listing
gsutil ls -lr gs://fc-74ac2103-b976-4792-9857-9deb489ba149 > terra_files_2020-05-15.txt
# Size
awk '{sum+=$1/(1024^4);} END{print sum;}' terra_files_2020-05-15.txt 
# 35.488 TB


# Select files >1Mb
awk '{ if (NF==3 && $1 > 1024^2) { print } }' terra_files_2020-05-15.txt > terra_files_2020-05-15_1mb.txt
# 19747 files
# Size:
awk '{sum+=$1/(1024^4);} END{print sum;}' terra_files_2020-05-15_1mb.txt
# 35.4876 TB

# Run terra_files_to_keep.Rmd to generate list of files to keep from the sample and sample_set metadata tables

# Get list of files >1mb that are NOT listed in the metadata tables
cat terra_files_2020-05-15_1mb.txt | grep -v -f terra_files_to_keep.txt > terra_files_to_delete_2020-05-15.txt
# Size 
awk '{sum+=$1/(1024^4);} END{print sum;}' terra_files_to_delete_2020-05-15.txt
# 30.5809 TB

# Move delete file list to the bucket so we can run delete from a Google cloud VM
# gsutil cp terra_files_to_delete_2020-05-15.txt gs://fc-74ac2103-b976-4792-9857-9deb489ba149/
awk '{print $3}' terra_files_to_delete_2020-05-15.txt | gsutil -m rm -I


# Storage before and after delete
# Before delete
gsutil du -ahs gs://fc-74ac2103-b976-4792-9857-9deb489ba149
# 35.49 TiB   
# $1014.51/month

# After delete
# 4.91 TiB
# $140.28/month