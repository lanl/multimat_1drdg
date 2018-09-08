#
# This script cleans all the .dat and .plt files in the testcases directory
#

count=1

cd ./testcases

for d in ./*
do
  echo $count : $d
  cd "$d" && rm *.dat .DS_Store
  cd ..
  count=$((count+1))
done

echo Cleanup complete.
