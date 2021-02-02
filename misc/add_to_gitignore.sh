# Script for adding large files to gitignore
cat .gitignore > .gitignore_tmp
echo ".Rproj.user 
.Rhistory 
.RData 
.Ruserdata" >> .gitignore_tmp
find . -type f -size +2M >> .gitignore_tmp
sed -i 's/^\.\///' .gitignore_tmp  
awk '!seen[$0]++' .gitignore_tmp > .gitignore
rm .gitignore_tmp
