# Script for adding large files to gitignore
find . -type f -size +10M >> .gitignore
