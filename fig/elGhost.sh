

for f in $(ls *.pdf); do
  
  pdffile="$f"
    
  echo "$pdffile $f"
  ghostscript -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE -dQUIET -dBATCH -sOutputFile="../$f" "$f"  

done
