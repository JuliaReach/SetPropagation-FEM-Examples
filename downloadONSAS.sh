

if [[ -d "ONSAS" ]]; then
  echo "pulling ONSAS repo"
  cd ONSAS
  git pull
  cd ..
else
  echo "cloning ONSAS repo"
  git clone https://github.com/ONSAS/ONSAS.m.git ONSAS
fi
