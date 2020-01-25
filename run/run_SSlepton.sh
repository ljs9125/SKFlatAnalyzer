#DATA
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2016 -n 4 &
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2017 -n 4 &
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2018 -n 4 &

python python/SKFlat.py -a SSlepton -i SingleMuon -y 2016 -n 4 --userflags RunNI &
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2017 -n 4 --userflags RunNI &
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2018 -n 4 --userflags RunNI &
#MC
