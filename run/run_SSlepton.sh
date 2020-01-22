#DATA
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2016 -n 20 &
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2017 -n 20 &
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2018 -n 20 &

python python/SKFlat.py -a SSlepton -i SingleMuon -y 2016 -n 20 --userflags RunNI &
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2017 -n 20 --userflags RunNI &
python python/SKFlat.py -a SSlepton -i SingleMuon -y 2018 -n 20 --userflags RunNI &
#MC
