import json
import sys

utg2reads={}
with open(sys.argv[1]) as fr:
    for line in fr:
        if line.startswith('a'):
            a=line.split()
            read=':'.join(a[3].split(':')[:-1])
            if a[1] in utg2reads:
                utg2reads[a[1]]=utg2reads[a[1]]+' '+read
            else:
                utg2reads[a[1]]=read

with open('utg2reads.json','w') as fw:
    fw.write(json.dumps(utg2reads))


