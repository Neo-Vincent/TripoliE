from __future__ import print_function
from TripoliE import TripoliE

import threading
import re
class TripoliHandle(threading.Thread):  
    def __init__(self,config):
        super(TripoliHandle,self).__init__()
        print(config)
        self.Tripoli=TripoliE(**config)
    def run(self):
        self.Tripoli.run()
def runJobs(configure='tripoli.conf'):
    f=open(configure)
    conf=f.read()
    f.close()
    conf=re.sub(r'//.*','',conf)
    conf=re.sub(r'/\*[\w\W]*?\*/','',conf)
    config=eval(conf)
    for i in config:
        job=TripoliHandle(i)
        job.setDaemon(True)
        job.start()
    job.join()
        
        
runJobs()
