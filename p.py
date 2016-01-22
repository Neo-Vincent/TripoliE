from os import listdir
from postProcess import TripoliPost
def postData(cdir='.'):
    af=listdir(cdir)
    import pathlib
    postfiles=[file for file in af if pathlib.PurePath(file).suffix=='.txt']
    postfiles.sort()
    postdata=[]
    for i in postfiles:
        a=TripoliPost(i)
        postdata.append((i,a.T,a.generateTable()))
        print(i)
    postdata.sort()
    names=['H','L','M']
    for file in names:
        f=open(file+'.dat','w')
        context=str()
        for i in postdata:
            if i[0].find(file)!=-1:
                context+=i[2]
        f.write(context)
        f.close()
        
postData()
