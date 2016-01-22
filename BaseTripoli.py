#!/usr/bin/python
# -*- coding: utf-8 -*-
# Tripoli 基类，用来处理Tripoli输出数据
from __future__ import print_function
import re
import time

class BaseTripoli(object):
    def __init__(self,postFile=None):
        #科学计数法
        self._science_num=r'\d\.?\d*[eE][\-\+]\d+'
        #Tripoli输出文件
        self._postFile=postFile
        #本程序的日志文件
        self._log=(self._postFile or 'Myproject')+'.log'
        #出现错误休眠时间
        self._dt=3
    
    def _rmComment(self,content):
        #去除Tripoli的注释
        cont=re.sub(r'//.*','',content)
        cont=re.sub(r'/\*[\w\W]*?\*/','',cont)
        return cont
    def writeLog(self,*args):
        log=open(self._log,mode='a')
        print(time.strftime('%m-%d,%H:%M:%S')+':',*args, sep=' ', end='\n', file=log)

    ##更改数据格式
    def _V_R(self):
    #creat a self.total attribute to store the result in the format:
    # self.total[volume][reaction]=list of the result of the reation and in the voulume
    # self.score[volume][reaction][isotope]=score
    #volume and reaction type
        self.volume=list({i['volume'] for i in self._dl})
        self.reaction=list({i['reaction'] for i in self._dl})
        self.isotopes=list({i['nucleus'] for i in self._dl})
        self.volume.sort()
        self.reaction.sort()
        self.isotopes.sort()
        self.total={i:{j:[dd for dd in self._dl \
                          if dd['volume']==i and j==dd['reaction']] \
                       for j in self.reaction} for i in self.volume}
        self.score={i:{j:{dd['nucleus']:dd['score'] for dd in self._dl 
                          if dd['volume']==i and j==dd['reaction'] }
                       for j in self.reaction}\
                    for i in self.volume}
        self.spectrum={i:{j:{dd['nucleus']:dd['spectrum'] for dd in self._dl 
                          if dd['volume']==i and j==dd['reaction'] }
                       for j in self.reaction}\
                    for i in self.volume}
    #处理Tripoli输出信息
    def maketable(self,filename=None):
        if filename:
            self._postFile=filename
        f=open(self._postFile)
        content=f.read()
        f.close()
        science_num=self._science_num
        d_l=[]
        try:
            
            #pattern for the input info
            input_pattern=r'checking association of compositions and volumes :  ok([\w\W]+)\s+Loading response functions'
            self._input=re.findall(input_pattern,content)[0]
            self._input=self._rmComment(self._input)
            #content for lastest batch
            pattern1=r'(RESULTS ARE GIVEN FOR SOURCE INTENSITY[\W\w]+?simulation time)'
            content1=re.findall(pattern1,content)[-1]
            
            #block
            mod=r'(RESPONSE FUNCTION : REACTION\s+[\w\W]+?number of batches used:\s\d+\s\d\.\d+e[\-\+]\d+\s\d\.\d+e[\-\+]\d+)\s+'
            block=re.findall(mod,content1)
            #pattern
            p_nucleus=r'reaction on nucleus : ([A-Z]+\d+)\s'
            p_volume=r'num of volume : (\d+)\s'
            p_reaction=r'reaction consists in codes : \n\t\t(\d+)'
            p_batch=r'number of batches used: (\d+)\t'
            p_score2=r'\s+ENERGY INTEGRATED RESULTS\s+number of first discarded batches :\s+\d+\s+number of batches used: \d+\s+('+science_num+r')\s+('+science_num+r')'
            p_score=r'lethargy\s+'+ science_num+r'\s\-\s'+science_num+r'\t('+science_num+ r')\t'+science_num
            p_sigma=r'lethargy\s+'+ science_num+r'\s\-\s'+science_num+r'\t'+science_num+ r'\t('+science_num+r')'
            p_spectrum=r'SPECTRUM RESULTS([\w\W]+)?ENERGY INTEGRATED RESULTS'
            p_sp=r'('+science_num+r')\s+\-\s+('+science_num+r')\s+('+science_num+r')\s+('+science_num+r')\s+('+science_num+r')\s+'
            #dictionary
            for blck in block:
                d={}
                if blck.find('H1_H2O')!=-1:
                    d['nucleus']='H1_H2O'
                else:
                    d['nucleus']=re.findall(p_nucleus ,blck)[0]
                d['batch']=int(re.findall(p_batch ,blck)[0])
                d['volume']=int(re.findall(p_volume,blck)[0])
                d['score']=float(re.findall(p_score2,blck)[0][0])
                d['reaction']=int( re.findall(p_reaction,blck)[0])
                d['sigma']=float( re.findall(p_score2,blck)[0][1])
                sp_blck=re.findall(p_spectrum,blck)[0]
                spectrum=[]
                sp=re.findall(p_sp,sp_blck)
                for i in sp:
                    thesp={}
                    sp_re=i
                    thesp['group']=(float(sp_re[0]),float(sp_re[1]))
                    thesp['score']=float(sp_re[2])
                    thesp['sigma']=float(sp_re[3])
                    thesp['score//lethargy']=float(sp_re[4])
                    spectrum.append(thesp)
                d['spectrum']=spectrum
                d_l.append(d)
        except IndexError:
            print(self._postFile)
            print("INDEXERROR: is the file a coorect Tripoli RESULT file?")
            time.sleep(self._dt)
            self.maketable()
            return self.max_sigma
        Sigma={dd['sigma'] for dd in d_l}
        self._dl=d_l
        self.max_sigma=max(Sigma)
        self._V_R()
        return self.max_sigma

    def parsing_input(self):
        #parsing the input infomation
        #creat the attribute self.comps which store the all compostions in a list
        #self.geoCom store the geoCom info

        if not hasattr(self,'_input'):
            self.maketable()
        #get the relation between compositions and volums
        self._input=self._rmComment(self._input)
        geoCom_pattern=r'GEOMCOMP([\w\W]+?)END_GEOMCOMP'
        geostr=re.findall(geoCom_pattern,self._input)[0]
        geoCom_line=r'(\w+)\s+(\d+\s+[\d\s]+)'
        geoC_0=re.findall(geoCom_line,geostr)
        self.geoCom=[]
        for i in geoC_0:
            theGeo={'name':i[0]}
            nums=re.findall(r'(\d+)',i[1])
            theGeo['num']=int(nums[0])
            theGeo['volume']=int(nums[1]) # there is only one volume
            self.geoCom.append(theGeo)

        #获取composition部分
        #get composition section
        comps_pattern=r'COMPOSITION\s+\d+\s+([\w\W]+?)END_COMPOSITION'
        comps_str=re.findall(comps_pattern,self._input)[0]
        comps_0=comps_str.split('POINT_WISE')[1:]
        isotope_pattern=r'([A-Z]+\-??\w*?)\s+'+r'(\d\.?\d*[eE]{0,1}[\-\+]{0,1}\d+)'
        self._comps=[]
        self._mat_name=set()
        ##self._comps 数据格式 list，元素是一个dict，每个dict里面是一种材料的所有信息
        for i in comps_0:
            theComp={'name':re.findall(r'\d+\s+(\w+)\s+\d+?',i)[0]}
            self._mat_name.add(theComp['name'])
            i=i.replace(theComp['name'],'')
            isotopes=re.findall(isotope_pattern,i)
            for isotope in isotopes:
                theComp[isotope[0]]=float(isotope[1])
            self._comps.append(theComp)
        for i in range(len(self.geoCom)):
            for j in self._comps:
                if self.geoCom[i]['name']==j['name']:
                    self.geoCom[i].update(j)

        #reformat the geo_comps in the format self.concentration[vol][isotope]=concentration
        self.concentration={}
        for i in self.geoCom:
            this_v={}
            the_vol=-1
            for info in i:
                if info=='volume':
                    the_vol=i[info]
                else:
                    this_v[info]=i[info]
            self.concentration[the_vol]=this_v

        #get the vol number with the material name
        self._mat_vol={}
        for i in self.geoCom:
            for mat in self._mat_name:
                if i['name']==mat:
                    self._mat_vol[mat]=i['volume']
                    self._mat_vol[i['volume']]=mat

        #自动分析成分，获取燃料的体积
        self._comV={self._mat_vol[i['name']]  for i in self._comps\
                    if ('U235' in i or 'U238' in i  or 'PU239' in i) }
        #get the real volume in cm3
        #获取体积的真实体积大小
        real_V_pattern=r'VOLSURF([\w\W]+?)END_VOLSURF'
        try:
            vol_str=re.findall(real_V_pattern,self._input)[0]
            vol_pattern=r'(\d+)\s+(\d+\.\d+)?('+ self._science_num +')?\s+'
            vols=re.findall(vol_pattern,vol_str)
            self.realVol={int(vol[0]):float(vol[1]+vol[2]) for vol in vols}
        except IndexError:
            msg="无法找到体积参数，请手工添加体积参数在VOLSURF END_VOLSURF模块中\n",\
                "Cannot find the volume paramter, please add the paramter into the [VOLSURF END_VOLSURF] module manually.\n",\
                "The program will abort!!\n",\
                "Entre any key to exit!!\n"
            print(msg)
            self.writeLog(msg)
            exit(0)
