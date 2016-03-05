#!/usr/bin/python
#-*-coding:utf-8-*-
from __future__ import print_function
#Tripoli evoluation calclation
#2015-1 -18 V0.2.1
from BaseTripoli import BaseTripoli
import os,sys
import re
import time
import lancejob


    
class TripoliE(BaseTripoli):
    _post_fix='.txt'
    isNewTable=False
    dt=3.
    #dt=0.1

    def __init__(self,project=None,ifile=None,file=None,deltaBu=None,density=None,maxSigma=None,T=None):
        super(TripoliE,self).__init__()
        if T:
            self.T=T
        self.jobId='NON'
        self.jobName='NewJob'
        self.setProjectName(project)
        self._log=self._proName+'.log'
        f=open(self._log,mode='w')
        f.writelines(time.strftime('%m-%d,%H:%M:%S'))
        f.close()
        
        self.setT0File(ifile)
        if T:
            self.T=T
        else:
            self.T=0
        if file:
            self.setTFile(file)
            self.writeLog("begin compute with T=",T)
        self.setDeltaBu(deltaBu)
        self.setDensity(density)
        self.setMaxSigma(maxSigma)

        
    def maketable(self,filename=None):
        self.readConfig()
        super(TripoliE,self).maketale(filename)
        

            


##########################################################
    ##设置项目的基本属性
    def setMaxSigma(self,sgm=None):
        self._err_sigma=sgm or 0.8
    def setDeltaBu(self,bu=None):# unit GWd/tU
        self._dBu=(bu or 3.)*1.e9*24.*3600./1000.  #unit j/kg
    def setDensity(self,den=None):#unit kg/m3
        self._rho=(den or 11007.22976)*1.e-6   #unit kg/cm3
    def setProjectName(self,name=None):
        self._proName=name or "MyProject"
    def setT0File(self,file):
        if file:
            self.T=0
            self._modelName=self._proName+'_base.tri'
            self._ifilename=self._proName+'_T'+str(self.T)+'.tri'
            self._postFile=self._ifilename+self._post_fix
            f=open(file)
            c=f.read()
            f.close()
            f=open(self._modelName,mode='w')
            f.write(c)
            f.close()
            f=open(self._ifilename,mode='w')
            f.write(c)
            f.close()
    def setTFile(self,file):
        self._modelName=self._proName+'_base.tri'
        self._ifilename=self._proName+'_T'+str(self.T)+'.tri'
        self._postFile=self._ifilename+self._post_fix
        self.maketable(file)
        c=self._input
        f=open(self._modelName,mode='w')
        f.write(c)
        f.close()


###########################################################

    def parsing_input(self):
        super(TripoliE,self).parsing_input()
        self._computeDeltaM()
        
    def _computeDeltaM(self):
        self._delta_M=0.0
        for v in self._comV:
            self._delta_M+=self.realVol[v]*self._rho
        return self._delta_M
###########################################################
    ##处理output部分
    #if self.input dosent existe, then creat it
    def isMakeTable(self):
        if self.isNewTable:
            self.maketable()
            self.isNewTable=False
        if not hasattr(self,'_input'):
            self.maketable()

    ##计算生成新input文件
    def _compute_burnup(self):
        self.isMakeTable()
        if not hasattr(self,'realVol'):
            self.parsing_input()
        energie=0.0
        for v in self.volume:
            taux_reaction=0.0
            for isotop in self.isotopes:
                taux_reaction+=self.score[v][33][isotop]
            energie+=taux_reaction*self.realVol[v]
        energie=energie*200*1.e6*1.6*1.e-19  #每次反应发出200MeV, 这里单位为w
        #计算烧完deltaBu 所需要的时间
        if not hasattr(self,'_delta_M'):
            self._computeDeltaM()
        self._deltaT=self._dBu*self._delta_M/energie #unit seconde
        return self._deltaT
    #计算下一个时刻点的核素丰度
    def _compute_concentration(self):
        '''
        compute model: we suppose that the reaction rate is a constant
        N5=N5-(5R52+5R33)*T
        N6=N6-(6R52+6R33-5R52)*T
        N7=N7-(7R52+7R33-6R52)*T
        NU8=NU8-(U8R52+U8R33)*T
        N8=N8-(8R52+8R33-7R52)*T
        N9=N9-(9R52+9R33-U8R52)*T
        N40=N40-(40R52+40R33-9R52)*T
        N41=N41-(41R52+41R33-40R52)*T
        N42=N42-(42R52+42R33-42R52)*T
        '''
        self.isMakeTable()
        self._compute_burnup()
        if not hasattr(self,'_comps'):
            self.parsing_input()
        #newComps format: self.newComps[vol][isotope]=concentration
        self._newComs=dict()
        delta_bu=self._deltaT*pow(10,-24)
        for v in self.volume:
            the_v=dict()
            for isotope in self.isotopes:
                try:
                    if isotope=='U235':
                        the_v[isotope]=self.concentration[v][isotope]-(self.score[v][33][isotope]+self.score[v][52][isotope])*delta_bu
                    if isotope=='U236':
                        the_v[isotope]=self.concentration[v][isotope]-(self.score[v][33][isotope]+self.score[v][52][isotope]-self.score[v][52]['U235'])*delta_bu
                    if isotope=='NP237':
                        the_v[isotope]=self.concentration[v][isotope]-(self.score[v][33][isotope]+self.score[v][52][isotope]-self.score[v][52]['U236'])*delta_bu
                    if isotope=='PU238':
                        the_v[isotope]=self.concentration[v][isotope]-(self.score[v][33][isotope]+self.score[v][52][isotope]-self.score[v][52]['NP237'])*delta_bu
                    if isotope=='U238':
                        the_v[isotope]=self.concentration[v][isotope]-(self.score[v][33][isotope]+self.score[v][52][isotope])*delta_bu
                    if isotope=='PU239':
                        the_v[isotope]=self.concentration[v][isotope]-(self.score[v][33][isotope]+self.score[v][52][isotope]-self.score[v][52]['PU238'])*delta_bu
                    if isotope=='PU240':
                        the_v[isotope]=self.concentration[v][isotope]-(self.score[v][33][isotope]+self.score[v][52][isotope]-self.score[v][52]['PU239'])*delta_bu
                    if isotope=='PU241':
                        the_v[isotope]=self.concentration[v][isotope]-(self.score[v][33][isotope]+self.score[v][52][isotope]-self.score[v][52]['PU240'])*delta_bu
                    if isotope=='PU242':
                        the_v[isotope]=self.concentration[v][isotope]-(self.score[v][33][isotope]+self.score[v][52][isotope]-self.score[v][52]['PU241'])*delta_bu
                except KeyError:
                    self.writeLog('vol=',v,'isotope=',isotope)
            self._newComs[v]=the_v
    #生成新的组分模块
    #默认了燃料温度为924k，慢化剂温度为574k
    #这样是为了保证数据库中拥有相应温度的数据
    def _generate_GeoCom(self):
        self._compute_concentration()
        com_begin='COMPOSITION '+str(len(self.concentration))+' \n'
        com_end='END_COMPOSITION \n'
        com_new='\n'
        for vol in self.concentration:
            if vol in self.volume:
                #燃料温度在这里修改：924
                com_new+='POINT_WISE 924 \t'+self.concentration[vol]['name']+' '+\
                    str(len(self.concentration[vol])-2)+'\n'
                for isotope in self.concentration[vol]:
                    if isotope!='name' and isotope!='num':
                        if isotope in self._newComs[vol]:
                            com_new+=isotope+'\t'+str(self._newComs[vol][isotope])+'\n'
                        else :
                            com_new+=isotope+'\t'+str(self.concentration[vol][isotope])+'\n'
            else:
                #慢化即温度在这里修改
                com_new+='POINT_WISE 574 \t'+self.concentration[vol]['name']+' '+\
                    str(len(self.concentration[vol])-2)+'\n'
                for isotope in self.concentration[vol]:
                    if isotope!='name' and isotope!='num':
                        com_new+=isotope+'\t'+str(self.concentration[vol][isotope])+'\n'
        comSection=com_begin+com_new+com_end
        return comSection
    #生成下一时刻的input文件名
    def _creat_ifilename(self):
        self.T=self.T+1
        self._ifilename=self._proName+'_T'+str(self.T)+'.tri'
        self._postFile=self._ifilename+self._post_fix
        return self._ifilename
    def generate_Inputfile(self):
        self.readConfig()
        self.writeLog("Now, max error sigma is:",self._err_sigma,
                      "deltaBU is:",self._dBu)
        comSection=self._generate_GeoCom()
        f=open(self._modelName)
        base=f.read()
        f.close()
        new=re.sub(r'COMPOSITION([\w\W]+?)END_COMPOSITION',comSection,base)
        self.writeLog('the current filename is :',self._postFile)
        T=self.T
        self._creat_ifilename()
        if T!=self.T-1:
            self.writeLog("Error!!")
        f=open(self._ifilename,mode='w')
        f.write(new)
        f.close()
        return new
##########################################################
    ##运行部分
    def run(self):
        try:
            self.isMakeTable()
            self.writeLog(self._postFile+'\n','Max sigma=',self.max_sigma)
            self.writeResult()
            if self.max_sigma<self._err_sigma:
                self.writeLog('lance a new job!')
                self.lancheNewJob()
                try:
                    self.writeResult()
                except IOError:
                    pass
        except IOError:
            self.lancheNewJob()
            pass
        print(self._ifilename)
        initial=time.time()
        time.sleep(3)
        while time.time()-initial <self.dt+10:
            if time.time()-initial>self.dt:
                try:
                    self.writeLog("try to make a table!")
                    self.maketable()
                    self.writeResult()
                    self.writeLog("make a table successfully !")
                    self.writeLog("max sigma is ",self.max_sigma)
                    if self.max_sigma<self._err_sigma:
                        self.lancheNewJob()
                        initial=time.time()
                    else:
                        initial=time.time()
                except IOError:
                    initial=time.time()
                    self.writeLog("wait for the Tripoli result!")
                    pass
    def lancheNewJob(self):
        if self.T!=0:
            self.generate_Inputfile()
        if self.T==0:
            try:
                self.generate_Inputfile()
            except IOError:
                pass
        self.writeLog('generate a new inputfile OK!\nThe new file is:'+self._ifilename)
        self.lastJobId=self.jobId
        self.lastJobName=self.jobName
        self.jobId,self.jobName=lancejob.generate(['','-d',self._ifilename])
        self.writeLog("jobId is "+str(self.jobId)+"\tjobName is :"+str(self.jobName))
        ##删除上一次作业
        qstat=os.popen('qstat').read()
        self.writeLog(qstat)
        if qstat.find(self.lastJobId)!=-1:
            qdel=os.popen('qdel '+self.lastJobId+'\n')
            self.writeLog(qdel)
        self.isNewTable=True
    def readConfig(self,configure='tripoli.conf'):
        f=open(configure)
        conf=f.read()
        f.close()
        conf=re.sub(r'//.*','',conf)
        conf=re.sub(r'/\*[\w\W]*?\*/','',conf)
        config=eval(conf)
        for i in config:
            if i['project']==self._proName:
                self.setDeltaBu(i['deltaBu'])
                self.setMaxSigma(i['maxSigma'])
                break
        
##########################################################
    ##结果输出部分
    #输出一般结果
    def writeResult(self):
        self.isMakeTable()
        table_title=self._proName+'\tBU='+str(self._dBu*self.T/(1.e9*24.*3600./1000))+'\n'+'isotope\tvol_number\treaction_type\tscore\tsigma\tbatch_nb\n'
        table=''
        for vol in self.total:
            for reaction in self.total[vol]:
                for isotope in self.total[vol][reaction]:
                    table+=isotope['nucleus']+'\t'+str(isotope['volume'])+'\t'+\
                        str(isotope['reaction'])+'\t'+str(isotope['score'])+'\t'+\
                        str(isotope['sigma'])+'\t'+str(isotope['batch'])+'\n'
        f=open(self._ifilename+'.res',mode='w')
        f.write(table_title+table)
        f.close()
        return(table)
    
