#!/usr/bin/python
# -*- coding: utf-8 -*-
# Tripoli 数据 处理类
from __future__ import print_function
#For parsing the tripoli4 output file
#2015-1 -18 V0.2.1
from BaseTripoli import BaseTripoli
import os,sys
import re
import time


    
class TripoliPost(BaseTripoli):
    _post_fix='.txt'
    #dt=0.1

    def __init__(self,file=None,deltaBu=None,density=None):
        super(TripoliPost,self).__init__()
        self._postFile=file
        self._log=file+'.log'
        f=open(self._log,mode='w')
        f.writelines(time.strftime('%m-%d,%H:%M:%S'))
        f.close()
        self._dBu=deltaBu or 3.0 #unit MWd/tu
        self._rho=density or 11007.22976*1.e-6 #unit kg/cm3
        try:
            self.T=int(re.findall(r'T(\d+).',file)[0])
        except IndexError:
            self.T=0

    def parsing_input(self):
        super(TripoliPost,self).parsing_input()
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
    def _computeLBU(self):
        self._compute_burnup()
        self.localBu={}
        for v in self.volume:
            taux_reaction=0.0
            for isotop in self.isotopes:
                taux_reaction+=self.score[v][33][isotop]
            energie=taux_reaction*self.realVol[v] 
            energie=energie*200*1.e6*1.6*1.e-19 #unit w
            self.localBu[v]=energie*self._deltaT/self._delta_M
    def generateTable(self):
        '''
        Table 的形式
        time:  T0
        volume   v1  v2  v3  v4  v5  v6
        reaction rate33 
        reaction rate 52
        concentation
        localBu
        '''
        if not hasattr(self,'_input'):
            self.maketable()
        if not hasattr(self,'realVol'):
            self.parsing_input()
        if not hasattr(self,'_deltaT'):
            self._computeLBU()
        table_title=self._postFile+'\tT='+str(self.T)+'\tFule density:\t'+str(self._rho)+'kg/cm3\n'
        table=str()
        volume=list(self.volume)
        rows=len(volume)
        table+='volume\tnumber\t'
        table+='%d\t'*rows%tuple(volume)+'\n'
        table+='volume\t cm3\t'
        table+='%E\t'*rows%tuple([self.realVol[vol] for vol in volume])+'\n\n'
        for reaction in self.reaction:
            table+='reation '+str(reaction)
            for isotope in self.isotopes:
                table+='\t'+isotope+'\t'+'%E\t'*rows%tuple([self.score[vol][reaction][isotope] for vol in volume])+'\n'
        table+='\nconcentration'
        for isotope in self.isotopes:
            table+='\t'+isotope+'\t'+'%E\t'*rows%tuple(self.concentration[vol][isotope] for vol in volume)+'\n'
        table+='\nLocalBU\tMWd/U\t'+'%E\t'*rows%tuple([self.localBu[vol] for vol in volume])+'\n\n\n\n'
        return table_title+table
            


    
