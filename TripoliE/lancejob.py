#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import re

tripoli_path = r'/home/yuancenxi/tripoli/Tripoli/CODE/bin/linux-intel-64/static_tripoli4 '
input_path = r"/home/Tripoli/input/"
utput_path = r"/home/Tripoli/output/"
cea_path = "/home/yuancenxi/tripoli/Tripoli/Env/t4path.ceav5"
endl_path = "/home/yuancenxi/tripoli/Tripoli/Env/t4path.endl"
cwd = os.getcwd()


def lance(content):
    file = open('./work.temp', mode='w+')
    file.write(content)
    file.close()
    print('sub a job\n')
    exe = os.popen("qsub work.temp")
    return exe.read()


def generate(argv):
    name = 'defaultname'
    lan = "english"
    para = ' '
    lib = cea_path
    mode = "NJOY"
    pe = '2'
    host = 'default_host'
    hostname = ' '
    outputf = 'defaultopt'
    for i in range(1, len(argv)):
        if argv[i] == "-d":
            inputf = argv[i + 1]
            if name == "defaultname":
                lst = inputf.rfind('/')
                if lst != -1:
                    name = inputf[lst + 1:]
                else:
                    name = inputf
            if outputf == 'defaultopt':
                outputf = inputf + '.txt'
        if argv[i] == "-o":
            outputf = argv[i + 1]
        if argv[i] == "-l":
            lan = argv[i + 1]
        if argv[i] == "-p":
            para = argv[i + 1] + ' -t bsd '
        if argv[i] == "-c":
            lib = argv[i + 1]
        if argv[i] == "-s":
            mode = argv[i + 1]
        if argv[i] == "-N":
            name = argv[i + 1]
        if argv[i] == "-pe":
            pe = argv[i + 1]
        if argv[i] == "-host":
            host = argv[i + 1]
    parameters = 'mpirun -np  ${NSLOTS}  ' + tripoli_path + ' -d ' + inputf + ' -i ' + pe + ' -s ' + mode + ' -c ' + lib + ' -o ' + outputf + ' -l ' + lan + para + ' –TABPROB  \n'
    if host != 'default_host':
        hostname = "#$ -l hostname=compute-0-" + host + '\n'
    jobname = '#$ -N ' + name + '\n'
    mutliP = '#$ -pe mpich ' + pe + '\n'
    preEviroment = '#! /bin/bash  \n' \
                   '#$ -S /bin/bash \n' \
                   '#$ -cwd \n' \
                   '#$ -V  \n' \
                   '#$ -o my.opt  \n' \
                   '#$ -e my.err  \n'
    content = preEviroment + jobname + mutliP + hostname + parameters
    info = lance(content)
    print(info)
    pattern = r'Your job (\d+) ([\w\W]+?)has been submitted'
    ##jobId,jobName=re.findall(pattern,info)[0]
    # print(jobId,jobName)
    return '0', '1'


if __name__ == "__main__":
    generate(sys.argv)
