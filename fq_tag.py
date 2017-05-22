#!/usr/bin/env python3
# coding: utf-8

import os,sys,re,hashlib
import gzip
import getopt


def xopen_fq(file,mode='r'):
    #sub of check_fastq
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    else:
        return open(file,mode)

def Usage():
    print("""
        --indir --id --fq1 --fq2
""")
    exit(1)

def PAR(argv):
    Pardic = {}
    opts, args = getopt.getopt(argv[1:], 'h',['indir=','id=','fq1=','fq2='])
    for o,a in opts:
        if o in ('-h'):
            print("indir id fq1 fq2")
            exit()
        elif o in ('--indir'):
            Pardic.update({'indir':a})
        elif o in('--id'):
            Pardic.update({'id':a})
        elif o in ('--fq1'):
            Pardic.update({'fq1':a})
        elif o in ('--fq2'):
            Pardic.update({'fq2':a})
        else:
            print('unhandled option')
            exit(1)
    return Pardic

def find_sample(Pardic):
    fq1 = []
    all = os.listdir(Pardic['indir'])
    img1 = re.compile(Pardic['id']+".*R1.fastq.gz$"+"|"+Pardic['id']+".*R1.*fastq")
    for i in all:
        if  img1.findall(i):
            fq1.append(i)
    fq2 = [i.replace("R1","R2") for i in fq1]
    fq1 = [os.path.join(Pardic['indir'],i) for i in fq1]
    fq2 = [os.path.join(Pardic['indir'],i) for i in fq2]
    return fq1,fq2

def table(List):
    li_dic = {}
    for i in List:
        if i in li_dic:
            li_dic[i] += 1
        else:
            li_dic.update({i:1})
    return li_dic

def get_main_index(li_dic):
    big = max(li_dic.values())
    for i in li_dic:
        if li_dic[i] == big:
            index = i
            return index
        

def get_index(argv,index_num):
    Pardic = PAR(argv)
    fq1s,fq2s = find_sample(Pardic)
    fq1out = Pardic['fq1']
    fq2out = Pardic['fq2']
    stat = open(os.path.join(os.path.dirname(fq1out),os.path.basename(fq1out).split("_")[0]+".index.stat"),"w")
    i = 0
    index1 = [];index2=[]
    for fq1 in fq1s:
        fq2 = fq1.replace("R1.fastq","R2.fastq")
        if fq1.endswith(".gz"):
            gz = True
        else:
            gz = False

        L1 = []; L2 = []
        for line1,line2 in zip(xopen_fq(fq1),xopen_fq(fq2)):
            if gz:
                line1 = line1.rstrip().decode('utf-8')
                line2 = line2.rstrip().decode('utf-8')
            else:
                line1 = line1.rstrip()
                line2 = line2.rstrip()

            if line1.startswith("@"):
                i += 1
                if i == 1:
                    if "+" in line1.split(None)[1].split(":")[3]:
                        is_adapter = True
                    else:
                        is_adapter = False

                if i < index_num: continue

                if i == index_num+index_num:
                    break
                if is_adapter:
                    line1_index = line1.split(None)[1].split(":")[3].split("+")[0]
                    line2_index = line2.split(None)[1].split(":")[3].split("+")[0]
                else:
                    line1_index = line1.split(None)[1].split(":")[3]
                    line2_index = line2.split(None)[1].split(":")[3]

                index1.append(line1_index)
                index2.append(line2_index)
        table_index1 = table(index1)
        table_index2 = table(index2)
        print("{index_num}-2X{index_num} reads in fq1".format(index_num=index_num),file=stat)
        print(table_index1,file=stat)
        print("{index_num}-2X{index_num} reads in fq2".format(index_num=index_num),file=stat)
        print(table_index2,file=stat)
        
        index1 = get_main_index(table_index1)
        print("fq1 index",index1,file=stat)
        index2 = get_main_index(table_index2)
        print("fq2 index",index2,file=stat)
        if index1 != index2:
            print("fq1 main index is not same with fq2",file=sys.stderri)
            exit(1)

        stat.close()
        return is_adapter,index1
            

def check_fastq(argv,is_adapter,index):
    Pardic = PAR(argv)
    fq1s,fq2s = find_sample(Pardic)
    fq1out = Pardic['fq1']
    fq2out = Pardic['fq2']
    
    f1 = open(fq1out,"w")
    f2 = open(fq2out,"w")
    f1_drop = open(fq1out+".drop","w")
    f2_drop = open(fq2out+".drop","w")

    is_drop = False
    is_proper1 = True
    is_proper2 = True

    for fq1 in fq1s:
        fq2 = fq1.replace("R1","R2")
        if fq1.endswith(".gz"):
            gz = True
        else:
            gz = False

        L1 = []; L2 = []; DROP = {'is_drop':None,'is_proper1':None,'is_proper2':None}
        for line1,line2 in zip(xopen_fq(fq1),xopen_fq(fq2)):

            if gz:
                line1 = line1.rstrip().decode('utf-8')
                line2 = line2.rstrip().decode('utf-8')
            else:
                line1 = line1.rstrip()
                line2 = line2.rstrip()

            if line1.startswith("@"):
                if  L1:
                    if len(L1[1]) != len(L1[3]):
                        DROP['is_proper1'] = False
                    else:
                        DROP['is_proper1'] = True
                    if len(L2[1]) != len(L2[3]):
                        DROP['is_proper2'] = False
                    else:
                        DROP['is_proper2'] = True
                    
                    if is_adapter:
                        line1_index = L1[0].split(None)[1].split(":")[3].split("+")[0]
                        line2_index = L2[0].split(None)[1].split(":")[3].split("+")[0]
                    else:
                        line1_index = L1[0].split(None)[1].split(":")[3]
                        line2_index = L2[0].split(None)[1].split(":")[3]

                    if line1_index == index and line2_index == index:
                        DROP['is_drop'] = False
                    else:
                        DROP['is_drop'] = True

                    tag1 = L1[1][9:19]
                    tag2 = L2[1][9:19]
                    tag = '@' + tag1 + tag2 + '|'
                    L1[0] = tag + L1[0][1:]
                    L2[0] = tag + L2[0][1:]
                    L1_s = '\n'.join(L1);L2_s = '\n'.join(L2)
                    #print(DROP,index)
                    if DROP['is_drop'] or not DROP['is_proper1'] or not DROP['is_proper2']:
                        print(L1_s,file = f1_drop)
                        print(L2_s,file = f2_drop)
                        #print(L1_s,"drop")
                        #print(L2_s,"drop")
                    else:
                        print(L1_s,file = f1)
                        print(L2_s,file = f2)
                        #print(L1_s,"real")
                        #print(L2_s,"real")
                    L1 = []
                    L2 = []
                    #DROP = {'is_drop':False,'is_proper1':True,'is_proper2':True}
 
                L1.append(line1)
                L2.append(line2)
            else:
                L1.append(line1)
                L2.append(line2)
#last line   
     
        if len(L1[1]) != len(L1[3]):
            DROP['is_proper1'] = False
        else:
            DROP['is_proper1'] = True
        if len(L2[1]) != len(L2[3]):
            DROP['is_proper2'] = False
        else:
            DROP['is_proper2'] = True
        if is_adapter:
            line1_index = L1[0].split(None)[1].split(":")[3].split("+")[0]
            line2_index = L2[0].split(None)[1].split(":")[3].split("+")[0]
        else:
            line1_index = L1[0].split(None)[1].split(":")[3]
            line2_index = L2[0].split(None)[1].split(":")[3]

        if line1_index == index and line2_index == index:
            DROP['is_drop'] = False
        else:
            DROP['is_drop'] = True

        tag1 = L1[1][9:19]
        tag2 = L2[1][9:19]
        tag = '@' + tag1 + tag2 + '|'
        L1[0] = tag + L1[0][1:]
        L2[0] = tag + L2[0][1:]
        L1_s = '\n'.join(L1);L2_s = '\n'.join(L2)
        if DROP['is_drop'] or not DROP['is_proper1'] or not DROP['is_proper2']:
            print(L1_s,file = f1_drop)
            print(L2_s,file = f2_drop )
        else:
            print(L1_s,file = f1)
            print(L2_s,file = f2)
 



if __name__ == '__main__':
    #try:
    index_num = 2000
    is_adapter,index = get_index(sys.argv,index_num)    
    check_fastq(sys.argv,is_adapter,index)
    #except:
    #    Usage()
