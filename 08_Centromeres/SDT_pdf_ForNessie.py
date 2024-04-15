import matplotlib.pyplot as plt
import pandas as pd
import argparse
import math
HuaTu = argparse.ArgumentParser(description="根据Nessie结果作图")
HuaTu.add_argument("-i","--input_nessie",required=True,help="input Nessieresult")
HuaTu.add_argument("-o","--output_png",required=True,help="out png")
HuaTu.add_argument("-Low","--highLimit",type=float,required=True,help="y_Lowest")   #Y值下限，一般是0.5  0.85 0.9
HuaTu.add_argument("-Xmax","--XmaxNume",type=float,required=True,help="W_Longest")   #X轴的最长值，是染色体长度/10000结果
HuaTu.add_argument("-Xparts","--Xmeanparts",type=int,default=15,help="x刻度拆分份数")   #X轴的间距，默认是5
HuaTu.add_argument("-Pl","--plotlongth",type=int,default=30,help="图形长度") #为了展示更细节的东西，可以将图形拉长
HuaTu.add_argument("-size","--pointsize",default=0.01,help="点的大小") #参考值8步频0.00001  1000步频1
HuaTu.add_argument("-Xid","--XaryChrname",default="Chr",help="横坐标前缀") 


args=HuaTu.parse_args()
chr = pd.read_csv(args.input_nessie, sep='\t', header=None)
Ylimit=args.highLimit
Xmax=math.ceil(args.XmaxNume)
Xparts=args.Xmeanparts
Plotlong=args.plotlongth
Outpng=args.output_png
size=float(args.pointsize)
Xid=args.XaryChrname

Xunit=Xmax//Xparts
#resule_nessie1000_E_CHY_Notitle.out
plt.figure(figsize=(Plotlong, 10))  #长宽
x=chr.iloc[:,0]             #x
y=chr.iloc[:,1]             #y 
plt.scatter(x,y,s=size)   # 画图并设置点大小 0.00001
plt.ylim(Ylimit,1)              #y值的范围 
plt.xticks(range(0,Xmax,Xunit),fontsize="24")
plt.yticks(fontsize="24")
#plt.xticks([ 0,   5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65, 70,  75,  80,  85,  90,  95, 100, 105, 110, 115, 120, 125, 130, 135,140, 145, 150])    #x的刻度
plt.xlabel(Xid+" (Mb)", fontsize="28",labelpad=30)              #x的单位
plt.ylabel("Entropy", fontsize="28",labelpad=30)                 #Y的单位

#plt.show()
plt.savefig(Outpng, dpi=600)
