from matplotlib import pyplot as plt
import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
#https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html
Breslauer= mt.DNA_NN1
Santalucia= mt.DNA_NN4 #same as benchling
Tm_Opt =Santalucia

def reverse_complement(Seq):
    Rv=''
    for bp in Seq.upper():
        if bp == 'A':
            Rv+='T'
        if bp == 'T' or bp == 'U':
            Rv+='A'
        if bp == 'C':
            Rv+='G'
        if bp == 'G':
            Rv+='C'
    assert len(Rv)==len(Seq), "There was an error reversing, Check Input Sequence"
    return Rv[::-1]


def make_anneal(Seq, direction, Tm): 
    '''defines the anneal part of a primer given it direction, a seq to amplify,and an annealing TM
    seq sequence:
    direction: F or R
    Tm: melting temp'''
    if direction == 'R':
        Seq= reverse_complement(Seq) #takes the reversed complement of Seq for reverse primer 
    
    Primer = ''
    bp=0
    Primer += Seq[bp]
    while mt.Tm_NN(Primer.upper(), nn_table=Tm_Opt, check=False) < Tm:
        bp+=1
        Primer+=Seq[bp]
    return Primer, mt.Tm_NN(Primer, nn_table=Tm_Opt, check=False)


def User_Assebmly(Parts, Tm_target, User_Overhangs):
    '''design a user assembly based on parts 
    Tm_target is the TM used as a target to design primers 
    if Tm_target == 0 specified primer sequences will be used as a reference
    if no primer sequence is input a Tm_target is necessary
    User_Overhangs is the dictionnary storing the user-defined overhangs'''
    #finds the average metling temp of the available primers (this will facilitate PCRs)
    Tms=[]
    for part in Parts:
        if len(part['Primer_up']) !=0:
            Tm = mt.Tm_NN(part['Primer_up'].upper(), nn_table=Tm_Opt, check=False)
            part['Primer_up_Tm']= Tm
            Tms.append(Tm)
        if len(part['Primer_down']) !=0:
            Tm = mt.Tm_NN(part['Primer_down'].upper(), nn_table=Tm_Opt, check=False)
            part['Primer_down_Tm']= Tm
            Tms.append(Tm)
    if Tm_target ==0: 
        Tm_target = np.mean(Tms)
        print('Unpsecified Target Tm, will use average TM of provided primers: ', Tm_target)
        if len(Tms) ==0:
                print('Error, if no primers are provided you must input at Tm_target ')
    else:
        print('Will use input Tm_target for design: ', Tm_target)
    

    # puts all part in CAP in the correct orientation:
    for part in Parts:
        part['seq'] = part['seq'].upper()
            

    #create the annealing parts of primers
    for part in Parts:
        if len(part['Primer_up'])==0:
            part['Primer_up'], part['Primer_up_Tm']= make_anneal(part['seq'], 'F', Tm_target )
        else: #capitalize primer
            part['Primer_up']= part['Primer_up'].upper()
        if len(part['Primer_down'])==0:
            part['Primer_down'], part['Primer_down_Tm']= make_anneal(part['seq'], 'R', Tm_target )
        else: #capitalize primer
            part['Primer_down']= part['Primer_down'].upper()

    #add the User sites for primers:
    
    for part in Parts:
        #define the over_hang:
        Over_hang_left = User_Overhangs[part['Primer_up_U']]
        #the Right overhang needs to be the complement:
        Over_hang_right = User_Overhangs[part['Primer_down_U']]
        Over_hang_right = reverse_complement(Over_hang_right)[0:-1]+'u'
        #stick the overhangs to the correct primers:
        
        if part['need_reverse']==0: #the seq is in the right orientation
            part['Primer_up']= Over_hang_left + part['Primer_up']
            part['Primer_down']= Over_hang_right + part['Primer_down']
            
        elif part['need_reverse']==1: #the seq needs inversion
            part['Primer_up']= Over_hang_right + part['Primer_up']
            part['Primer_down']= Over_hang_left + part['Primer_down']
    return Parts


def run_gel(Sizes, Ladder, Names=[]):
    '''Draws a gel given expected size of the bands, 
    a ladder (a list with the sizes within the ladder)
    if the names are provided, will label the lanes'''
    fig, ax = plt.subplots(figsize=(15,5))
   # ax.set_yscale('symlog')
    x_pos=5

    for lane in range(len(Sizes)):
        #place a ladder every 8 bands
        if lane%8 ==0:
            ax.hlines(y= [np.log10(l) for l in Ladder], xmin=[x_pos for i in Ladder], xmax=[x_pos + 6 for i in Ladder],
                      linewidth=[3,3,3,3,3,6,3,3,3,1], color='r')
            x_pos+=8   
        #plot band
        
        band_size= np.log10(Sizes[lane]) #log10 add some non linearity to the plot
        ax.hlines(y=band_size, xmin=x_pos, xmax=x_pos+6, linewidth=2, color='k')
        
        #annotatate band if Names provided
        if len(Names)>0:
            ax.annotate(Names[lane], (x_pos,band_size +0.4),rotation=60)
        x_pos+=8
    ax.set_ylim([-1, 1.5])
    plt.show()