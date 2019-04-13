from Bio import Entrez, SeqIO

import os
import sys
import multiprocessing
import datetime

import re
import math
import random

from email_validator import validate_email, EmailNotValidError

sys.setrecursionlimit(10000)    #再起深度の上限

import collections
Result = collections.namedtuple("Result", "i1, i2, alignment, distance")
Matrices = collections.namedtuple("Matrices", "alignment, distance")
MS = collections.namedtuple("MS", "score, direction")

## 一致を+3pt、不一致を-1pt、ギャップを-2ptとする ##
SC = lambda nucleotide1, nucleotide2: 3 if nucleotide1==nucleotide2 else -1
GP = -2

def main():
    Entrez.email = input_email()
    
    show_help()
    rrn16S=[]
    descriptions_list=[]
    matrices = Matrices([],[])
    clustering=None
    while True:
        term_input_list=[]
        while True:
            print('')
            term_input = input("input: ")
            if re.search(r'(\w+\([\d\s,]*\))', term_input):
                end = handle_command(term_input, descriptions_list, matrices, clustering)
                if end: break
                continue
            term_input_list.append(term_input)
        
        for term_input in term_input_list:
            gbids=[]
            gbids = search_gbids(term_input)

            random.shuffle(gbids)
            
            if len(gbids)==0: continue
            
            for gbid in gbids:
                record = fetch_record(gbid)
                seq = extract_16SrRNA(record)
                if seq: break
            
            if seq == None:
                print('')
                print(term_input+" doesn't heve 16S-rRNA data.")
                continue
            
            description = record.description.replace(', complete genome','')
            descriptions_list.append(description)
            rrn16S.append(seq)

            print('')
            print(str(len(rrn16S)-1)+'. '+description)
            print("ID: "+gbid)
            print("16S rRNA sequence: ")
            print(seq)
        
        if len(rrn16S)>1:
            matrices, clustering = add_alignment(matrices, rrn16S, len(matrices.alignment)+1)
            show_clustering(matrices.distance, clustering)

#################
# Handle inputs #
#################

def input_email():
    while True:
        print('')
        email = input("Input your email: ")
        try:
            v = validate_email(email)
            return v["email"]
        except EmailNotValidError as e:
            print(str(e))
            continue
    
def handle_command(term_input, descriptions_list, matrices, clustering):
    command=re.findall(r'(\w+)', term_input)[0]
    arguments=re.findall(r'(\d+)', (re.findall(r'(\([\d\s,]*\))', term_input)[0]))
    try:
        if   command=='nametable':
            show_nametable(descriptions_list)
        elif command=='sequence':
            show_sequence(rrn16S, int(arguments[0]))
        elif command=='alignment':
            show_alignment(matrices, int(arguments[0]), int(arguments[1]))
        elif command=='cluster':
            show_clustering(matrices.distance, clustering)
        elif command=='end':
            return True
        elif command=='save':
            save_text(descriptions_list, matrices.distance)
        elif command=='exit':
            sys.exit()
        else:
            show_help()
        return False
    except:
        show_help()
        return False
    
############################
# Fetch 16S-rRNA sequences #
############################
    
def search_gbids(term_input):
    term_ = correct_spell(term_input)
    
    handle = Entrez.esearch(db="nucleotide", term=term_+"[orgn] AND complete genome[title]", idtype="acc", retMax='1000')
    record = Entrez.read(handle)
    gbids = record["IdList"]
    
    if len(gbids)==0:
        print('')
        print("\""+term_+"\" does not exist")
    return gbids
    
def correct_spell(term_):
    handle = Entrez.espell(term=term_)
    record = Entrez.read(handle)
    
    if record["CorrectedQuery"]!='':
        term_=record["CorrectedQuery"]
        print('')
        print("Showing results for \""+term_+"\"")
        
    return term_
        
def fetch_record(id_):
    handle = Entrez.efetch(db="nucleotide", id=id_, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record

def extract_16SrRNA(record):
    for feature in record.features:
        if feature.type != 'rRNA': continue
        if feature.qualifiers['product']==['16S ribosomal RNA']:
            seq=feature.location.extract(record.seq)
            return seq
    return None

#####################
# Execute alignment #
#####################

def add_alignment(matrices, seqs, len_seqs_b):
    len_seqs = len(seqs)
    canceled = False
    jobs = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    concurrency=multiprocessing.cpu_count()
    create_processes(seqs, jobs, results, concurrency)
    add_jobs(len_seqs, len_seqs_b, jobs)
    
    try:
        jobs.join()
    except KeyboardInterrupt:
        canceled = True
    
    for i in range(len_seqs_b-1):
        for _ in range(len_seqs - len_seqs_b):
            matrices.alignment[i].append(None)
            matrices.distance[i].append(None)
    for _ in range(len_seqs - len_seqs_b):
        matrices.alignment.append([None for _ in range(len_seqs-1)])
        matrices.distance.append([None for _ in range(len_seqs-1)])
    
    while not results.empty():
        result=results.get_nowait()
        matrices.alignment[result.i1][result.i2-1]  = result.alignment
        matrices.distance[result.i1][result.i2-1]   = result.distance
    
    clustering = neighbor_joining(matrices.distance)
    
    return matrices, clustering

def create_processes(seqs, jobs, results, concurrency):
    for _ in range(concurrency):
        process = multiprocessing.Process(target=worker, args=(seqs, jobs, results))
        process.daemon = True
        process.start()    
    
def worker(seqs, jobs, results):
    while True:
        try:
            si1, si2 = jobs.get()
            try:
                alignment, distance = align(seqs[si1], seqs[si2])
                result = Result(si1, si2, alignment, distance)
                results.put(result)
            except:
                 pass   
        finally:
            jobs.task_done()

def add_jobs(len_seqs, len_seqs_b, jobs):
    for i2 in range(len_seqs - len_seqs_b):
        for i1 in range(len_seqs_b + i2):
            jobs.put((i1, len_seqs_b + i2))

def align(seq1, seq2):
    memo = [[None for _ in range(len(seq2)+1)] for __ in range(len(seq1)+1)]
    memo[0][0] = MS(0, 'S')
    for si1 in range(len(seq1)): memo[si1+1][0] = MS((si1+1)*GP, 'H')
    for si2 in range(len(seq2)): memo[0][si2+1] = MS((si2+1)*GP, 'V')
    
    alignment=['' for _ in range(3)]
    si1=len(seq1)
    si2=len(seq2)
    length_of_alignment=0
    noMatch=0
    while si1!=0 or si2!=0:
        length_of_alignment+=1
        if   maxScore(si1, si2, seq1, seq2, memo).direction=='D':
            alignment[0] += seq1[si1-1]
            alignment[1] += '|' if seq1[si1-1]==seq2[si2-1] else ' '
            alignment[2] += seq2[si2-1]
            si1-=1
            si2-=1
        elif maxScore(si1, si2, seq1, seq2, memo).direction=='V':
            alignment[0] += '-'
            alignment[1] += ' '
            alignment[2] += seq2[si2-1]
            si2-=1
            noMatch+=1
        elif maxScore(si1, si2, seq1, seq2, memo).direction=='H':
            alignment[0] += seq1[si1-1]
            alignment[1] += ' '
            alignment[2] += '-'
            si1-=1
            noMatch+=1
    for i in range(3):
        alignment[i]=alignment[i][::-1]
    
    f = noMatch/length_of_alignment
    distance = abs(3/4*math.log(1-4*f/3))
    
    return alignment, distance

def maxScore(mi1, mi2, seq1, seq2, memo):
    if memo[mi1][mi2]: return memo[mi1][mi2]
    
    scoreList=[]
    scoreList.append(maxScore(mi1-1, mi2-1, seq1, seq2, memo).score + SC(seq1[mi1-1], seq2[mi2-1]))
    scoreList.append(maxScore(mi1  , mi2-1, seq1, seq2, memo).score + GP)
    scoreList.append(maxScore(mi1-1, mi2  , seq1, seq2, memo).score + GP)
    
    direction = scoreList.index(max(scoreList))
    if   direction == 0: direction = 'D'
    elif direction == 1: direction = 'V'
    elif direction == 2: direction = 'H'
    
    memo[mi1][mi2] = MS(max(scoreList), direction)
    return memo[mi1][mi2]

################
# Make cluster #
################

def neighbor_joining(distance_matrix):
    dm = [dr[:] for dr in distance_matrix]
    node_list=[i for i in range(len(distance_matrix)+1)]
    while len(node_list)>2:
        r=[]
        for node1 in range(len(node_list)):
            r_i = 0
            for node2 in range(len(node_list)):
                if   node1<node2: r_i += dm[node1][node2-1]
                elif node1>node2: r_i += dm[node2][node1-1]
            r_i /= len(node_list)-2
            r.append(r_i)
        min_drr=None
        min_drr_i1=None
        min_drr_i2=None
        for di1, dr in enumerate(dm):
            node1 = di1
            for di2, d in enumerate(dr[di1:]):
                node2 = di1+di2+1
                drr = d - r[node1]- r[node2]
                if min_drr==None or drr < min_drr:
                    min_drr = drr
                    min_drr_i1 = di1
                    min_drr_i2 = di1+di2+1
        
        node_list.append((node_list[min_drr_i1], node_list[min_drr_i2]))
        node_list.pop(min_drr_i2)
        node_list.pop(min_drr_i1)
        
        dm.append([None for _ in range(len(dm))])
        for dri, dr in enumerate(dm):
            if dri in {min_drr_i1, min_drr_i2}: continue
            d_mi1 = dr[min_drr_i1-1] if dri<min_drr_i1 else dm[min_drr_i1][dri-1]
            d_mi2 = dr[min_drr_i2-1] if dri<min_drr_i2 else dm[min_drr_i2][dri-1]
            d_mk = (d_mi1 + d_mi2 - dm[min_drr_i1][min_drr_i2-1])/2
            dr.append(d_mk)
        dm.pop(min_drr_i2)
        dm.pop(min_drr_i1)
        for dr in dm:
            dr.pop(min_drr_i2-1)
            dr.pop(min_drr_i1-1) if min_drr_i1>0 else dr.pop(min_drr_i1-0)
        
    clustering = (node_list[0],node_list[1])
    return clustering

#######################
# Output informations #
#######################

def show_clustering(matrix, clustering):
    print('')
    print('{:>3}'.format('') + '|' + 
              '|'.join('{:>6}'.format(str(mi2+1)) for mi2 in range(len(matrix))))
    for mi1 in range(len(matrix)):
        print('{:>3}'.format(str(mi1)) + '|' + 
              '|'.join('{:>6}'.format(str('{:.3f}'.format(matrix[mi1][mi2])) if matrix[mi1][mi2]!=None else '') 
                       for mi2 in range(len(matrix))))
    
    print('')
    print(clustering)

def show_nametable(descriptions_list):
    print('')
    for di, description in enumerate(descriptions_list):
        print('{:>3}'.format(str(di)) + '. ' + description)

def show_sequence(seqs, si):
    print('')
    print(seqs[si])

def show_alignment(matrices, si1, si2):
    print('')
    WIDTH=os.get_terminal_size().columns
    alignment = matrices.alignment[si1][si2-1]
    sai=0
    while sai<len(alignment[0]):
        length = min(len(alignment[0])-sai, WIDTH)
        for i in range(3):
            print(alignment[i][sai:sai+length])
        print('')
        sai+=length
    
    similarity = 1 - math.tanh(matrices.distance[si1][si2-1]*4/3)
    print('Similarity: '+str(similarity))
    
def show_help():
    print("""
Usage:
    1. input: <academic name>
    2. input: <command(arg)>
    
Commands:
    nametable()             Show nametable
    sequence(int)           Show DNA-sequence
    alignment(int, int)     Show alignment result
    cluster()               Show clustering result
    end()                   Stop inputting and start analizing
    save()                  Save to text file(same directory)
    exit()                  Exit this application""")

#####################
# Save to text file #
#####################
    
def save_text(descriptions_list, matrix):
    text=[]
    for di, description in enumerate(descriptions_list):
        text.append('{:>3}'.format(str(di)) + '. ' + description)
    text.append('')
    
    for mi1 in range(len(matrix)):
        for mi2 in range(len(matrix)):
            if matrix[mi1][mi2]==None: continue
            similarity = 1 - math.tanh(matrix[mi1][mi2]*4/3)
            text.append(str(similarity)+', '+str(mi1)+'-'+str(mi2+1))
            
    filename=datetime.datetime.now().strftime('%y%m%d%H%M%S')+'.txt'
    with open(os.getcwd()+'/'+filename,'w') as f:
        f.write('\n'.join(text))
    
if __name__ == "__main__":
    main()