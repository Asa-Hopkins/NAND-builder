from dimod.generators import and_gate, or_gate, xor_gate, halfadder_gate, fulladder_gate
from dimod.binary import Binary
from dimod import BinaryQuadraticModel
from dwave.preprocessing.lower_bounds import roof_duality
from dwave.preprocessing.composites import FixVariablesComposite
from dwave.system import LeapHybridSampler
from dwave.system import DWaveSampler, EmbeddingComposite
import dwave.inspector
import math
import pickle
import dimod

gates = [] #BQM is sum of these
out = 2 #Reserve variables 0 and 1

def nand_gate(a,b,c):
    temp = and_gate(a,b,c)
    temp.flip_variable(c)
    return temp

#Choose gate
#ac - 2acd + bbc + b + 2bcd - 2bd + d
#-2acd + ac + 2bcd - bc - 2bd + b + d
#gate = dimod.BinaryPolynomial({'a':1, 'd':1, ('b','c'):1, ('b','c','d'):-2, ('a','c'):-1, ('a','c','d'):2, ('a','d'):-2}, dimod.BINARY)
#dimod.make_quadratic(gate,1,dimod.BINARY)

def test2():
    gate = dimod.BinaryPolynomial({'a':1, 'd':1, ('b','c'):1, ('b','c','d'):-2, ('a','c'):-1, ('a','c','d'):2, ('a','d'):-2}, dimod.BINARY)
    q = dimod.make_quadratic(gate,3,dimod.BINARY)
    print(q)
    return
    for a in range(0,2):
        for b in range(0,2):
            for c in range(0,2):
                for d in range(0,2):
                    print(gate.energies({'a':a,'b':b,'c':c,'d':d}))
                    for e in [d*c,1-d*c]:
                        print(q.energies({'a':a,'b':b,'c':c,'d':d,'c*d':e,'d*c':e}))
                    print("")

def choose2(a,b,c):
    global out
    global gates
    # A gate that outputs [a,b][c]
    #gates.extend([and_gate(a,c,str(out)),and_gate(b,c,str(out+1)),or_gate(str(out),str(out+1),str(out+2))])
    #gates[-3].flip_variable(c)
    #output = BinaryQuadraticModel({str(out+1): 1.0, c: 0.0, str(out): 4.5, a: 1.0, b: 0.0}, {(c, str(out+1)): 1.5, (str(out), str(out+1)): -3.0, (str(out), c): -3.0, (a, str(out+1)): -2.0, (a, c): -1.0, (a, str(out)): 2.0, (b, c): 1.0, (b, str(out)): -2.0}, 0.0, 'BINARY')
    #output = BinaryQuadraticModel({c: 0.0, str(out+1): 1.0, str(out): 3.0, a: 1.0, b: 0.0}, {(str(out+1), c): 1.0, (str(out), c): -2.0, (str(out), str(out+1)): -2.0, (a, c): -1.0, (a, str(out+1)): -2.0, (a, str(out)): 2.0, (b, c): 1.0, (b, str(out)): -2.0}, 0.0, 'BINARY')
    output = BinaryQuadraticModel({str(out+1): 1.0, c: 0.0, str(out): 9.0, a: 1.0, b: 0.0}, {(c, str(out+1)): 3.0, (str(out), str(out+1)): -6.0, (str(out), c): -6.0, (a, str(out+1)): -2.0, (a, c): -1.0, (a, str(out)): 2.0, (b, c): 1.0, (b, str(out)): -2.0}, 0.0, 'BINARY')
    gates.append(output)
    out+=2
    if out-1 == 9:
        print(a,b,c,output)
    return str(out-1)

def choose(l1,c):
    # l is a list of inputs to be chosen from, c is a list of binary values that choose from them
    l = l1[:]
    levels = [l]
    depth = 0
    while len(levels[-1])!=1:
        levels.append([])
        length = len(levels[-2])
        for i in range(0,length-1,2):
            #print(levels[-2][i],levels[-2][i+1],c[depth])
            levels[-1].append(choose2(levels[-2][i],levels[-2][i+1],c[depth]))
        if length&1: #If there's an odd number in this level, promote the last entry to the next level
            levels[-1].append(levels[-2][-1])
        depth += 1
    return levels[-1][0]


def rchoose(l1,c):
    # l is a list of inputs to be chosen from, c is a list of binary values that chooses from them
    #print(c)
    l = l1[:]
    #print(l,c)
    levels = [l]
    depth = -1
    while len(levels[-1])!=1:
        depth += 1
        levels.append([])
        length = len(levels[-2])
        for i in range(0,length-1,2):
            #print(levels[-2][i],levels[-2][i+1],c[depth])
            levels[-1].append(choose2(levels[-2][i],levels[-2][i+1],c[depth]))
        if length&1: #If there's an odd number in this level, promote the last entry to the next level
            levels[-1].append(levels[-2][-1])
    index = c[depth]
    while levels[depth][index].isnumeric():
        #print(depth,index,levels)
        choice = levels[depth][index]
        depth -= 1
        if levels[depth][-1] == choice: #If the choice is in this level too, that means it was promoted and so we skip to the next layer 
            index = index * 2
        else:
            index = index * 2 + c[depth]
    return levels[depth][index]

def model(num_gates,input_bits, output_bits, truth_table):
    global gates
    gates = []
    #truth_table is a list of (input,output) pairs
    choice1_bits = [] #For choosing first input of each gate
    choice2_bits = [] #For choosing second input of each gate
    
    for i in range(0,num_gates):
        choice1_bits.append([])
        choice2_bits.append([])
        for j in range(0,math.ceil(math.log2(input_bits + i))):
            #Each gate can take input from all previous inputs and gates, so input_bits + i
            choice1_bits[-1].append(str(f'c1_{i}_{j}'))
            choice2_bits[-1].append(str(f'c2_{i}_{j}'))
    copy = 0
    for i in truth_table:
        choices = []
        for j in range(0,input_bits):
            choices.append(str((i[0] >> j)&1))
        print(choices, i[1])
        for j in range(0,num_gates):
            gates.append(nand_gate(choose(choices,choice1_bits[j]),choose(choices,choice2_bits[j]),f'NAND{j}_{copy}'))
            choices.append(f'NAND{j}_{copy}')
        for j in range(0,output_bits):
            bit = str((i[1] >> j)&1)
            print(bit)
            choice_out = [f'cout_{j}_{k}' for k in range(0,math.ceil(math.log2(input_bits + num_gates)))]
            output = choose(choices,choice_out)
            print(output)
            gates.append(BinaryQuadraticModel({output: 1.0, bit:1.0}, {(output,bit):-2.0}, 0.0, 'BINARY')) #Enforces equality between out{j} and the output bit
            print(gates[-1])
        copy += 1

bqm = 0
sampleset = 0
sample = 0
def submit():
    global bqm
    global sampleset
    global sample
    global gates
    input_bits = 2
    output_bits = 2
    num_gates = 5
    model(num_gates,input_bits,output_bits,[(0,0),(1,1),(2,1),(3,2)]) #Half adder circuit, best known solution is 5 gates
    bqm = sum(gates)
    bqm.fix_variables({'0':0,'1':1}) #Fix '0' and '1'
    #return
    fix = roof_duality(bqm, strict = False)[1] #This lets around 10% of values be fixed
    bqm.fix_variables(fix)
    #sampler = DWaveSampler()
    #sampler_embedded = EmbeddingComposite(sampler)
    #sample = sampler_embedded.sample(bqm, num_reads=5000,label='NAND circuit')
    
    #sampler = DWaveSampler()
    #print(fix)
    #f = open('sample.set','rb')
    #sampleset = pickle.load(f)
    sampler = LeapHybridSampler()
    #sample = sampler.sample(bqm)
    sample = sampler.sample(bqm,label='NAND circuit')
    sampleset = list(sample.data())[0][0]
    #return
    f = open('sample.set','wb')
    try:
        pickle.dump(sample,f)
    except:
        pickle.dump(sampleset,f)
    f.close()
    d = sampleset#list(sampleset.data())[0][0] #dict of variable values
    inputs = []
    for i in range(0,input_bits):
        inputs.append(f'input{i}')
    for i in range(0,num_gates):
        c1 = []
        c2 = []
        for j in range(0,math.ceil(math.log2(len(inputs)))):
            c1.append(int(d[f'c1_{i}_{j}']))
            c2.append(int(d[f'c2_{i}_{j}']))
        
        print(f'NAND{i} = {rchoose(inputs,c1)} NAND {rchoose(inputs,c2)}')
        inputs.append(f'NAND{i}')
    for i in range(0,output_bits):
        cout = []
        for j in range(0,math.ceil(math.log2(len(inputs)))):
            cout.append(int(d[f'cout_{i}_{j}']))
        print(f'OUT{i} = {rchoose(inputs,cout)}')

def test():
    for i in gates:
        temp = {}
        for j in i.variables:
            if j in solution:
                temp[j] = solution[j]
        #print(i)
        i.fix_variables(temp)
        #print(i, roof_duality(i)[1])
        #input("enter to continue")
        solution.update(roof_duality(i)[1])
##    for i in gates:
##        print(i)
##        if i.energies(solution)!=0:
##            print(i.energies(solution))

# Hard-coded solution for reference
##NAND0 = input0 NAND input1
##NAND1 = input0 NAND NAND0
##NAND2 = input1 NAND NAND0
##NAND3 = NAND0 NAND NAND0
##NAND4 = NAND1 NAND NAND2
##Out0 = NAND3
##Out1 = NAND4
solution = {'0':0, '1':1,
'c1_0_0':0, 'c2_0_0':1,
'c1_1_0':0, 'c1_1_1':0, 'c2_1_0':0, 'c2_1_1':1,
'c1_2_0':1, 'c1_2_1':0, 'c2_2_0':0, 'c2_2_1':1,
'c1_3_0':0, 'c1_3_1':1, 'c1_3_2':0, 'c2_3_0':0, 'c2_3_1':1, 'c2_3_2':0,
'c1_4_0':0, 'c1_4_1':0, 'c1_4_2':1, 'c2_4_0':1, 'c2_4_1':1, 'c2_4_2':0,
'cout_0_0':0, 'cout_0_1':1, 'cout_0_2':1,
 'cout_1_0':1, 'cout_1_1':0, 'cout_1_2':1}

##f = open('sample.set','rb')
##solution = pickle.load(f)
##f.close()

