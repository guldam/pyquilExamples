"""
implementation of the swap test and dotproduct
Refer <Understanding Quantum Support Vector Machine> written by Petter Witteck
https://peterwittek.com/understanding-quantum-svms.html
"""
import numpy as np
from pyquil.quil import Program
from pyquil.gates import H,CSWAP,X,MEASURE,CNOT
from numpy import linalg as LA
from pyquil.paulis import ID, sX, sY, sZ,term_with_coeff,PauliSum,PauliTerm,exponential_map,exponentiate_commuting_pauli_sum
from pyquil.api import WavefunctionSimulator
import math
from pyquil import get_qc

qc = get_qc('9q-qvm')
wf_sim = WavefunctionSimulator()

def estimateZ(A,B,N,t):
    absA = LA.norm(np.asarray(A),1)
    absB = LA.norm(np.asarray(B),1) 

    A = 0.5*(absA+absB)
    B = 0.5*(absA-absB)
    I = PauliTerm("I",0,coefficient=A)
    Z = PauliTerm("Z",0,coefficient=B)
    H = (I+Z)*PauliTerm("X",1)
    exponential = exponentiate_commuting_pauli_sum(H)

    value=0
    for i in range(N):
        phi = Program(X(0),H(0))
        ro = phi.declare('ro', 'BIT', 1)
        phi += exponential(t)
        phi += MEASURE(1, ro[0])

        if(wf_sim.wavefunction(phi).amplitudes[0]==0):
            value +=1
    return math.sqrt(2*(value/N))/(t)



def swap_test_program(ancilla,registerA,registerB):
    swapTest = Program()
    swapTest += H(ancilla)
    for a in reversed(registerA):
        for b in registerB:
            swapTest += CSWAP(ancilla,a,b)
    swapTest = Program()

def run_swap_test(programA,programB,num_measurement,quantum_resources,classical_memory,ancilla):
    registerA = list(programA)
    registerB = list(programB)
    
    swap_test_program = Program()
    swap_test_program += programA + programB
    swap_test_program += swap_test_program(ancilla,registerA,registerB)
    ro = swap_test_program.declare('ro', 'BIT', 1)
    swap_test_program.measure(ancilla,[classical_memory])
    
    results = quantum_resources.run(swap_test_program,classical_memory,trials=num_measurement)
    probabiliy_of_one = np.mean(results)
 
    return 1-2*probabiliy_of_one




def dotProduct(A,B):
    """
    calculate dot product between two boolean vectors of same length 
    """
    absA = LA.norm(np.asarray(A),1)
    absB = LA.norm(np.asarray(B),1) 

    a = 0.5*(absA+absB)
    b = 0.5*(absA-absB)
    I = PauliTerm("I",0,coefficient=a)
    Z = PauliTerm("Z",0,coefficient=b)
    Hamil = (I+Z)*PauliTerm("X",1)
    exponential = exponentiate_commuting_pauli_sum(Hamil)

    zeros = 0
    ones = 0
    while(True):
        #print('while..')
        phi = Program()
        phi.inst(X(0))
        phi.inst(H(0))
        phi.inst(exponential(0.1))
        ro = phi.declare('ro', 'BIT', 1)
        phi += MEASURE(1, ro[0])
        if(wf_sim.wavefunction(phi).amplitudes[0]==0):
            psi = Program(H(1))
            for index in range(len(B)):
                if B[index]==1:
                    psi = psi +CNOT(1,index+2)
            psi = psi+X(1)
            for index in range(len(A)):
                if A[index]==1:
                    psi = psi+CNOT(1,index+2)
            psi = psi+X(0)

            ancilla = Program(H(2+len(A)))
            for index in range((len(B))):
                ancilla += CSWAP(2+len(A),index,index+1)  
            ancilla += H(2+len(A))
    
            swaptest = phi + psi + ancilla

            ro2 = swaptest.declare('ro2', 'BIT', 1)
            swaptest += MEASURE(2+len(A),ro2[0])
            for i in range(2**(3+len(A))):
                if(wf_sim.wavefunction(swaptest).amplitudes[i]!=0):
                    if(i<(2**(2+len(A)))):
                        zeros+=1
                    else:
                        ones+=1
                    break
            print('measure success')
            print('zeros+ones = ',zeros+ones)
            if(zeros+ones>1000):
                print('zeors and ondes are',zeros,ones)
                return (zeros-ones)/(zeros+ones)
                break
            else:
               continue
        else:
            continue

dotProduct([1,1,1,1,1],[0,1,0,0,0])
### x^T y = 0.5z(1-Z(P0 - P1))

    


   


