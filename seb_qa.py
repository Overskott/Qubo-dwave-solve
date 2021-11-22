import dimod
from dwave.system import LeapHybridSampler
from math import log2, floor


W = 70
N = 7
c = [35, 85, 30, 50, 70, 80, 55]
w = [12, 27, 11, 17, 20, 10, 15]

A = max(c)                                          # Lagrange multiplier (scaling constant)
#M = floor(log2(W))                                  # Less than or equal to base_2 log of max_weights
B = 1


#(https://docs.ocean.dwavesys.com/en/stable/docs_dimod/reference/quadratic.html#binary-quadratic-models)

# Initialize BinaryQuadraticModel with a Bianary variable type
vartype = dimod.Vartype.BINARY
BinaryQuadraticModel =   dimod.BinaryQuadraticModel(vartype)



for a in range(N):
    X_a         = ('x' + str(a) )                                   # X_a Variable in the quadratic model. 
    X_a_bias    = ( ( A*(w[a]**2)) - (B * c[a]) )              # Linear bias for the X_a_bias variable.
    BinaryQuadraticModel.set_linear(X_a, X_a_bias)


for n in range(W):
    Y_n         = ('y' + str(n) )                                   # Y_n Variable in the quadratic model. 
    Y_n_bias    = ( A * ( (n**2) - 1) )                             # Linear bias for the Y_n_bias variable.
    BinaryQuadraticModel.set_linear(Y_n, Y_n_bias)


for a in range(N):
    for n in range(W):
        X_a     = ('x' + str(a) )
        Y_n     = ('y' + str(n) )
        W_a_bias = (-2 * A * n * w[a])                              # Quadratic bias for the X_a / Y_n variables.
        BinaryQuadraticModel.set_quadratic(X_a, Y_n, W_a_bias )



for i in range(N):
    for j in range(i+1, N):
        X_i     = ('x' + str(i) )                                           
        X_j     = ('x' + str(j) )
        WW_bias = (2 * A * w[i] * w[j])                             # Quadratic bias for the X_ij variables.
        
        BinaryQuadraticModel.set_quadratic(X_i, X_j, WW_bias ) 



for i in range(W):
    for j in range(i+1, W):
        Y_i     = ('y' + str(i) )
        Y_j     = ('y' + str(j) )
        IJ_bias = (2 * A * (1 + i*j) )                                # Quadratic bias for the Y_ij variable

        BinaryQuadraticModel.set_quadratic(Y_i, Y_j, IJ_bias )

#https://dwave-systemdocs.readthedocs.io/en/latest/reference/samplers.html



DWaveSampler = LeapHybridSampler() 
sampleset = DWaveSampler.sample(BinaryQuadraticModel)
sample = sampleset.first.sample
energy = sampleset.first.energy
selected_items = []
weight_sum = 0
value_sum = 0

for varname, value in sample.items():
    if value and varname.startswith('x'):
        selected_items.append(int(varname[1:]))
selected_items = sorted(selected_items)
for i in selected_items:
    weight_sum += w[i]
    value_sum += c[i]


print("Energy:  ", energy, "  Items:  ", selected_items, "  Total weight:  ", weight_sum, "  Total value:  ", value_sum)