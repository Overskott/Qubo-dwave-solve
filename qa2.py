import dimod
from dwave.system import LeapHybridSampler
from math import log2, floor


W = 70
N = 7
v = [35, 85, 30, 50, 70, 80, 55]
w = [12, 27, 11, 17, 20, 10, 15]

A = max(v)                                          # Lagrange multiplier (scaling constant)
M = floor(log2(W))                                  # Less than or equal to base_2 log of max_weights
K_i = [2**i for i in range(M)]                      # K_i = 2**i, i∈[0, M−1]  
K_m = [W + 1 - 2**M]                                # K_m = W −(2**m −1)
k = K_i + K_m                                       # Lucas ?


#(https://docs.ocean.dwavesys.com/en/stable/docs_dimod/reference/quadratic.html#binary-quadratic-models)

# Initialize BinaryQuadraticModel with a Bianary variable type
vartype = dimod.Vartype.BINARY
BinaryQuadraticModel =   dimod.BinaryQuadraticModel(vartype)



for i in range(N):
    X_i         = ('x' + str(i) )                                   # X_i Variable in the quadratic model. 
    X_bias      = (A * (w[i]**2) - v[i])                            # Linear bias for the X_i variable.
    BinaryQuadraticModel.set_linear(X_i, X_bias)

    for j in range(i + 1, N):
        X_i     = ('x' + str(i) )                                           
        X_j     = ('x' + str(j))
        WW_bias = (2 * A * w[i] * w[j])                             # Quadratic bias for the X_ij variables.
        
        BinaryQuadraticModel.set_quadratic(X_i, X_j, WW_bias ) 



for i in range(M + 1):
    Y_i         = ('x' + str(i))                                    # Y_i Variable in the quadratic model. 
    Y_bias      = (A * (w[i]**2) - v[i])                            # Linear bias for the Y_i variable

    BinaryQuadraticModel.set_linear(Y_i, Y_bias)

    for j in range(i + 1, M + 1):
        Y_i     = ('y' + str(i))
        Y_j     = ('y' + str(j))
        kk_bias = (2 * A * k[i] * k[j])                             # Quadratic bias for the Y_ij variable

        BinaryQuadraticModel.set_quadratic(Y_i, Y_j, kk_bias )


for i in range(N):
    for j in range(M + 1):
        X_i     = ('x' + str(i))
        Y_j     = ('y' + str(j))
        WK_bias = (-2 * A * w[i] * k[j])                            # Quadratic bias for the WK_ij variables.

        BinaryQuadraticModel.set_quadratic(X_i, Y_j,WK_bias )
 


#https://dwave-systemdocs.readthedocs.io/en/latest/reference/samplers.html

DWaveSampler = LeapHybridSampler() 
sampleset = DWaveSampler.sample(BinaryQuadraticModel)
sample = sampleset.first.sample
energy = sampleset.first.energy
selected_items = []
weight_sum = 0
value_sum = 0

for varname, value in sample.items():
    if value and varname.startswith('x'): # x*
        selected_items.append(int(varname[1:]))
selected_items = sorted(selected_items)
for i in selected_items:
    weight_sum += w[i]
    value_sum += v[i]


print("Energy:  ", energy, "  Items:  ", selected_items, "  Total weight:  ", weight_sum, "  Total value:  ", value_sum)