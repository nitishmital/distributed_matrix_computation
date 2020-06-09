######## Minimum Bandwidth Regenerating code ###############
## Author : Nitish Mital

## Implementation of my paper "Practical Functional Regenerating Codes for Broadcast Repair of Multiple Nodes", ISIT Paris, 2019. The code implementation resembles centralized repair more closely because all repair packets are collected in one matrix (Yx, Yy). Theoretically, centralized repair and broadcast repair are equivalent.

######## Generates a random message, encodes and stores it in n nodes using the code proposed in the paper.
######## Contains a function repair() that takes as arguments lists of helper node indices and newcomer node indices, and implements the repair procedure proposed in the paper.

## This code works only if k and d are divisible by r. An update for other values will be posted later.
################################################################################################

import random

############################
########## Function that repairs the newcomers using the scheme in the paper ###########################

def repair(helpers, newcmrs):
    tx_pkt_ind = [[0 for i in range(len(newcmrs))] for j in range(len(helpers))] # indices of the transmitted packets from each helper node
    for i in range(len(helpers)):
       ind = [t for t in range(n) if t not in [helpers[i]]]  # labels corresponding to all other nodes except the self for the stored packets in each node 
       for j in range(len(newcmrs)):
          tx_pkt_ind[i][j] = ind.index(newcmrs[j])

    Yx = matrix([[nodes[helpers[i]][tx_pkt_ind[i][j],0] for i in range((h-1)*r,h*r)] for j in range(r) for h in range(1,d/r+1)])  ## repair matrix for the evaluation points
    Yy = matrix([[nodes[helpers[i]][tx_pkt_ind[i][j],1] for i in range((h-1)*r,h*r)] for j in range(r) for h in range(1,d/r+1)]) ## repair matrix for the coded message symbols 

    Yx_repair = Yx*C_gen_repair[:,r:]  # Matrix of repair packets
    Yy_repair = Yy*C_gen_repair[:,r:]
    c=0
    for i in newcmrs:
        nodes[i][:d,0] = Yx_repair[:,c]
        nodes[i][:d,1] = Yy_repair[:,c]
        nodes[i][d:,0] = C_gen_storage[:,d:].transpose() * nodes[i][:d,0]
        nodes[i][d:,1] = C_gen_storage[:,d:].transpose() * nodes[i][:d,1]
        c = c+1 
###############################################
######## Function to verify the linear independence of a list of points in extension field ##########################
######## Returns True is they are linearly independent, False if not #########################

def verify_independence(vpoints):
    lnth=len(vpoints)
    v=[[base_k(0)]*l for i in range(lnth)]
    for i in range(lnth):
        v[i] = vector(vpoints[i])
    return((base_k^l).linear_dependence([v[j] for j in range(lnth)]) == [])
################################################

n=8
k=4
d=4
r=2
subpacketization = int(k/2*(2*d - k + r))
q=11
base_k=GF(q)
l=d^2  # size of extension chosen so that there are enough linearly independent points in the finite extension field
k.<x>=GF(q^l) 

Frob=k.frobenius_endomorphism();
S.<y>=k['y',Frob]
msg=[0]*subpacketization
for i in range(subpacketization):
   msg[i] = k.random_element()   # list of randomly picked points which act as the message symbols

####### Construct linearized polynomial ################
fy=0
for i in range(subpacketization):
    fy=fy+msg[i]*y^i

####### snippet to choose l=d^2 linearly independent evaluation points ##########
points=[[0,0] for i in range(l)]
v=[[base_k(0)]*l for i in range(l)]

while(True):
    for i in range(l):
        points[i][0] = k.random_element()   # list of randomly picked points which act as the evaluation points
        points[i][1] = fy(points[i][0])
        v[i] = vector(points[i][0])
    print("initial evaluation points are linearly independent - ", (base_k^l).linear_dependence([v[j] for j in range(l)]) == [])
    if ((base_k^l).linear_dependence([v[j] for j in range(l)]) == []):
        break
##############################################


base_field_elements = [i for i in base_k if i not in [0]]

eval_pts = random.sample(base_field_elements, 2*r)
C_repair = codes.GeneralizedReedSolomonCode(eval_pts,r)
C_gen_repair = C_repair.systematic_generator_matrix()

eval_pts = random.sample(base_field_elements, n-1)
C_storage = codes.GeneralizedReedSolomonCode(eval_pts,d)
C_gen_storage = C_storage.systematic_generator_matrix()
##############################################


########## Store coded packets/symbols in the first d nodes ############
nodes = [None]*n
for i in range(n):
    ####### initializing n-1 placeholders in each node. Each row is of type (x,y), where x is the evaluation point, and y is the coded packet (evaluation of fy on x)
    nodes[i] = zero_matrix(k,n-1,2) 

for i in range(d):
   nodes[i][:d,:] = matrix(points[i*d:(i+1)*d])   # d packets in node i
   nodes[i][d:,:] = C_gen_storage[:,d:].transpose() * nodes[i][:d,:]  # storing the parity bits

##############################################

########## Fill in the remaining n-d nodes r nodes at a time by using the repair scheme #############
for i in range((n-d)/r):
    newcmrs = range(d+r*i,d+r*(i+1))  # node indices of r newcomers
    helpers = range(d)  # nodes indices of d helper nodes
    repair(helpers,newcmrs)



