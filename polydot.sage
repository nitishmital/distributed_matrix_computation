############################################
## Author: Nitish Mital

## Implementation of the paper "Polynomial Codes: an Optimal Design for High-Dimensional Coded Matrix Multiplication" - Qian Yu, Mohammad Ali Maddah-Ali, Salman Avestimehr.

## This code simulates source nodes, computing cluster nodes with N workers, and the user wanting to compute a particular function of the input matrices generated on the source nodes using the computing cluster. Redundancy is introduced using polynomial codes proposed in the paper, such that the result can be obtained by the user from any K1K2 out of N workers that return their computations first.

############################################

import random

#######################################################
#######################################################

class Cluster(object):
   def __init__(self,N):
      self.N=N  # number of worker nodes in the cluster
      self.nodes=[]
      for i in range(N):
	 self.nodes.append(Worker()) # instantiating a list of worker nodes

########################################################
########################################################


class Worker(object):
     def __init__(self): 
         self.mats=[] # list of encoded matrices stored in this worker. Each element is a tuple [matrix, N,K,T]. 

#############################
     def multiply(self,ind1,ind2):  # ind1 and ind2 are the indexes of the matrices that are multiplied
         self.mats.append([self.mats[ind1][0]*self.mats[ind2][0],self.mats[ind1][1],self.mats[ind1][2]])
#############################
     def add(self,ind1,ind2):       # ind1 and ind2 are the indexes of the matrices that are added
         self.mats.append([self.mats[ind1][0] + self.mats[ind2][0],self.mats[ind1][1],self.mats[ind1][2]])


########################################################
########################################################     


class Source(object):
     def __init__(self,size):
         self.size=size
         self.mat=zero_matrix(k,self.size[0],self.size[1])
         #self.encoded_mat=matrix(k,[])

#############################
# Function generates a random matrix as the source data
#############################
     def gen_data(self): 
         self.mat=random_matrix(k,self.size[0],self.size[1])

#############################
#Function left-encodes the matrix stored in the source instance
#############################
     def matdot_rencode(self, N,K1,K):    # encodes the right-matrix (i.e. matrix B)

         s_key=[self.size[0],self.size[1]/K]  # size of the keys
         sub_ind=[i*s_key[1] for i in range(1,K)]
         self.mat.subdivide([],sub_ind)
         self.encoded_mat=zero_matrix(k,self.size[0],self.size[1]/K*N)
         sub_ind=[s_key[1]*i for i in range(1,N)]
         self.encoded_mat.subdivide([],sub_ind)
         
         for i in range(N):
             submat=zero_matrix(k,self.size[0],self.size[1]/K)
             for j in range(K):
             	submat=submat + self.mat.subdivision(0,j)*eval[i]^(K1*j)
             if i==0 :
                temp_mat = submat
             else:
                temp_mat=block_matrix([[temp_mat,submat]])
         self.encoded_mat=self.encoded_mat + temp_mat
         sub_ind=[s_key[1]*i for i in range(1,N)]
         self.encoded_mat.subdivide([],sub_ind)
 
#############################
     def matdot_lencode(self, N,K):  # encodes the left-matrix (i.e., matrix A)

         s_key=[self.size[0]/K,self.size[1]]
         sub_ind=[i*s_key[0] for i in range(1,K)]
         self.mat.subdivide(sub_ind,[])
         self.encoded_mat=zero_matrix(k,s_key[0]*N,s_key[1])
         sub_ind=[s_key[0]*i for i in range(1,N)]
         self.encoded_mat.subdivide(sub_ind,[])
         
         for i in range(N):
             submat=zero_matrix(k,self.size[0]/K,self.size[1])
             for j in range(K):
             	submat=submat + self.mat.subdivision(j,0)*eval[i]^(j)
             if i==0 :
                temp_mat = submat
             else:
                temp_mat=block_matrix([[temp_mat],[submat]])
         
         self.encoded_mat=self.encoded_mat + temp_mat
         sub_ind=[s_key[0]*i for i in range(1,N)]
         self.encoded_mat.subdivide(sub_ind,[])

###############################
#Function sends the encoded matrices to the worker nodes
###############################
     def uploadtonodes(self,clstr,K):
         N=len(clstr.nodes)
         subdiv=self.encoded_mat.subdivisions()
         if subdiv[0]==[]:
             for i in range(N):
                 clstr.nodes[i].mats.append([self.encoded_mat.subdivision(0,i),N,K])
         else:
             for i in range(N):
                 clstr.nodes[i].mats.append([self.encoded_mat.subdivision(i,0),N,K])

####################################################################
# The main program code acts as the user:
####################################################################

q=11^3
N=11     # N > K1^2 and N>K2^2
K1=3 
K2=3
k.<x>=GF(q) # Fixes a Finite Field of size q
field_elements=[i for i in k] #if i not in [0]]  # enumeration of non-zero field elements
eval=random.sample(field_elements,N)   # list of randomly picked evaluation points

###################################################
## This section generates the source matrices, encodes them at the source, and uploads the encoded shares to the computing cluster
###################################################
dim1=[60,15]
dim2=[15,60]
A=Source([dim1[0],dim1[1]])
A.gen_data()
A.matdot_lencode(N,K1)
B=Source([dim2[0],dim2[1]])
B.gen_data()
B.matdot_rencode(N,K1,K2)

comp_cluster=Cluster(N)
A.uploadtonodes(comp_cluster,K1)
B.uploadtonodes(comp_cluster,K2)


######################################################
##### This section illustrates the computation of matrix product A*B.
######################################################
for i in range(N):
	comp_cluster.nodes[i].multiply(0,1)  ## computes shares of A*B on each server

############ The user interpolates on the received matrices to obtain a polynomial whose (K-1)-th term is the desired matric product  #############

fast_servers=random.sample(range(N),K1*K2)   # the servers that return the values first, i.e., non-stragglers, simulated by randomly sampling from the N servers.
R=PolynomialRing(k,'y')
C_poly=zero_matrix(R,comp_cluster.nodes[0].mats[2][0].nrows(),comp_cluster.nodes[0].mats[2][0].ncols())
C=zero_matrix(k,dim1[0],dim2[1])
for i in range(comp_cluster.nodes[0].mats[2][0].nrows()):
    for j in range(comp_cluster.nodes[0].mats[2][0].ncols()):
        list_points=[(0,0)]*(K1*K2)  
        # at the (i,j)-th position of the resultant matrices, list_points[] forms the list of the evaluation points of the received matrices. By interpolating on the points in list_points, we will obtain the polynomial C(x) at the (i,j)-th position
        for l in range(K1*K2):
           list_points[l] = (eval[fast_servers[l]],comp_cluster.nodes[fast_servers[l]].mats[2][0][i][j]) 
        f = R.lagrange_polynomial(list_points)
        C_poly[i,j]=f

for i in range(C_poly.nrows()):
    for j in range(C_poly.ncols()):
        coeff=C_poly[i,j].coefficients(sparse=False)
        for l2 in range(K2):
            for l1 in range(K1):
               C[C_poly.nrows()*(l1)+i,C_poly.ncols()*(l2)+j]=coeff[l2*K1+l1]

           
                
verify = C - A.mat*B.mat
#print(verify)  # verify must be an all zero matrix if the result is correct.

