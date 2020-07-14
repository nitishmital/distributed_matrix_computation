############################################
## Author: Nitish Mital
## Implementation of the paper "Secure Distributed Matrix Computation with Discrete Fourier Transform" - Nitish Mital, Cong Ling, Deniz Gunduz. 
## Open access - https://arxiv.org/abs/2007.03972

## This code simulates distributed source nodes, a computing cluster with N servers, and the user wanting to compute a particular function of the input matrices generated on the source nodes using the computing cluster. 

## Class Cluster emulates a cluster of N computing servers. It has cluster level functions like inter-server communication.

## Class Worker emulates a computing server. It has server level functions like matrix operations, encoding and decoding.

## Class Source emulates the data sources. It has source node level functions like generation of data, encoding, and uploading the data to the computing cluster.

## The fixed part of the code is the classes and functions, while the main section consists of a demonstration of how to use the classes to perform desired matrix computations and operations.

## This code implements matrix multiplication, addition, conversion of shares. The matrix exponentiation and matrix inversion algorithms from the paper are not included in this code. Also, unlike the paper, we employ a naive discrete Fourier transform algorithm, that is, simple evaluation of a polynomial on the roots of unity, instead of applying the fast Fourier transform algorithm.

############################################

#######################################################
#######################################################

class Cluster(object):
   def __init__(self,N):
      self.N=N  # number of worker nodes in the cluster
      self.nodes=[]
      for i in range(N):
	 self.nodes.append(Worker()) # instantiating a list of worker nodes
#############################  
   '''
   def left_to_right_share(self,index_mat,K2,T2):
      for i in range(self.N):
         self.nodes[i].right_encode(index_mat,self.N,K2,T2)
      self.communication_phase(index_mat+1)
      for i in range(self.N):
         self.nodes[i].left_interpolate(index_mat+2, self.N,1)
   ''' 
#############################
   '''
   def right_to_left_share(self,index_mat,K2,T2):
      for i in range(self.N):
         self.nodes[i].left_encode(index_mat,self.N,K2,T2)
      self.communication_phase(index_mat+1)
      for i in range(self.N):
         self.nodes[i].right_interpolate(index_mat+2,self.N,1)
   '''
#############################
   def communication_phase(self,index_mat):
      subdiv = self.nodes[0].mats[index_mat][0].subdivisions()
      if subdiv[1]==[]:
         for i in range(N):
            recv_i=self.nodes[0].mats[index_mat][0].subdivision(i,0)
            for j in range(1,N):
                  recv_i = block_matrix([[recv_i, self.nodes[j].mats[index_mat][0].subdivision(i,0)]])
            self.nodes[i].mats.append([recv_i,self.nodes[0].mats[index_mat][1],self.nodes[0].mats[index_mat][2],self.nodes[0].mats[index_mat][3]])

      else:
         for i in range(N):
            recv_i=self.nodes[0].mats[index_mat][0].subdivision(0,i)
            for j in range(1,N):
                 recv_i = block_matrix([[recv_i], [self.nodes[j].mats[index_mat][0].subdivision(0,i)]])
            self.nodes[i].mats.append([recv_i,self.nodes[0].mats[index_mat][1],self.nodes[0].mats[index_mat][2],self.nodes[0].mats[index_mat][3]])
              

########################################################
########################################################





class Worker(object):
     def __init__(self): 
         self.mats=[] # list of encoded matrices stored in this worker. Each element is a tuple [matrix, N,K,T]. 

#############################
     def multiply(self,ind1,ind2):  # ind1 and ind2 are the indexes of the matrices that are multiplied
         self.mats.append([self.mats[ind1][0]*self.mats[ind2][0],self.mats[ind1][1],1,self.mats[ind1][3]])
#############################
     def add(self,ind1,ind2):       # ind1 and ind2 are the indexes of the matrices that are added
         self.mats.append([self.mats[ind1][0] + self.mats[ind2][0],self.mats[ind1][1],self.mats[ind1][2],self.mats[ind1][3]])
#############################
     def gen_lkeys(self,s_key, T): 
         R=random_matrix(k,s_key[0],T*s_key[1])
         sub_ind=[s_key[1]*i for i in range(1,T)]
         R.subdivide([],sub_ind)
         return(R)
#############################
     def gen_rkeys(self,s_key,T):
         R=random_matrix(k,T*s_key[0],s_key[1])
         sub_ind=[s_key[0]*i for i in range(1,T)]
         R.subdivide(sub_ind,[])
         return(R)
#############################
     def left_encode(self, index_mat, N, K, T):    # outputs encoded matrix with subdivision markers
         #List of N^th roots of unity
         roots=[k.zeta(N)^i for i in range(N)]
         size=[self.mats[index_mat][0].nrows(),self.mats[index_mat][0].ncols()]
         s_key=[size[0],size[1]/K]  # size of the keys
         R=self.gen_lkeys(s_key, T)
         sub_ind=[i*s_key[1] for i in range(1,K)]
         local_mat=copy(self.mats[index_mat][0])
         local_mat.subdivide([],sub_ind)
         encoded_mat=zero_matrix(k,size[0],size[1]/K*N)
         sub_ind=[i*s_key[1] for i in range(1,N)]
         encoded_mat.subdivide([],sub_ind)
         
         for i in range(N):
             submat=zero_matrix(k,size[0],size[1]/K)
             for j in range(K):
             	submat=submat + local_mat.subdivision(0,j)*roots[i]^j
             if i==0 :
                temp_mat = submat
             else:
                temp_mat=block_matrix([[temp_mat,submat]])
         encoded_mat=encoded_mat + temp_mat
         encoded_mat.subdivide([],sub_ind)
         for i in range(N):
             submat=zero_matrix(k,size[0],size[1]/K)
             for j in range(K,K+T):
             	submat=submat + R.subdivision(0,j-K)*roots[i]^j
             if i==0 :
                temp_mat = submat
             else:
                temp_mat=block_matrix([[temp_mat,submat]])
         
         encoded_mat=encoded_mat + temp_mat
         encoded_mat.subdivide([],sub_ind)
         self.mats.append([encoded_mat,N,K,T])

#############################
     def right_encode(self, index_mat, N, K, T):    # outputs encoded matrix with subdivision markers
         #List of N^th roots of unity
         roots=[k.zeta(N)^i for i in range(N)]
         size=[self.mats[index_mat][0].nrows(),self.mats[index_mat][0].ncols()]
         s_key=[size[0]/K,size[1]]  # size of the keys
         R=self.gen_rkeys(s_key, T)
         sub_ind=[i*s_key[0] for i in range(1,K)]
         local_mat=copy(self.mats[index_mat][0])
         local_mat.subdivide(sub_ind,[])
         encoded_mat=zero_matrix(k,size[0]/K*N,size[1])
         sub_ind=[i*s_key[0] for i in range(1,N)]
         encoded_mat.subdivide(sub_ind,[])
         
         for i in range(N):
             submat=zero_matrix(k,size[0]/K,size[1])
             for j in range(K):
             	submat=submat + local_mat.subdivision(j,0)*roots[i]^(-j)
             if i==0 :
                temp_mat = submat
             else:
                temp_mat=block_matrix([[temp_mat],[submat]])

         encoded_mat=encoded_mat + temp_mat
         encoded_mat.subdivide(sub_ind,[])

         for i in range(N):
             submat=zero_matrix(k,size[0]/K,size[1])
             for j in range(K+T,K+2*T):
             	submat=submat + R.subdivision(j-K-T,0)*roots[i]^(-j)
             if i==0 :
                temp_mat = submat
             else:
                temp_mat=block_matrix([[temp_mat],[submat]])
         
         encoded_mat=encoded_mat + temp_mat
         encoded_mat.subdivide(sub_ind,[])
         self.mats.append([encoded_mat,N,K,T])

############################
## Computes the first K1 coefficients of inverse FFT of right shares. 
## K1 : the paritioning of the original right-encoding
## This function is used by the function right_to_left_share() in Cluster class.
############################ 
     def right_interpolate(self,index_mat,N,K1):  
       #List of N^th roots of unity
       roots=[k.zeta(N)^i for i in range(N)]
       skip=self.mats[index_mat][0].nrows()/N
       sub_ind=[i*skip for i in range(1,N)]
       local_mat=copy(self.mats[index_mat][0])
       local_mat.subdivide(sub_ind,[])
       interpolated_mat=[]
       for i in range(K1):
          temp_mat=zero_matrix(k, local_mat.subdivision(0,0).nrows(),local_mat.subdivision(0,0).ncols())
          for j in range(N):
             temp_mat=temp_mat + local_mat.subdivision(j,0)*roots[j]^(i)/N
          interpolated_mat.append(temp_mat)
       length_list=len(interpolated_mat)
       interpolated_matrix=interpolated_mat[0]
       for i in range(1,length_list):
          interpolated_matrix=block_matrix([[interpolated_matrix],[interpolated_mat[i]]])
       #print(interpolated_matrix.nrows(),interpolated_matrix.ncols())
       self.mats.append([interpolated_matrix,self.mats[index_mat][1],self.mats[index_mat][2],self.mats[index_mat][3]])

#############################
## Computes the first K1 coefficients of inverse FFT of left shares. 
## K1 : the paritioning of the original left-encoding
## This function is used by the function left_to_right_share() in Cluster class.
#############################
     def left_interpolate(self,index_mat,N,K1):  
       #List of N^th roots of unity
       roots=[k.zeta(N)^i for i in range(N)]
       skip=self.mats[index_mat][0].ncols()/N
       sub_ind=[i*skip for i in range(1,N)]
       local_mat=copy(self.mats[index_mat][0])
       local_mat.subdivide([],sub_ind)
       interpolated_mat=[]
       for i in range(K1):
          temp_mat=zero_matrix(k, local_mat.subdivision(0,0).nrows(),local_mat.subdivision(0,0).ncols())
          for j in range(N):
             temp_mat=temp_mat + local_mat.subdivision(0,j)*roots[j]^(-i)/N
          interpolated_mat.append(temp_mat)
       length_list=len(interpolated_mat)
       interpolated_matrix=interpolated_mat[0]
       for i in range(1,length_list):
          interpolated_matrix=block_matrix([[interpolated_matrix,interpolated_mat[i]]])
       self.mats.append([interpolated_matrix,self.mats[index_mat][1],self.mats[index_mat][2],self.mats[index_mat][3]])

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
#Function generates the encryption keys for left-encoding
#############################
     def gen_lkeys(self,s_key, T): 
         R=random_matrix(k,s_key[0],T*s_key[1])
         sub_ind=[s_key[1]*i for i in range(1,T)]
         R.subdivide([],sub_ind)
         return(R)
#############################
#Function generates the encryption keys for right-encoding
#############################
     def gen_rkeys(self,s_key,T):
         R=random_matrix(k,T*s_key[0],s_key[1])
         sub_ind=[s_key[0]*i for i in range(1,T)]
         R.subdivide(sub_ind,[])
         return(R)
#############################
#Function left-encodes the matrix stored in the source instance
#############################
     def left_encode(self, N,K, T):
         #List of N^th roots of unity
         roots=[k.zeta(N)^i for i in range(N)]
         s_key=[self.size[0],self.size[1]/K]  # size of the keys
         R=self.gen_lkeys(s_key, T)
         sub_ind=[i*s_key[1] for i in range(1,K)]
         self.mat.subdivide([],sub_ind)
         self.encoded_mat=zero_matrix(k,self.size[0],self.size[1]/K*N)
         sub_ind=[s_key[1]*i for i in range(1,N)]
         self.encoded_mat.subdivide([],sub_ind)
         
         for i in range(N):
             submat=zero_matrix(k,self.size[0],self.size[1]/K)
             for j in range(K):
             	submat=submat + self.mat.subdivision(0,j)*roots[i]^j
             if i==0 :
                temp_mat = submat
             else:
                temp_mat=block_matrix([[temp_mat,submat]])
         self.encoded_mat=self.encoded_mat + temp_mat
         sub_ind=[s_key[1]*i for i in range(1,N)]
         self.encoded_mat.subdivide([],sub_ind)
         for i in range(N):
             submat=zero_matrix(k,self.size[0],self.size[1]/K)
             for j in range(K,K+T):
             	submat=submat + R.subdivision(0,j-K)*roots[i]^j
             if i==0 :
                temp_mat = submat
             else:
                temp_mat=block_matrix([[temp_mat,submat]])
         
         self.encoded_mat=self.encoded_mat + temp_mat
         sub_ind=[s_key[1]*i for i in range(1,N)]
         self.encoded_mat.subdivide([],sub_ind)
 
###############################
#Function right-encodes the matrix stored in the source instance
###############################
     def right_encode(self, N,K, T):
         #List of N^th roots of unity
         roots=[k.zeta(N)^i for i in range(N)] 
         s_key=[self.size[0]/K,self.size[1]]
         R=self.gen_rkeys(s_key, T)
         sub_ind=[i*s_key[0] for i in range(1,K)]
         self.mat.subdivide(sub_ind,[])
         self.encoded_mat=zero_matrix(k,s_key[0]*N,s_key[1])
         sub_ind=[s_key[0]*i for i in range(1,N)]
         self.encoded_mat.subdivide(sub_ind,[])
         
         for i in range(N):
             submat=zero_matrix(k,self.size[0]/K,self.size[1])
             for j in range(K):
             	submat=submat + self.mat.subdivision(j,0)*roots[i]^(-j)
             if i==0 :
                temp_mat = submat
             else:
                temp_mat=block_matrix([[temp_mat],[submat]])
         
         self.encoded_mat=self.encoded_mat + temp_mat
         sub_ind=[s_key[0]*i for i in range(1,N)]
         self.encoded_mat.subdivide(sub_ind,[])
         for i in range(N):
             submat=zero_matrix(k,self.size[0]/K,self.size[1])
             for j in range(K+T,K+2*T):
             	submat=submat + R.subdivision(j-K-T,0)*roots[i]^(-j)
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
     def uploadtonodes(self,clstr,K,T):
         N=len(clstr.nodes)
         subdiv=self.encoded_mat.subdivisions()
         if subdiv[0]==[]:
             for i in range(N):
                 clstr.nodes[i].mats.append([self.encoded_mat.subdivision(0,i),N,K,T])
         else:
             for i in range(N):
                 clstr.nodes[i].mats.append([self.encoded_mat.subdivision(i,0),N,K,T])





####################################################################
# The main program code acts as the user:
####################################################################

q=7^3
N=18
T=6
K=N-2*T
gamma=3 # number of source matrices
k.<x>=GF(q) # Fixes a Finite Field of size q
# Note: q-1 must be divisble by N, so that all N^th roots of unity lie in this finite field

###################################################
## This section generates the source matrices, encodes them at the source, and uploads the encoded shares to the computing cluster
###################################################

A=Source([60,60])
A.gen_data() # Source A generates its data
A.left_encode(N,N-2*T,T)
B=Source([60,60])
B.gen_data()  # Source B generates its data
B.right_encode(N,N-2*T,T)
C=Source([60,60])
C.gen_data()  # Source C generates its data
C.right_encode(N,N-2*T,T)

comp_cluster=Cluster(N)
A.uploadtonodes(comp_cluster,N-2*T,T)
B.uploadtonodes(comp_cluster,N-2*T,T)
C.uploadtonodes(comp_cluster,N-2*T,T)


######################################################
##### This section illustrates the computation of matrix product A*B*C, and experimentation with other operations.
######################################################
for i in range(N):
	comp_cluster.nodes[i].multiply(0,1)  ## computes shares of A*B on each server

#########
# right_to_left_share, to generate (N,N-2T,T) left-shares of A*B from (N,1,T) shares of A*B
#########
for i in range(N):
   comp_cluster.nodes[i].left_encode(3,N,K,T)
comp_cluster.communication_phase(4)
for i in range(N):
   comp_cluster.nodes[i].right_interpolate(5,N,1)
#########

for i in range(N):
	comp_cluster.nodes[i].multiply(6,2)  ## multiplies the left-shares of A*B with right-shares of C

######### Demonstration ###########
# right_to_left_share, to generate (N,N-T,T) left-shares of A*B*C from (N,1,T) shares of A*B*C
#########
for i in range(N):
   comp_cluster.nodes[i].left_encode(7,N,N-T,T)
comp_cluster.communication_phase(8)
for i in range(N):
   comp_cluster.nodes[i].right_interpolate(9,N,1)
#########
################# The servers may transmit these (N,N-T,T) shares of A*B*C to the user to achieve a download cost of N/(N-T). ########################

######### Demonstration ###########
# left_to_right_share, to generate (N,1,T) shares of A*B*C from (N,N-T,T) left-shares of A*B*C
#########
for i in range(N):
   comp_cluster.nodes[i].right_encode(10,N,1,T)
comp_cluster.communication_phase(11)
for i in range(N):
   comp_cluster.nodes[i].left_interpolate(12,N,12)
#########

############ If the user receives (N,1,T) shares of the result, then the user simply averages the received matrices. The download cost in this case is N. #############

C_multiply=comp_cluster.nodes[0].mats[13][0]/N
for i in range(1,N):
	C_multiply=C_multiply + comp_cluster.nodes[i].mats[13][0]/N

verify=C_multiply - A.mat*B.mat*C.mat   # print verify to check that the reconstructed result is correct. Should be a zero matrix if the implementation is correct


## Following section commented out.
'''  
#####################################################
## This section demonstrates the computation of the matrix product B*A. The working of the functions right_to_left_share() and left_to_right_share() is illustrated through this example.
#####################################################

#### Convert left-shares of A to right-shares

for i in range(N):
    comp_cluster.nodes[i].right_encode(0,N,K,T)
comp_cluster.communication_phase(3) 
for i in range(N):
    comp_cluster.nodes[i].left_interpolate(4,N,K)

#### Convert right-shares of B to left-shares

for i in range(N):
    comp_cluster.nodes[i].left_encode(1,N,K,T)
comp_cluster.communication_phase(6) 
for i in range(N):
    comp_cluster.nodes[i].right_interpolate(7,N,K)
########################

for i in range(N):
    comp_cluster.nodes[i].multiply(8,5)  ## multiplies the left-shares of A*B with right-shares of C

C_multiply=comp_cluster.nodes[0].mats[9][0]/N
for i in range(1,N):
	C_multiply=C_multiply + comp_cluster.nodes[i].mats[9][0]/N

verify=C_multiply - B.mat*A.mat   # print verify. Should be a zero matrix if the implementation is correct

'''
