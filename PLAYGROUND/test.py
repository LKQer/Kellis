import numpy as np

# sigmoid function
def nonlin(x,deriv=False):
    if(deriv==True):
        return x*(1-x)
    return 1/(1+np.exp(-x))
    
# input dataset
X = np.array([  [0,0,1,0,1,0,0,0],
                [0,1,1,0,0,0,1,1],
                [1,0,0,1,1,0,1,0],
                [1,0,1,0,1,0,1,0],
                [0,0,0,1,0,1,0,0],
                [1,0,0,1,0,0,1,1],
                [0,1,1,0,0,1,0,1],
                [0,1,0,1,0,1,0,1],
                 ])
    
# output dataset            
y = np.array([[0,1,1,0,0,1,1,0]]).T

# seed random numbers to make calculation
# deterministic (just a good practice)
np.random.seed(2)

# initialize weights randomly with mean 0
syn0 = 2*np.random.random((8,1)) - 1
syn1 = 2*np.random.random((8,1)) - 1

for iter in xrange(10000):

    # forward propagation
    l0 = X
    l1 = nonlin(np.dot(l0,syn0))
    l2 = nonlin(np.dot(l0,syn1))


    # Acts as a max layer
    final = np.ndarray(shape=(8,1))
    for i in range(len(l1)):
        if l1[i] > l2[i]:
            l2[i] = y[i]
            final[i] = l1[i]
        else:
            l1[i] = y[i]
            final[i] = l2[i]
    # print l1, l2

    # how much did we miss?
    l1_error = y - l1
    l2_error = y - l2

    # multiply how much we missed by the 
    # slope of the sigmoid at the values in l1
    l1_delta = l1_error * nonlin(l1,True)
    l2_delta = l2_error * nonlin(l2,True)

    # update weights
    syn0 += np.dot(l0.T,l1_delta)
    syn1 += np.dot(l0.T,l2_delta)

print "Output After Training:"
print l1
print l2
print final
print syn0
print syn1