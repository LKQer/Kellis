import numpy as np

# sigmoid function
def nonlin(x,deriv=False):
    if(deriv==True):
        return x*(1-x)
    return 1/(1+np.exp(-x))
    
def update_syn1(syn0):
    return [syn0[1], syn0[0], syn0[3], syn0[2], syn0[5], syn0[4]]

# # input dataset
# X = np.array([  [0,0,1,0],
#                 [0,1,1,0],
#                 [1,0,0,1],
#                 [1,1,0,0],
#                 [0,0,1,1],
#                 [0,1,0,1],
#                 [1,0,1,0],])
    
# # output dataset            
# y = np.array([[0,1,1,0,0,0,0]]).T

X = np.array([  [0,0,1,0,0,0],
                [0,1,1,0,1,0],
                [1,0,0,1,0,0],
                [1,0,1,0,1,1],
                [0,0,0,0,0,1],
                [1,0,0,0,1,0]
                ])
    
# output dataset            
y = np.array([[0,1,1,0,0,0]]).T

# seed random numbers to make calculation
# deterministic (just a good practice)
np.random.seed(2)

# initialize weights randomly with mean 0
syn0 = 2*np.random.random((6,1)) - 1
# syn0[0] = -5
# syn0[1] = 5
# syn0[2] = 5
# syn0[3] = -5
# syn0[4] = -0.25
# syn0[5] = -0.25
# syn1 = update_syn1(syn0)
syn1 = 2*np.random.random((6,1)) - 1

for iter in xrange(10000):

    # forward propagation
    l0 = X
    l1 = nonlin(np.dot(l0,syn0))
    l2 = nonlin(np.dot(l0,syn1))
    l = l1
    source = []
    for i in range(len(l)):
        l[i] = max(l1[i], l2[i])
        if l[i] == l1[i]:
            source.append('l1')
        else:
            source.append('l2')

    # how much did we miss?

    l1_error = y - l
    if iter % 1000 == 0:
        print "Error:" + str(np.mean(np.abs(l1_error)))
    l2_error = y - l

    for i in range(len(source)):
        if source[i] != 'l1':
            l1_error[i] = 0
        if source[i] != 'l2':
            l2_error[i] = 0
    # print source
    # print l1_error
    # print l2_error

    # multiply how much we missed by the 
    # slope of the sigmoid at the values in l1
    l1_delta = l1_error * nonlin(l1,True)
    l2_delta = l2_error * nonlin(l2,True)

    # update weights
    syn0 += np.dot(l0.T, l1_delta)
    syn1 += np.dot(l0.T, l2_delta)
    # syn1 = update_syn1(syn0)


print "Output After Training:"
print 'l1:', l1
print 'l2:', l2
print 'l:', l
print 's0:', syn0
print 's1:', syn1