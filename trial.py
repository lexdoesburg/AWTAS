import numpy as np

# parameter space
xlim = [0.01, 0.2]	# phi range
ylim = [1e-16, 1e-12]	# k range


def transf0(theta):
	"""transform parameter vector to a new base adapted for automated calibration

    Args:
		theta (list): parameter vector in original base
    """
	return [(theta[0]-xlim[0])/(xlim[1]-xlim[0]), (theta[1]-ylim[0])/(ylim[1]-ylim[0])]
	
def transf1(X):
	"""transform parameter vector to their original base

    Args:
		X (list): parameter vector in adapted base
    """
	return [X[0]*(xlim[1]-xlim[0])+xlim[0], X[1]*(ylim[1]-ylim[0])+ylim[0]]


params = np.array([0.2, 1e-12])
paramsf0 = transf0(params)
print(paramsf0)
paramsf1 = transf1(paramsf0)
print(paramsf1)

print(params-paramsf1)