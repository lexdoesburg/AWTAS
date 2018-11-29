
#####################################################################################
#
# Function library for plotting various steps in steepest descent lab
#
# 	Functions:
#		plot_dr0: plot initial parameter vector and steepest descent direction
#		plot_step0: plot first step in descent
#		plot_dr1: plot first step and next descent direction
#		plot_steps: plot all steps in descent
# 		plot_steps_wellbore: plot all given steps and corresponding production histories
#
#####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import csv
from wellbore_model import *

J_scale = 0.2
	
def plot_dr0(obj, theta0, dr0):
	"""Plot initial parameter vector and steepest descent direction, overlaid on contours of objective function.

    Args:
		obj (callable): objective function
		theta0 (list): initial parameter vector
		dr0 (list): initial steepest descent direction
    """
	# create plotting grid
	N = 501								# number grid points
	x = np.linspace(-1, 1., N)			# vector in x direction
	y = np.linspace(-1., 1., N)			# vector in y direction
	xv, yv = np.meshgrid(x, y)			# generate meshgrid
	
	# compute objective function on grid
	Z = np.zeros(xv.shape) 				# empty matrix
	for i in range(len(y)):					# for each x value
		for j in range(len(x)):					# for each y value
			Z[i][j] = obj([xv[i][j], yv[i][j]])		# compute objective function
	
	# plotting
	plt.clf()
	ax1 = plt.axes()			# create axes
	ax1.contourf(xv, yv, Z, 21, alpha = .8)		# contours of objective function
	ax1.scatter(theta0[0], theta0[1], color='k', s = 20.)		# show parameter steps
	ax1.arrow(theta0[0], theta0[1], -J_scale*dr0[0], -J_scale*dr0[1])	# show descent direction
	
	# plotting upkeep
	ax1.set_xlim(-1.,1.)
	ax1.set_ylim(-1.,1.)
	ax1.set_xlabel('x')
	ax1.set_ylabel('y')
	ax1.set_aspect('equal')
	
	# save and show
	plt.savefig('lab2_plot1.png', bbox_inches = 'tight')
	plt.show()
	
def plot_step(obj, theta0, dr, theta1):
	"""Plot first step in descent, overlaid on contours of objective function.

    Args:
		obj (callable): objective function
		theta0 (list): initial parameter vector
		dr0 (list): initial steepest descent direction
		theta1 (list): parameter vector at end of first step
    """
	# create plotting grid
	N = 501								# number grid points
	x = np.linspace(-1, 1., N)			# vector in x direction
	y = np.linspace(-1., 1., N)			# vector in y direction
	xv, yv = np.meshgrid(x, y)			# generate meshgrid
	
	# compute objective function on grid
	Z = np.zeros(xv.shape) 				# empty matrix
	for i in range(len(y)):					# for each x value
		for j in range(len(x)):					# for each y value
			Z[i][j] = obj([xv[i][j], yv[i][j]])		# compute objective function

	# plotting
	plt.clf()
	ax1 = plt.axes()		# create axes
	ax1.contourf(xv, yv, Z, 21, alpha = .8)		# contours of objective function
	# show parameter steps
	ax1.scatter([theta0[0], theta1[0]], [theta0[1], theta1[1]], color='k', s = 20.)
	ax1.plot([theta0[0], theta1[0]], [theta0[1], theta1[1]], color='k', linestyle = '--')
	# show descent direction
	ax1.arrow(theta0[0], theta0[1], -J_scale*dr[0], -J_scale*dr[1])
	
	# plotting upkeep
	ax1.set_xlim(-1.,1.)
	ax1.set_ylim(-1.,1.)
	ax1.set_xlabel('x')
	ax1.set_ylabel('y')
	ax1.set_aspect('equal')
	
	# save and show
	plt.savefig('lab2_plot2.png', bbox_inches = 'tight')
	plt.show()
	
def plot_dr1(obj, theta0, dr0, theta1, dr1):
	"""Plot first step and next descent direction, overlaid on contours of objective function.

    Args:
		obj (callable): objective function
		theta0 (list): initial parameter vector
		dr0 (list): initial steepest descent direction
		theta1 (list): parameter vector at end of first step
		dr0 (list): descent direction at end of first step
    """
	# create plotting grid
	N = 501								# number grid points
	x = np.linspace(-1, 1., N)			# vector in x direction
	y = np.linspace(-1., 1., N)			# vector in y direction
	xv, yv = np.meshgrid(x, y)			# generate meshgrid
	
	# compute objective function on grid
	Z = np.zeros(xv.shape) 				# empty matrix
	for i in range(len(y)):					# for each x value
		for j in range(len(x)):					# for each y value
			Z[i][j] = obj([xv[i][j], yv[i][j]])		# compute objective function

	# plotting
	plt.clf()
	ax1 = plt.axes()		# create axes
	ax1.contourf(xv, yv, Z, 21, alpha = .8)		# contours of objective function
	# show parameter steps
	ax1.scatter([theta0[0], theta1[0]], [theta0[1], theta1[1]], color='k', s = 20.)
	ax1.plot([theta0[0], theta1[0]], [theta0[1], theta1[1]], color='k', linestyle = '--')
	# show descent directions
	ax1.arrow(theta0[0], theta0[1], -J_scale*dr0[0], -J_scale*dr0[1])
	ax1.arrow(theta1[0], theta1[1], -J_scale*dr1[0], -J_scale*dr1[1])
	
	# plotting upkeep
	ax1.set_xlim(-1.,1.)
	ax1.set_ylim(-1.,1.)
	ax1.set_xlabel('x')
	ax1.set_ylabel('y')
	ax1.set_aspect('equal')
	
	# save and show
	plt.savefig('lab2_plot3.png', bbox_inches = 'tight')
	plt.show()
	
def plot_steps(obj, theta_all, dr_all):
	"""Plot all steps in descent, overlaid on contours of objective function.

    Args:
		obj (callable): objective function
		theta_all (list): list of parameter vector updates during descent
		dr_all (list): list of descent directions
    """
	# create plotting grid
	N = 501								# number grid points
	x = np.linspace(-1, 1., N)			# vector in x direction
	y = np.linspace(-1., 1., N)			# vector in y direction
	xv, yv = np.meshgrid(x, y)			# generate meshgrid
	
	# compute objective function on grid
	Z = np.zeros(xv.shape) 				# empty matrix
	for i in range(len(y)):					# for each x value
		for j in range(len(x)):					# for each y value
			Z[i][j] = obj([xv[i][j], yv[i][j]])		# compute objective function

	# plotting
	plt.clf()
	ax1 = plt.axes()		# create axes
	ax1.contourf(xv, yv, Z, 21, alpha = .8)		# contours of objective function
	# show parameter steps
	ax1.scatter([theta[0] for theta in theta_all[1:-1]], [theta[1] for theta in theta_all[1:-1]], color='k', linestyle = '--')
	ax1.scatter(theta_all[0][0], theta_all[0][1], color='b', linestyle = '--', label = 'Initial values')
	ax1.scatter(theta_all[-1][0], theta_all[-1][1], color='g', linestyle = '--', label = 'Final values')
	ax1.plot([theta[0] for theta in theta_all], [theta[1] for theta in theta_all], color='k', linestyle = '--')
	# show decent directions
	for i in range(len(theta_all)-1):
		ax1.arrow(theta_all[i][0], theta_all[i][1], -J_scale*dr_all[i][0], -J_scale*dr_all[i][1])
	
	# plotting upkeep
	ax1.set_xlim(-1.,1.)
	ax1.set_ylim(-1.,1.)
	ax1.set_xlabel('x')
	ax1.set_ylabel('y')
	ax1.set_aspect('equal')
	
	# save and show
	plt.savefig('lab2_plot4.png', bbox_inches = 'tight')
	plt.show()

def plot_steps_wellbore(theta_all, r):
	"""Plot 
	(left) all steps in descent, overlaid on contours of objective function.
	(upright) enthalpy production history, data and modelled for the initial
	and the last two iterations
	(downright) pressure production history, data and modelled for the initial
	and the last two iterations

    Args:
		theta_all (list): list of parameter vector updates during descent
		r (float): objective function value calculated during the last iteration
    """

	# get the iteration number
	i = len(theta_all)-1
	
	# set parameter space
	xlim = [90., 125.]
	ylim = [0.05, .5]

	# create plotting grid
	N = 41									# number grid points
	x = np.linspace(xlim[0], xlim[1], N)	# vector in x direction
	y = np.linspace(ylim[0], ylim[1], N)	# vector in y direction
	xv, yv = np.meshgrid(x, y)				# generate meshgrid
	
	# load pre-calculated objective function on grid
	Z = np.loadtxt(open('contour.csv',"rb"),delimiter=',')
	
	# load production data
	t_data = []						# initiate time vector
	h_data = []						# initiate enthalpy vector
	with open('SKG9D_enth.csv', 'rb') as csvfile:	# open enthalpy production history csv file
		reader=csv.reader(csvfile, delimiter=',')	# open reader
		for row in reader:						# go through each row
			t_data.append(float(row[0]))			# save time
			h_data.append(float(row[1]))			# save enthalpy
	p_data = []						# initiate pressure vector
	with open('SKG9D_press.csv', 'rb') as csvfile:	# open pressure production history csv file
		reader=csv.reader(csvfile, delimiter=',')	# open reader
		for row in reader:						# go through each row								
			p_data.append(float(row[1]))			# save pressure
			
	# compute pressure and enthalpy profiles for the initial and last two iterations
	list_t0, list_h0, list_p0 = wellbore_model(theta_all[0])	# initial guess
	if i > 0:	# check if more than one guess
		list_t2, list_h2, list_p2 = wellbore_model(theta_all[i])	# iteration i
	if i > 1:	# check if more than two guesses
		list_t1, list_h1, list_p1 = wellbore_model(theta_all[i-1])	# iteration i-1

	# create figure
	fig = plt.figure(figsize = [20., 10.])
	
	# left plot (objective function contour)
	ax1 = plt.subplot2grid((2,4), (0,0), colspan=2, rowspan=2)	# open subplot
	ax1.contourf(xv, yv, Z, 100)				# plot contour
	ax1.scatter(theta_all[0][0], theta_all[0][1], color='b', s = 60)	# initial guess
	if i > 0:
		ax1.plot([theta[0] for theta in theta_all[:i+1]], [theta[1] for theta in theta_all[:i+1]], color='k', linestyle = '--')	
		ax1.scatter(theta_all[-1][0], theta_all[-1][-1], color='r', s=60)	# iteration i
	if i > 1:
		ax1.scatter(theta_all[-2][0], theta_all[-2][-1], color='r', s=60, alpha=.4)	# iteration i-1
	if  i > 2:			# plot following iterations if they exist
		ax1.scatter([theta[0] for theta in theta_all[1:-2]], [theta[1] for theta in theta_all[1:i-1]], color='k', s = 60)
	
	# plotting upkeep
	ax1.set_xlabel('Reservoir pressure [bars]')
	ax1.set_ylabel('Reservoir steam fraction [%]')
	ax1.set_xlim(xlim)
	ax1.set_ylim(ylim)
	props = dict(facecolor='w', alpha=1., pad = 10.)
	ax1.text(92., .48, 'Iteration '+str(i)+': P0 = '+str(round(theta_all[-1][0], 2))+' bars, '+str(round(theta_all[-1][1]*100.,2))+'% steam, Objective function = '+str(round(r, 2)), bbox=props)	
		
	# upright plot (enthalpy production history)
	ax2 = plt.subplot2grid((2,4), (0,2), colspan=2)	# open subplot
	ax2.plot(t_data, h_data, c='k', marker = 'o', label = 'Data')
	ax2.plot(list_t0, list_h0, c='b', marker = 'o', label = 'Initial guess')
	if i > 0:
		ax2.plot(list_t2, list_h2, c='r', marker = 'o', label = 'Iteration '+str(i))
	if i > 1:
		ax2.plot(list_t1, list_h1, c='r', marker = 'o', label = 'Iteration '+str(i-1), alpha = .4)
		
	# plotting upkeep
	ax2.set_xlabel('Time [days]')
	ax2.set_ylabel('Enthalpy [kJ/kg]')
	ax2.set_xlim(0., 140.)
	ax2.set_ylim(1400., 2800.)
	ax2.legend(loc=2, framealpha=0.)

	# downright plot (pressure production history)
	ax3 = plt.subplot2grid((2,4), (1,2), colspan=2)	# open subplot
	ax3.plot(t_data, p_data, c='k', marker = 'o')
	ax3.plot(list_t0, list_p0, c='b', marker = 'o')
	if i > 0:
		ax3.plot(list_t2, list_p2, c='r', marker = 'o')
	if i > 1:
		ax3.plot(list_t1, list_p1, c='r', marker = 'o', alpha = .4)	
		
	# plotting upkeep
	ax3.set_xlabel('Time [days]')
	ax3.set_ylabel('Pressure [bars]')
	ax3.set_xlim(0., 140.)
	ax3.set_ylim(0., 120.)	
	
	# Save and show
	plt.tight_layout()
	plt.savefig('lab2_plot5_iteration_'+str(i)+'.png',bbox_inches='tight')
	# plt.show()


	
	