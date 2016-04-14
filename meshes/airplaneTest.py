#!/usr/bin/env python

import sys
sys.path.append("/opt/bempp/install/bempp/python/")
#sys.path.append("/opt/fb/bempp/python/")

from bempp.lib import *
import numpy as np
import bempp.lib as lib
import scipy.special as sc
import numpy.polynomial.legendre as leg
import time

accuracyOptions = lib.createAccuracyOptions()
accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2)
accuracyOptions.singleRegular.setRelativeQuadratureOrder(2)
quadStrategy = lib.createNumericalQuadratureStrategy("float64", "complex128", accuracyOptions)

options = lib.createAssemblyOptions()
options.switchToAca(lib.createAcaOptions())
context = lib.createContext(quadStrategy, options)

grid = lib.createGridFactory().importGmshGrid("triangular", "airplane.msh")
pconsts = lib.createPiecewiseConstantScalarSpace(context, grid)

k = 16
#k = 64; #16; #0.16

# ansatz: direction of incident plane wave: assume x here
lhsOp = lib.createHelmholtz3dSingleLayerBoundaryOperator(context, pconsts, pconsts, pconsts, k, "SLP")

# Create a grid function representing the Dirichlet trace of the incident wave
def uIncData(point):
	x, y, z = point
	return np.exp(1j * k * x)

def incident(point):
	x, y, z = point
	res = 0.0*x + 0.0j*y + 0*z
	for pt in range(0,x.size):
		res[pt] = np.exp(1j*k*x[pt]) 
	return res

def evalInc(point):
    x, y, z = point
    if x.size == 1:
	return uIncData(point)
    res = 0.0*x + 0.0j*y + 0.0*z
    for pt in range(0,x.size):
	res[pt] = uIncData([x[pt], y[pt], z[pt] ]) 
    return res

uInc = lib.createGridFunction(context, pconsts, pconsts, evalInc)
rhs = -uInc

# PART 4: Discretize and solve the equations ###################################
solver = lib.createDefaultIterativeSolver(lhsOp)
params = lib.defaultGmresParameterList(1e-8)
solver.initializeSolver(params)
# Solve the equation
solution = solver.solve(rhs)
print solution.solverMessage()

# PART 5: Extract the solution #################################################
sol = solution.gridFunction()
print "************** k = ", k, " **********************"


slPot = lib.createHelmholtz3dSingleLayerPotentialOperator(context, k)
dlPot = lib.createHelmholtz3dDoubleLayerPotentialOperator(context, k)
evalOptions = lib.createEvaluationOptions()

endpl = 0.5;
radpl = 0.1;
hwwing = 0.05;
angl = 0.3;
wings = 0.3;

xe = radpl*np.cos(angl)+wings;
ye = 0.0; #radpl*Sin(angl);
ze = 0.0; # : Middle of wing plane
potRes = slPot.evaluateAtPoints(sol, [[xe],[ye],[ze]], evalOptions)
negInci = -uIncData([xe,ye,ze])
print [[xe],[ye],[ze]], negInci, " = pt and neg of inc wave, approx and potential =", potRes, ", error = ", np.abs(potRes-negInci)


xe, ye, ze = radpl*np.cos(angl), radpl*np.sin(angl), hwwing+endpl/2;
print "Error on airplane body = ", np.abs(slPot.evaluateAtPoints(sol, [[xe],[ye],[ze]], evalOptions) + uIncData([xe,ye,ze]))

#quit()

nPointsX = 301
nPointsY = 301
#xx, yy, zz = np.mgrid[-22.5:22.5:nPointsX*1j, -8.45:17:nPointsY*1j, -0.7:-0.7:1j]
#xx, yy, zz = np.mgrid[-2.1:2.1:nPointsX*1j, 0.0:0.0:1j,  -4.5:4.5:nPointsY*1j]
xx, yy, zz = np.mgrid[-1.75:1.75:nPointsX*1j, 0.0:0.0:1j,  -1.2:1.2:nPointsY*1j]

if True:
	from bempp import visualization as vis
	points = np.vstack((xx.ravel(), yy.ravel(), zz.ravel()))
	#uActor = vis.scalarDataOnRegularGridActor(points, vals, (233,26), transformation="abs")
	#vals = (- slPot.evaluateAtPoints(sol, points, evalOptions) + incident(points) )
	vals = (slPot.evaluateAtPoints(sol, points, evalOptions) + incident(points) )
	#vals = slPot.evaluateAtPoints(sol, points, evalOptions)
	uActor = vis.scalarDataOnRegularGridActor(points, vals, (nPointsX, nPointsY), transformation="abs")
	legendActor = vis.legendActor(uActor)
	gridActor = vis.gridActor(grid)
	vis.plotTvtkActors([uActor, gridActor, legendActor])
	#from tvtk.api import tvtk as tva
	#line = tva.LineSource(point1=(0, 0, 0), point2=(4, 0, 0))
	#line_mapper = tva.PolyDataMapper(input=line.output)
	#line_actor = tva.Actor(mapper=line_mapper)
	#vis.plotTvtkActors([uActor, gridActor, legendActor,line_actor])

