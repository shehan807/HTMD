import numpy as np
import scipy.integrate as integrate
import sys

#************************************************
# this script reads a energy derivative output file from a thermodynamic
# integration simulation (OpenMM), and integrates dE_dlambda from 0 to 1
# for the particular lambda scaling
#************************************************

dEdlambdafile=sys.argv[1]

with open(dEdlambdafile) as f:
   data = f.readlines()

# parse data.  average dE_dlambda values for each lambda value.  We might have multiple Hamiltonian scalings (repulsive shells) in
# the same output file
lambda_values=[]  # stores values of lambda used in TI
dEdlambda_data=[] # stores dE_dlambda data
scale_data=[]   # temporary list
lambda_data=[]  # temporary list
for line in data:
    if "lambda =  1.0" in line:
       # start of lambda scaling, see if this is not the first output...
       try:
           # if this exists, we've filled up the list with previous data
           test=scale_data[0]
           # push last lambda data to list
           scale_data.append(lambdavalue_data)
           dEdlambda_data.append(scale_data)
           lambda_values.append(lambda_data)
           # clear lists for new data
           lambda_data=[]
           scale_data=[]
           lambdavalue_data=[]   # temporary list
       except: 
           lambdavalue_data=[]   # temporary list
       # append lambda value
       lambda_data.append( 1.0 )
    elif "lambda" in line:
       # push data list from last lambda value to data array
       scale_data.append(lambdavalue_data)
       # new list for data for this lambda value
       lambdavalue_data=[]
       # now get value of lambda and append
       lambdapoint=line.split()
       lambda_data.append(float(lambdapoint[4]))
    elif "step " in line:
       # push data to list
       datapoint=line.split()
       lambdavalue_data.append(float(datapoint[2]))

# here we've gotten to the end of the data, append the last lists...
scale_data.append(lambdavalue_data)
dEdlambda_data.append(scale_data)
lambda_values.append(lambda_data)

# now loop over dE_dlambda data and average.  Fill in list with averages
avg_dEdlambda=[]

for i in range(len(dEdlambda_data)):
    print( "TI data for Hamiltonian scaling simulation : " , i+1 )
    avg_shell=[]
    for j in range(len(dEdlambda_data[i])):
         # now average
         avg = np.mean( dEdlambda_data[i][j] )
         print( "<dE_dlambda> for lambda value " , lambda_values[i][j], avg )
         avg_shell.append(avg)
    avg_dEdlambda.append(avg_shell)


avg_dEdlambda[0]=avg_dEdlambda[0][::1]
lambda_values[0]=lambda_values[0][::1]
print(avg_dEdlambda)
print(lambda_values)
# now integrate
for i in range(len(dEdlambda_data)):
    print( "integral(0,1) dE_dlambda : " , i+1 )
    delta_E = integrate.simpson(avg_dEdlambda[i], lambda_values[i] )
    # negative sign because data is input lambda 1 ==> lambda 0,
    # and we want integral 0 ==> 1
    delta_E = -delta_E
    print( delta_E , " kJ/mol " )


