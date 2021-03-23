#dc_piece_theta2.mod. 
#Implements a piecewise linear constrained version of the problem using
#AMPL's built in piecewise linear capabilities. Includes voltage
#angles.


set NODES;      #nodes 
set GENERATORS;  #generations stations

set NODE_GEN within {NODES, GENERATORS};  #assign generators to nodes

param LOADS {NODES} >= 0; #demands at each of the nodes

param GCOST {GENERATORS};  #marginal cost of generation

param GCAP {GENERATORS}; #Capacity of generation

set LINK within {NODES, NODES}; #branches (arcs) in the network

set LINK2 within {NODES, NODES}; #sort of an AC only node-arc incidence martix
                                 #except only indicates the tail.  Used to calculate
                                 #voltage angles

#set LOOP within {NODES, NODES}; #nodes which form a loop (circuit)

param X {LINK}; #Line reactances

param R {LINK}; #Line resistances

param B {(i,j) in LINK} :=
        -(X[i,j])/((R[i,j])^2 + (X[i,j])^2); #admittance

param U {LINK}; # Line capacities (MW)


              # set up the piecewise linear function information:

param npieces {LINK}; # number of segments in the piecewise linear function
                      #  for each link

param breakpoint {(i,j) in LINK,1..npieces[i,j]-1}; #breakpoints for each of
                                                    #the loss functions

param slope {(i,j) in LINK,1..npieces[i,j]};      #slopes for each of
                                                  #the loss function segments
     # note that there is one less breakpoint than slope, this is because it is
     # assumed that zero is always a breakpoint



# variables
var P {(i,j) in LINK} >=0, <= U[i,j];     # power flow 
var FS {(i,j) in LINK} >=0, <= U[i,j];    # sending power
var FR {(i,j) in LINK} >=0, <= U[i,j];    # receiving power
var DISPATCH {j in GENERATORS} >=0, <=GCAP[j];  # power dispatched by a station
var GEN {NODES} >=0;                          # total power injection into a bus
var SHEDPOWER {NODES} >=0;                #amount of power shed at a node
var THETA {NODES} >=0;                    #voltage angle
var d {(i,j) in LINK} binary;       #indicates if there is any power flow 
                                    #between node i and j, used to prevent
                                    #cyclical flow
minimize total_cost:
         sum { j in GENERATORS } GCOST[j]*DISPATCH[j];

subject to GENERATION {i in NODES}:
           GEN[i] = sum {(i,j) in NODE_GEN} DISPATCH[j];
         
subject to USERS {i in NODES}:
           LOADS[i] = GEN[i] + sum {(i,j) in LINK} (FR[j,i] - FS[i,j]);

subject to VOLTAGE_ANGLES {(i,j) in LINK2}:
           P[i,j] - P[j,i] = -(THETA[i] - THETA[j])*B[i,j];

subject to RLOSSES {(i,j) in LINK}:
           FR[i,j] = P[i,j] - <<{p in 1..npieces[i,j]-1} breakpoint[i,j,p];
                               {p in 1..npieces[i,j]} slope[i,j,p]>>P[i,j]*0.5;

subject to SLOSSES {(i,j) in LINK}:
           FS[i,j] = P[i,j] + <<{p in 1..npieces[i,j]-1} breakpoint[i,j,p];
                               {p in 1..npieces[i,j]} slope[i,j,p]>>P[i,j]*0.5;

  #sum of the reactance/B in a loop must equal zero
#subject to CYCLE:
#           sum {(i,j) in LOOP} ((P[i,j]-P[j,i])/B[i,j]) = 0;

 # the followingn two constraints prevent cyclical flow
subject to NON_CYCLICAL_FLOW1 {(i,j) in LINK}:
           P[i,j] <= d[i,j]*U[i,j];
subject to NON_CYCLICAL_FLOW2 {(i,j) in LINK}:
           d[i,j] + d[j,i] <=1;

#subject to SHEDDING {i in NODES}:
#           SHEDPOWER[i] = -LOADS[i]+ GEN[i] + 
#                            sum {(i,j) in LINK} (FR[j,i] - FS[i,j])                                                  
