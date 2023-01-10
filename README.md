# ERFLOW

1. Construct network -- choose from C0-C4 or honeycomb

C0: constructCOS7ER.m
C1-C4: constructC.m
Honeycomb: constructhexgraph.m 


2. Compute flows -- specify driving mechanism: contraction of tubules/junctions/peripheral sheets/perinuclear sheets

pinches_pressureBC.m: tubules pinching
pinches_junctions.m: junctions pinching
pinches_ID: minor modification of pinches_junctions.m to investigate of contraction of peripheral sheets (also called "interdispersed sheets")
pinches_sheets.m: contraction of perinuclear sheets; note we have investigated perinuclear sheets only for the C1 network. To investigate this in other networks, we may define a list of sheet nodes as in line 346 of C1DATA.m. 


3. Track motion of Brownian particles subject to flows and process data

manyparticles_COS7.m: tracks particles and processes data


4. Plot data

analysis_velocities: plots histograms of instantaneous speeds, edge traversal speeds, and average edge traversals speeds
plotAETS: plots histogram of average edge traversal speed with a Gaussian fit; straightforward modification allows for lognormal fit, etc. 
