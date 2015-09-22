# CPMG_SIM

#An octave function, will simulate CPMG RD trace for simple two-state exchange or for the full CypA catalytic cycle as described in Chapter IV:

#function [dat,vcpmg,track_t,track] = FGP_cpmg_sim(K,Cyp_Conc=1e-3,FGP_Conc=6e-3,w_bound=[45,90],time_T2=20,ncyc=[2:2:20],SIMPLE_CPMG=0,dead_enzyme=0,Natoms=1e5)

#K = Microscopic rate constants defining the catalytic cycle: [kab, kba, kbc, kcb, kcd, kdc] If SIMPLE_CPMG==1, only kbc and kcb will be used
#Cyp_Conc = concentration of CypA in M
#FGP_Conc = concentration of peptide in M
#w_bound = chemical shift difference of trans and cis from unbound: [w_trans, w_cis] 
#time_T2 = delay time, in ms
#ncyc = number of refocusing pulses to be applied, can be a single value or vector of values
#SIMPLE_CPMG = flag, when set to 1, will simulate a simple 2-state exchange, using the values of kbc and kcb as the forward and reverse rate constants, repectively
#dead_enzyme = flag, when set to 1, will simulate data with kbc=kcb=0
#Natoms = number of atoms to simulate for each ncyc value
#dat = vector of simulated R2 values for each vcpmg (s^-1)
#vcpmg = cpmg refocusing frequencies (s^-1)
#track = state of atom #1 throughout the simulation
#track_t = time corresponding to each track value (ms)

