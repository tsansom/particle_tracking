# Particle Tracking in TxBLEND with Applications in the Trinity-San Jacinto and Mission-Aransas Estuaries

To support oyster reef restoration projects and other potential transport projects along the Texas coast, I developed model-level particle tracking functionality for TxBLEND, a two-dimensional hydrodynamic and salinity transport model used by the Texas Water Development Board (TWDB) in support of regional water planning and development of environmental flow regime recommendations. TxBLEND is a vertically-averaged, finite-element model which employs an unstructured, triangular element grid with linear basis functions. The model was designed to simulate water levels, water circulation, and salinity conditions in Texas estuaries, which are generally very shallow and wide bodies of water.

Particle tracking studies are often conducted using velocity components from model output, which typically have time scales of 30 minutes to one hour. The novelty of implementing particle tracking at the model-level is that the calculation of velocity, and subsequently particle position, happen on the model time-step. The typical time-step for TxBLEND model runs is two or three minutes. This allows us to capture smaller perturbations of the velocity field that would otherwise be lost to averaging when using hourly model output for example.

---
