from smoldyn import Simulation


model_fp = 'biosimulators_smoldyn/Andrews-Min1.txt'
simulation = Simulation.fromFile(model_fp)

simulation.runSim()


