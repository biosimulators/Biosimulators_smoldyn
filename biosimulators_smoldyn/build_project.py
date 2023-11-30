from smoldyn import Simulation


model_fp = 'biosimulators_smoldyn/model.txt'
simulation = Simulation.fromFile(model_fp)

simulation.runSim()


