import sys

sys.path.append('D:\Dropbox\Code\Python\kinetic\Pynetics')
from pynetics import model

#create micro kinetic model instance
m = model.KineticModel(setup_file='ir.mkm')

m.parser.parse_data()  # parse data from rel_energy.py
