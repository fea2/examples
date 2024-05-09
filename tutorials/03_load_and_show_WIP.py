import os
from pathlib import Path

import compas_fea2
from compas_fea2.model import Model
from compas_fea2.problem import Problem

from compas_fea2.units import units
units = units(system='SI_mm')

HERE = os.path.dirname(__file__)
DATA = os.sep.join(HERE.split(os.sep)[:-1]+['data'])
model_name = 'simple_frame'

mdl: Model = Model.from_cfm(Path(DATA).joinpath(model_name+'.cfm'))
prb: Problem = mdl.find_problem_by_name('00_simple_problem')

prb.show_deformed(scale_factor=1000)
prb.show_displacements_contour(draw_loads=1000, draw_bcs=1000, style='vector')