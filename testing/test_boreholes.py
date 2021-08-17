import numpy as np
import krige_tools as kt

strat = kt.Stratigraphy()

##strat.read_boreholes_from_file('bh_test.txt','.')
strat.read_boreholes_from_file('boreholes_all_20160308_isoGeotest.txt','.')
##strat.bounds = [[684910.0, 685110.0], [-256040.0, -255940.0]]
strat.grid_size = 8
strat.additional_points.append([5,[684912.73, 255900],350])
strat.create_original_layers()
##strat.regularize_layers()
strat.updir = 'y'
strat.write_layers(original=True,output='global')
strat.write_boreholes()

strat.check_stratigraphy()
