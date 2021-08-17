import numpy as np
import krige_tools as kt

prob = 'boreholes_all_20160308_isoGeotest'
prob = 'boreholes_all_20160707_isoGeotest'

strat = kt.Stratigraphy()

strat.layer_names = ['OK_Deckschicht',
                     'OK_JSchotter',
                     'OK_JSeeabl',
                     'OK_Moraene',
                     'OK_Schotter',
                     'OK_Seeabl']

##strat.read_boreholes_from_file('bh_test.txt','.')
strat.read_boreholes_from_file(prob+'.txt','data')
##strat.bounds = [[684910.0, 685110.0], [-256040.0, -255940.0]]
strat.grid_size = 8
strat.additional_points.append([5,[684912.73, 255950],340])
strat.additional_points.append([6,[684912.73, 255950],330])
strat.additional_points.append([5,[684952.73, 255930],340])
strat.additional_points.append([6,[684952.73, 255930],330])
strat.additional_points.append([6,[684970.73, 256270],380])
strat.add_points('boreholes_all_20160707_addPts.txt','data')
strat.create_original_layers()
    
##strat.regularize_layers()
strat.updir = 'z'
strat.write_layers(original=True,output='individual',prob=prob,pathname='pv')
strat.write_layers(original=True,output='global',prob=prob,pathname='pv')
strat.write_boreholes(prob=prob,pathname='pv')

strat.check_stratigraphy()
