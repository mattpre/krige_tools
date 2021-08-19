import numpy as np
import krige_tools as kt

prob = 'boreholes_wF5'

strat = kt.Stratigraphy()

strat.layer_names = [#'6a',
                     '6 retrait',
##                     '6e',
##                     '7sup',
                     '7 moraine',
                     '9',
                     '12',
                     '15 molasse']

##strat.read_boreholes_from_file('bh_test.txt','.')
strat.read_boreholes_from_file(prob+'.txt','.')
strat.bounds = [[-40,130], [-190,45,]]
strat.grid_size = 5
strat.additional_points.append([8,[-50,120],385])
strat.additional_points.append([9,[-50,120],380])
strat.additional_points.append([10,[-50,120],375])
strat.additional_points.append([8,[-50,-40],385])
strat.additional_points.append([9,[-50,-40],380])
strat.additional_points.append([10,[-50,-40],375])
strat.range[2] = 200
strat.range[3] = 200
strat.range[4] = 200
strat.create_layers()
    
##strat.regularize_layers()
strat.updir = 'z'
##strat.write_layers(output='individual',prob=prob,pathname='pv')
strat.write_layers(output='global',prob=prob,pathname='pv')
strat.write_boreholes(prob=prob,pathname='pv')

strat.check_stratigraphy()

strat.create_inclusions(boreholes=['F1','F1','F2','F4','F5'],
                        data=[[4.7,6,10],
                             [6.9,7.4,8],
                             [3.5,4.2,10],
                             [7.3,14,20],
                             [10,13.6,15]],name='6e')

##strat.create_inclusions(boreholes=['F1','F1','F2','F4','F5'],
##                        data=[[4.7,6,10],
##                             [6.9,7.4,8],
##                             [3.5,4.2,10],
##                             [7.3,14,20],
##                             [10,13.6,15]],name='6e')

    
