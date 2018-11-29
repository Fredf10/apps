'''
Created on Oct 5, 2016

@author: fredrik
'''



from classCreateSimpleApp import SetupApp

from bokeh.layouts import Row, column
from bokeh.io import curdoc


app = SetupApp()

for w in app.Widgetlist:
    w.on_change('value', app.update_data)


# Set up layouts and add to document
# inputs = column(children=app.Widgetlist)
# outputs = column(children=[app.plot_radius, app.plot_pressure])
# curdoc().add_root(Row(children=[inputs, outputs], width=1500))

inputs = column(children=app.Widgetlist)
FigRow1 = column(children=[app.plot_radius, app.plot_pressure])
FigRow2 = column(children=[app.plot_velocity_pos, app.plot_pressure_pos])
curdoc().add_root(Row(children=[inputs, FigRow1, FigRow2], width=2000))