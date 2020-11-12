# imports extern
import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output
# import tkinter
import matplotlib
# import plotly
# from plotly.offline import plot_mpl
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2  # , venn2_circles


# imports intern
import src.data_calculation as data
import src.argument_parser as pars


# start dash
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


# Values for chi² results
df_chi = pd.DataFrame({
        "Typ": ["Test Statistic", "P-Value", "Degree Of Freedom"],
        "Value": [data.chi_results[0], data.chi_results[1], data.freedom]
})


# --------------------------------------------
def percent_string(number):
    clean = round(number * 100, 2)
    return "~" + str(clean) + "%"


# visualisation of output data
# dataframe for output
df = pd.DataFrame({
    "Degree of Overlap in %": ["Quotient", "First File", "Second File", "Quotient", "First File", "Second File"],
    "Coverage": [data.both_files_lazy, data.first_file_lazy, data.second_file_lazy,
                 data.both_files_log, data.first_file_log, data.second_file_log],
    "Size": ["Actual Value", "Actual Value", "Actual Value",
             "Log 2 Value", "Log 2 Value", "Log 2 Value"]
})

# dataframe output regular
df_reg = pd.DataFrame({
    "Degree of Overlap in %": ["Intersection", "First File", "Second File"],
    "Coverage": [data.both_files_lazy, data.first_file_lazy, data.second_file_lazy],
    "Size": ["Actual Value", "Actual Value", "Actual Value"],
})
reg_over = ["First File", "Intersection", "Second File"]
reg_cov = [data.first_file_lazy, data.both_files_lazy, data.second_file_lazy]
reg_anot = [percent_string(data.first_file_lazy),
            percent_string(data.both_files_lazy),
            percent_string(data.second_file_lazy)]


# dataframe output log
df_log = pd.DataFrame({
    "Degree of Overlap in %": ["Intersection", "First File", "Second File"],
    "Coverage": [data.both_files_log, data.first_file_log, data.second_file_log],
    "Size": ["Log 2 Value", "Log 2 Value", "Log 2 Value"]
})
log_over = ["First File", "Intersection", "Second File"]
log_cov = [data.first_file_log, data.both_files_log, data.second_file_log]
log_anot = [percent_string(data.first_file_log),
            percent_string(data.both_files_log),
            percent_string(data.second_file_log)]


# writing results to output file
if pars.args.outfile is not None:
    filename = pars.args.outfile + ".csv"
    df.to_csv(filename, index=False)
# https://medium.com/@ageitgey/python-3-quick-tip-the-easy-way-to-deal-with-file-paths-on-windows-mac-and-linux-11a072b58d5f


# figure
# TODO change to two(normal and log) Venn-Diagram for better readability
fig = px.bar(df, x="Degree of Overlap in %", y="Coverage", color="Size", barmode="group")
# fig_reg = px.bar(df_reg, x="Degree of Overlap in %", y="Coverage", color="Size", barmode="group", labels="Coverage")
# fig_log = px.bar(df_log, x="Degree of Overlap in %", y="Coverage", color="Size", barmode="group", labels="Coverage")
fig_reg = go.Figure(data=[go.Bar(
    x=reg_over, y=reg_cov, text=reg_anot, textposition='auto'
)])
fig_log = go.Figure(data=[go.Bar(
    x=log_over, y=log_cov, text=log_anot, textposition='auto'
)])


v = venn2((data.len_one, data.len_two, data.len_inter), set_labels=('Group A', 'Group B'))
matplotlib.pyplot.savefig('VennDiagram.png')
# a hint, how to convert venn diagram so that it can be shown in browser like the other graphs
# https://stackoverflow.com/questions/49851280/showing-a-simple-matplotlib-plot-in-plotly-dash

# conversions fails
# v1 = plot_mpl(v)
# v2 = plotly.tools.mpl_to_plotly(v1, resize=True, strip_style=True, verbose=True)


# app layout
app.layout = html.Div(children=[

    # title for the webpage
    html.H1(children="Overlap between both BED-files", style={'text-align': 'center'}),

    html.Div(id='output_container', children=[]),

    # dcc.Graph(id='venn1', figure=v),

    # following layout error message
    # dcc.Graph(id='venn1', figure=matplotlib.pyplot.gcf()),
    # dcc.Graph(id='venn1', plt.gcf()),
    # plt.figure(),

    # just opens new window
    # plt.show(),
    # try impeding as image
    # html.Img(src=app.get_asset_url('/src/VennDiagram.png'), style={'height': '60%', 'width': '60%'}),
    # html.Br(),


    # graph
    html.H6("Actual percent values of the degree of overlap"),
    dcc.Graph(id='overlap_reg', figure=fig_reg),

    # html.H6("Degree of overlap adjusted with log 2 as to compensate for different file sizes"),
    # dcc.Graph(id='overlap_log', figure=fig_log),


    #dcc.Graph(id='overlap_files', figure=fig),
    #html.H6("Shown is 1. degree overlap between both files in reference to their "
    #        "combined length, 2. degree overlap of first file in reference to it's "
    #        "length and 3. degree overlap of second file in reference to it's length "),

    #dcc.Graph(
    #    id='overlap-graph',
    #    figure={
    #        'data': [
    #            {'x': [1, 2, 3], 'y': [data.both_files_lazy, data.first_file_lazy, data.second_file_lazy],
    #             'type': 'bar', 'name': 'Actual Value'},
    #            {'x': [1, 2, 3], 'y': [data.both_files_log, data.first_file_log, data.second_file_log],
    #             'type': 'bar', 'name': u'Log 2 Value'},
    #        ],
    #        'layout': {
    #            'title': 'Dash Data Visualization',
    #            'xaxis': {'type': 'text', 'title': "Coverage"},
    #            'yaxis': {'type': 'linear', 'title': "Degree of Overlap in %"}
    #        }
    #    }
    #),
    #html.H6("Shown is 1. degree overlap between both files in reference to their "
    #        "combined length, 2. degree overlap of first file in reference to it's "
    #        "length and 3. degree overlap of second file in reference to it's length "),

    html.Br(),

# input for alpha and degrees of freedom
    html.H6(children='Degree of Freedom: '),
    dcc.Slider(
        id='freedom_slider',
        min=0,
        max=2,
        value=data.freedom,
        # marks={'1': '1', '5': '5', '10': '10', '15': '15',  '20': '20'}
    ),
    html.H6(children='Alpha-Value: '),
    # slider for alpha value
    # dropdown for alpha value
    # TODO field for custom alpha value
    dcc.Dropdown(
        id='alpha_dropdown',
        options=[
            {'label': '0.01%', 'value': '0.0001'},
            {'label': '0.1%', 'value': '0.001'},
            {'label': '1%', 'value': '0.01'},
            {'label': '5%', 'value': '0.05'},
            {'label': '10%', 'value': '0.1'}

        ],
        value=data.alpha
    ),

    html.H6("Change these to see if the result for files being 'statistical significant different' changes"),

    html.Br(),

    # html.H3(children=validated_own, style={'text-align': 'left'}),
    html.H2(id='val_text', children=data.chi_text, style={'text-align': 'left'}),
    html.H4(id='stat_text', children='Test Statistic: ' + str(data.chi_results[0]), style={'text-align': 'left'}),
    html.H4(id='p_val', children='P-Value: ' + str(data.chi_results[1]), style={'text-align': 'left'}),
    html.H4(id='alp', children='Alpha: ' + str(data.alpha), style={'text-align': 'left'}),
    html.H4(id='freed', children='Degree Of Freedom: ' + str(int(data.freedom)), style={'text-align': 'left'})

])


# update when new freedom or alpha detected
@app.callback([
    Output('val_text', 'children'),
    # Output('stat_text', 'children'),  new_val[0],
    Output('p_val', 'children'),
    Output('alp', 'children'),
    Output('freed', 'children')
    ],
    [
    Input('alpha_dropdown', 'value'),
    Input('freedom_slider', 'value')
    ]
)
def update_text(alpha_dropdown, freedom_slider):

    new_chi = data.calc_chi(float(freedom_slider), float(alpha_dropdown))
    new_val = new_chi[0]  # chi² results
    new_text = new_chi[1]  # validation text
    p_text = 'P-Value: ' + str(new_val[1])
    a_text = 'Alpha: ' + str(alpha_dropdown)
    f_text = 'Degree Of Freedom: ' + str(freedom_slider)
    return new_text, p_text, a_text, f_text


def run_gui():
    app.run_server(debug=True)
# run app
# app.run_server(debug=True)  # debug=true   means update browser by code change
# see in browser http://127.0.0.1:8050/ ONLY development server
