# imports extern
import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px

# imports intern
import src.data_calculation as data
import src.argument_parser as pars


# start dash
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

validated = 'Values'
validated += data.val_str  # result if p-value is bigger or smaller than chi²-value
validated += 'statistical significant different with'
validated += ':'


# Values for chi² results TODO remove or display
df_chi = pd.DataFrame({
        "Typ": ["Test Statistic", "P-Value", "Degree Of Freedom"],
        "Value": [data.chi_results[0], data.chi_results[1], data.freedom]
})


# --------------------------------------------
# visualisation of output data

# dataframe for output
df = pd.DataFrame({
    "Overlap": ["Quotient", "First File", "Second File", "Quotient", "First File", "Second File"],
    "Coverage": [data.both_files_lazy, data.first_file_lazy, data.second_file_lazy,
                 data.both_files_log, data.first_file_log, data.second_file_log],
    "Size": ["Absolute", "Absolute", "Absolute",
             "Log 2", "Log 2", "Log 2"]
})


# writing results to output file
if pars.args.outfile is not None:
    filename = pars.args.outfile + ".csv"
    df.to_csv(filename, index=False)
#https://medium.com/@ageitgey/python-3-quick-tip-the-easy-way-to-deal-with-file-paths-on-windows-mac-and-linux-11a072b58d5f


#https://plotly.com/python/creating-and-updating-figures/
#https://dash.plotly.com/dash-core-components/graph

#https://dash.plotly.com/basic-callbacks
# figure
fig = px.bar(df, x="Overlap", y="Coverage", color="Size", barmode="group")

# app layout
app.layout = html.Div(children=[

    # title for the webpage
    html.H1(children="Overlap between both BED-files", style={'text-align': 'center'}),

    html.Div(id='output_container', children=[]),

    html.Br(),

    dcc.Graph(id='overlap_files', figure=fig),

    html.Br(),

    # html.H3(children=validated_own, style={'text-align': 'left'}),
    html.H2(children=validated, style={'text-align': 'left'}),
    html.H4(children='Test Statistic: ' + str(data.chi_results[0]), style={'text-align': 'left'}),
    html.H4(children='P-Value: ' + str(data.chi_results[1]), style={'text-align': 'left'}),
    html.H4(children='Alpha: ' + str(data.alpha), style={'text-align': 'left'}),
    html.H4(children='Degree Of Freedom: ' + str(int(data.freedom)), style={'text-align': 'left'})

])


def run_gui():
    app.run_server(debug=True)
# run app
# app.run_server(debug=True)  # debug=true   means update browser by code change
# see in browser http://127.0.0.1:8050/ ONLY development server
