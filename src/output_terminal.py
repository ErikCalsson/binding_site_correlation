# imports extern
import pandas as pd

# imports intern
import src.data_calculation as data
import src.argument_parser as pars

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


# start terminal output
def use_terminal():
    print(data.chi_text)
    print(str(data.chi_results[0]))
    print(str(data.chi_results[1]))
    print(str(data.alpha))
    print(str(int(data.freedom)))
    print(df)
