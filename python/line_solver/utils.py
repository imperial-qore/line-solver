import line_solver
from line_solver import GlobalConstants, VerboseLevel


import pandas as pd

def tget(df, *args):
    # If no arguments are provided, return the original dataframe
    if not args:
        return df

    # Initialize the mask to select all rows
    mask = pd.Series([True] * len(df))

    # Initialize the columns to select all columns
    columns = df.columns.tolist()

    # Default columns to always include
    default_columns = ['Station', 'JobClass']

    # Iterate through the arguments
    for arg in args:
        if arg in df.columns:
            # If the argument is a column name, set columns to default columns + this column
            columns = default_columns + [arg]
        else:
            # Otherwise, filter the rows where any column matches the argument
            mask = mask & (df.eq(arg).any(axis=1))

    # Apply the mask to filter rows and select the columns
    return df.loc[mask, columns].drop_duplicates()
