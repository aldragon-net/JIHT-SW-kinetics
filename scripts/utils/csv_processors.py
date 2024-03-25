import pandas as pd
from pathlib import PurePath
from configs.species_slices import BIN_X_SLICES, Y_SOOT_SLICES


def sum_columns(filepath,
                old_columns=['grid', 'velocity', 'T'],
                new_columns=None):
    df = pd.read_csv(filepath)
    new_df = df.copy()[old_columns]
    for sum_name, elements in new_columns.items():
        try:
            new_df[sum_name] = df[elements].sum(axis=1)
        except KeyError:
            new_df[sum_name] = 0
    p = PurePath(filepath)
    new_stem = p.stem + '_sums'
    new_df.to_csv(p.with_stem(new_stem), index=False)


def postprocess_creck_soot(path, flame):
    x_path = path / flame.label / f'{flame.label}_X.csv'
    y_path = path / flame.label / f'{flame.label}_y.csv'
    sum_columns(x_path, new_columns=BIN_X_SLICES)
    sum_columns(y_path, new_columns=Y_SOOT_SLICES)
