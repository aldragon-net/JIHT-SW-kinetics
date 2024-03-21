import pandas as pd
from pathlib import PurePath
from configs.species_slices import Y2_PAH_AND_SOOT, Y5_SOOT, BIN_X_SLICES, Y_SOOT_SLICES


def sum_columns(filepath, suffix,
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
    sum_columns(x_path, suffix='BINs', new_columns={'test_sum': ['X_O2', 'X_H2']})
    sum_columns(y_path, suffix='soot', new_columns={'test_sum': ['Y_O2', 'Y_H2']})
