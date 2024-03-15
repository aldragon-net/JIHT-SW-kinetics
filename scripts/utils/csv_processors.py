import pandas as pd
from pathlib import PurePath
from configs.species_slices import Y2_PAH_AND_SOOT, Y5_SOOT

def sum_columns(filepath,
                old_columns=['grid', 'velocity', 'T'],
                new_columns=[[]]):
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


sum_columns(
    'output/stagnation_flame_data/ethylene_base_Y.csv',
    new_columns={'Y2': Y2_PAH_AND_SOOT, 'Y5': Y5_SOOT}
)
