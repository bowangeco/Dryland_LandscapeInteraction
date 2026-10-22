"""Full-period BestE: author-supplied Python analysis executed on an HPC system using MPI."""

DATA_DIR = "data/"

import numpy as np
import pandas as pd
import os
import pyEDM
from mpi4py import MPI

def bestE_fun(row):
    lib = row.iloc[2:]
    if np.any(pd.isna(lib)):
        Ebest = np.nan
    else:
        E_range = range(2, 24)
        results = []
        for E in E_range:
            s = pyEDM.Simplex(dataFrame=pd.DataFrame({'time': range(1, 530), 'lib': lib}), lib=(1, 264), pred=(265, 529), E=E, columns='lib', target='lib')
            metrics = pyEDM.ComputeError(s.Observations, s.Predictions)
            results.append(metrics['MAE'])
        Ebest = E_range[int(np.nanargmin(results))] if np.isfinite(results).any() else np.nan
    return_df = pd.DataFrame({'x': row.iloc[0], 'y': row.iloc[1], 'Ebest': [Ebest]})
    return return_df

def process_tile(tile):
    df_path = DATA_DIR + 'GlobalLISMap1km/Results/4Neighborhood/'
    result_path = DATA_DIR + 'GlobalLISMap1km/Results/5BestE/'
    df = pd.read_csv(f'{df_path}Neigh_{tile}.csv')
    if df.shape[1] != 531:
        raise ValueError('Expected two coordinate columns and 529 full-period observations.')
    os.makedirs(result_path, exist_ok=True)
    Complete_row = df.notna().all(axis=1)
    df_complete = df[Complete_row]
    results = []
    for row_n in range(len(df_complete)):
        row_i = df_complete.iloc[row_n, :]
        result = bestE_fun(row_i)
        results.append(result)
    total_result = pd.concat(results, ignore_index=True) if results else pd.DataFrame(columns=['x', 'y', 'Ebest'])
    tile_result_df = df.iloc[:, 0:2].copy()
    tile_result_df.columns = ["x", "y"]
    tile_result_df['Ebest'] = np.nan
    tile_result_df.loc[Complete_row, :] = total_result.values
    tile_result_df.to_csv(os.path.join(result_path, f'BestE_{tile}.csv'), index=False)
if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    STL_path = DATA_DIR + 'GlobalLISMap1km/Results/3stl/'
    file_list = [f for f in os.listdir(STL_path) if f.startswith('Residual_') and f.endswith('.csv')]
    Tile_list = sorted([int(f.split('_')[1].split('.')[0]) for f in file_list])
    all_tiles = Tile_list
    base_tiles_per_node = len(all_tiles) // size
    remaining_tiles = len(all_tiles) % size
    start_idx = rank * (base_tiles_per_node + 1) if rank < remaining_tiles else rank * base_tiles_per_node + remaining_tiles
    end_idx = (rank + 1) * (base_tiles_per_node + 1) if rank < remaining_tiles else (rank + 1) * base_tiles_per_node + remaining_tiles
    my_tiles = all_tiles[start_idx:end_idx]
    for tile in my_tiles:
        process_tile(int(tile))
