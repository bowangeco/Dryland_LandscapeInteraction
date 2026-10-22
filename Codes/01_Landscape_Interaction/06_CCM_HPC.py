"""Full-period CCM: author-supplied Python analysis executed on an HPC system using MPI."""

DATA_DIR = "data/"

import numpy as np
import pandas as pd
from scipy.stats import norm
import os
from pyEDM import CCM
import pymannkendall as mk
from mpi4py import MPI

def ccm_fun(residual_ts, Ebest, east_ts, north_ts, south_ts, west_ts, x, y):
    lib = residual_ts
    if np.any(pd.isna(lib)) or pd.isna(Ebest):
        ccm_result = pd.DataFrame({'east': [np.nan], 'east_tp': [np.nan], 'east_sig': [np.nan], 'north': [np.nan], 'north_tp': [np.nan], 'north_sig': [np.nan], 'south': [np.nan], 'south_tp': [np.nan], 'south_sig': [np.nan], 'west': [np.nan], 'west_tp': [np.nan], 'west_sig': [np.nan]})
    else:
        test_df = pd.DataFrame({'time': range(1, 530), 'lib': lib, 'east': east_ts, 'north': north_ts, 'south': south_ts, 'west': west_ts})
        ccm_results = []
        for direction in ['east', 'north', 'south', 'west']:
            if test_df[direction].isna().any():
                ccm_results.append(pd.DataFrame({direction: [np.nan],
                                                f'{direction}_tp': [np.nan],
                                                f'{direction}_sig': [np.nan]}))
                continue
            lag_seq = range(-12, 1)
            lags_rho = []
            for lag in lag_seq:
                lag_df = pd.DataFrame({'time': range(1, 530 + lag), 'lib': test_df['lib'].iloc[-lag:529].values, 'target': test_df[direction].iloc[0:529 + lag].values})
                cmap = CCM(dataFrame=lag_df, E=int(Ebest), Tp=0, columns='lib', target='target', random=False, libSizes=str(len(lag_df)), showPlot=False, sample=0)
                rho = cmap.loc[0, 'lib:target']
                lags_rho.append(rho)
            rho_max = max(lags_rho)
            tp_max = lag_seq[np.argmax(lags_rho)]
            libSizes = [24, 30, 40, 50, 70] + list(range(100, len(test_df), 50)) + [len(test_df)]
            signif_ccm = CCM(dataFrame=test_df, E=int(Ebest), columns='lib', target=direction, libSizes=" ".join(map(str, libSizes)), random=False, showPlot=False)
            kend = mk.original_test(signif_ccm.iloc[:, 1])
            rho_Lmax = signif_ccm.iloc[-1, 1]
            rho_Lmin = signif_ccm.iloc[0, 1]
            N = max(signif_ccm['LibSize'])
            z = abs(0.5 * (np.log((1 + rho_Lmax) / (1 - rho_Lmax)) - np.log((1 + rho_Lmin) / (1 - rho_Lmin))) * (2 / (N - 3)) ** (-0.5))
            z_p = 1 - norm.cdf(z)
            significance = int(kend.h > 0 and kend.p < 0.05 and (z_p < 0.05))
            r_df = pd.DataFrame({direction: [rho_max], f'{direction}_tp': [tp_max], f'{direction}_sig': [significance]})
            ccm_results.append(r_df)
        ccm_result = pd.concat(ccm_results, axis=1)
    return_df = pd.concat([pd.DataFrame({'x': x, 'y': y, 'Ebest': [Ebest]}), ccm_result], axis=1)
    return return_df

def process_tile(tile):
    residual_path = DATA_DIR + 'GlobalLISMap1km/Results/4Neighborhood/'
    result_path = DATA_DIR + 'GlobalLISMap1km/Results/6CCM/'
    terminal_result_path = DATA_DIR + 'GlobalLISMap1km/Results/6CCM_terminal/'
    residual_df = pd.read_csv(f'{residual_path}Neigh_{tile}.csv')
    east_df = pd.read_csv(f'{residual_path}Neigh_{tile}_E.csv')
    north_df = pd.read_csv(f'{residual_path}Neigh_{tile}_N.csv')
    south_df = pd.read_csv(f'{residual_path}Neigh_{tile}_S.csv')
    west_df = pd.read_csv(f'{residual_path}Neigh_{tile}_W.csv')
    Ebest_df = pd.read_csv(f'{DATA_DIR}GlobalLISMap1km/Results/5BestE/BestE_{tile}.csv')
    if residual_df.shape[1] != 531:
        raise ValueError('Expected two coordinate columns and 529 full-period observations.')
    for frame in (east_df, north_df, south_df, west_df):
        if frame.shape != residual_df.shape or not np.allclose(
                frame.iloc[:, :2].to_numpy(dtype=float),
                residual_df.iloc[:, :2].to_numpy(dtype=float),
                rtol=0, atol=1e-8, equal_nan=True):
            raise ValueError('Neighbor rows do not match the focal coordinates or observation count.')
    if len(Ebest_df) != len(residual_df) or not np.allclose(
            Ebest_df.iloc[:, :2].to_numpy(dtype=float),
            residual_df.iloc[:, :2].to_numpy(dtype=float),
            rtol=0, atol=1e-8, equal_nan=True):
        raise ValueError('BestE rows do not match the focal coordinates.')
    os.makedirs(result_path, exist_ok=True)
    os.makedirs(terminal_result_path, exist_ok=True)
    Complete_row = residual_df.notna().all(axis=1)
    residual_df_complete = residual_df[Complete_row]
    east_df_complete = east_df[Complete_row]
    north_df_complete = north_df[Complete_row]
    south_df_complete = south_df[Complete_row]
    west_df_complete = west_df[Complete_row]
    Ebest_df_complete = Ebest_df[Complete_row]
    results = []
    for row_n in range(len(residual_df_complete)):
        x = residual_df_complete.iloc[row_n, 0]
        y = residual_df_complete.iloc[row_n, 1]
        residual_ts = residual_df_complete.iloc[row_n, 2:]
        Ebest = Ebest_df_complete.iloc[row_n, 2]
        east_ts = east_df_complete.iloc[row_n, 2:]
        north_ts = north_df_complete.iloc[row_n, 2:]
        south_ts = south_df_complete.iloc[row_n, 2:]
        west_ts = west_df_complete.iloc[row_n, 2:]
        result = ccm_fun(residual_ts, Ebest, east_ts, north_ts, south_ts, west_ts, x, y)
        results.append(result)
    total_result = pd.concat(results, ignore_index=True) if results else pd.DataFrame(columns=['x', 'y', 'Ebest', 'east', 'east_tp', 'east_sig', 'north', 'north_tp', 'north_sig', 'south', 'south_tp', 'south_sig', 'west', 'west_tp', 'west_sig'])
    total_result.to_csv(f'{terminal_result_path}OneTile{tile}_CCM_.csv', index=False)
    new_cols = ['east', 'east_tp', 'east_sig', 'north', 'north_tp', 'north_sig', 'south', 'south_tp', 'south_sig', 'west', 'west_tp', 'west_sig']
    ccm_df_tile = Ebest_df.copy()
    ccm_df_tile.columns = ["x", "y", "Ebest"]
    ccm_df_tile[new_cols] = pd.DataFrame([[np.nan] * len(new_cols)] * len(ccm_df_tile), index=ccm_df_tile.index)
    ccm_df_tile.loc[Complete_row, new_cols] = total_result.loc[:, new_cols].values
    ccm_df_tile.to_csv(os.path.join(result_path, f'CCM_{tile}.csv'), index=False)
if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    Tiles_df = pd.read_csv(f'{DATA_DIR}GlobalLISMap1km/Results/Valid_Tiles.csv')
    all_tiles = Tiles_df.iloc[:, 0].values
    base_tiles_per_node = len(all_tiles) // size
    remaining_tiles = len(all_tiles) % size
    start_idx = rank * (base_tiles_per_node + 1) if rank < remaining_tiles else rank * base_tiles_per_node + remaining_tiles
    end_idx = (rank + 1) * (base_tiles_per_node + 1) if rank < remaining_tiles else (rank + 1) * base_tiles_per_node + remaining_tiles
    my_tiles = all_tiles[start_idx:end_idx]
    for tile in my_tiles:
        process_tile(int(tile))
