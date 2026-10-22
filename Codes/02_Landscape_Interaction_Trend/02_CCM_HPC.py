"""Estimate rolling-window landscape interactions in Python on an HPC system using MPI."""

import os

import numpy as np
import pandas as pd
from mpi4py import MPI
from pyEDM import CCM
from scipy.stats import kendalltau, norm

DATA_DIR = "data/"
DIRECTIONS = ("east", "north", "south", "west")


def lagged_ccm(focal, neighbor, embedding):
    focal = np.asarray(focal)
    neighbor = np.asarray(neighbor)
    if pd.isna(embedding) or pd.isna(focal).any() or pd.isna(neighbor).any():
        return np.nan, np.nan, np.nan
    estimates = []
    for lag in range(-12, 1):
        series = pd.DataFrame(
            {
                "time": range(1, len(focal) + lag + 1),
                "lib": focal[-lag:],
                "target": neighbor[: len(neighbor) + lag],
            }
        )
        sizes = sorted(
            size
            for size in {11, 15, 20, 30, 40, 50, 70, 90, 110, len(series)}
            if int(embedding) < size <= len(series)
        )
        if len(sizes) < 3:
            continue
        skill = CCM(
            dataFrame=series,
            E=int(embedding),
            Tp=0,
            columns="lib",
            target="target",
            libSizes=" ".join(map(str, sizes)),
            random=False,
            showPlot=False,
        )
        rho = skill["lib:target"].to_numpy()
        valid = np.isfinite(rho)
        rho = rho[valid]
        if len(rho) < 3:
            continue
        tau, probability = kendalltau(np.arange(len(rho)), rho)
        n_library = skill.loc[valid, "LibSize"].max()
        bounds = np.clip([rho[0], rho[-1]], -1 + 1e-10, 1 - 1e-10)
        z = abs(np.diff(np.arctanh(bounds))[0]) * np.sqrt((n_library - 3) / 2)
        significant = int(tau > 0 and probability < 0.05 and norm.sf(z) < 0.05)
        estimates.append((np.max(rho), lag, significant))
    return max(estimates, key=lambda result: result[0]) if estimates else (np.nan, np.nan, 0)


def process_tile(tile):
    input_dir = DATA_DIR + "GlobalLISMap1km/Results/4Neighborhood/"
    best_e_dir = DATA_DIR + "GlobalLISMap1km_Trend_py/Results/BestE/"
    result_dir = DATA_DIR + "GlobalLISMap1km_Trend_py/Results/CCM/"
    os.makedirs(result_dir, exist_ok=True)
    total = pd.read_csv(f"{input_dir}Neigh_{tile}.csv")
    neighbors = {
        direction: pd.read_csv(f"{input_dir}Neigh_{tile}_{suffix}.csv")
        for direction, suffix in zip(DIRECTIONS, ("E", "N", "S", "W"))
    }
    complete = total.notna().all(axis=1)
    for year in range(19):
        columns = [0, 1] + list(range(2 + year * 23, 2 + (year + 5) * 23))
        window = total.iloc[:, columns]
        best_e = pd.read_csv(f"{best_e_dir}BestE_{tile}_year{year}.csv")
        output = best_e.copy()
        output.columns = ["x", "y", "Ebest"]
        for direction in DIRECTIONS:
            output[[direction, f"{direction}_tp", f"{direction}_sig"]] = np.nan
        for index in window.index[complete]:
            focal = window.loc[index].iloc[2:].to_numpy()
            embedding = best_e.loc[index, "Ebest"]
            for direction in DIRECTIONS:
                neighbor = neighbors[direction].iloc[index, columns].iloc[2:].to_numpy()
                estimate = lagged_ccm(focal, neighbor, embedding)
                output.loc[index, [direction, f"{direction}_tp", f"{direction}_sig"]] = estimate
        output.to_csv(f"{result_dir}CCM_{tile}_year{year}.csv", index=False)


def main():
    comm = MPI.COMM_WORLD
    tiles = pd.read_csv(DATA_DIR + "GlobalLISMap1km/Results/Valid_Tiles.csv")
    all_tiles = tiles.iloc[:, 0].to_numpy()
    for tile in all_tiles[comm.Get_rank() :: comm.Get_size()]:
        process_tile(int(tile))


if __name__ == "__main__":
    main()
