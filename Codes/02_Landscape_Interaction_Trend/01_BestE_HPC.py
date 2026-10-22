"""Estimate window-specific embedding dimensions in Python on an HPC system using MPI."""

import os

import numpy as np
import pandas as pd
import pyEDM
from mpi4py import MPI

DATA_DIR = "data/"


def bestE_fun(row):
    series = row.iloc[2:].to_numpy()
    best_e = np.nan
    if not np.any(pd.isna(series)):
        dimensions = range(2, int(np.sqrt(len(series))) + 1)
        frame = pd.DataFrame({"time": range(1, len(series) + 1), "lib": series})
        split = len(series) // 2
        errors = []
        for embedding in dimensions:
            prediction = pyEDM.Simplex(
                dataFrame=frame,
                lib=(1, split),
                pred=(split + 1, len(series)),
                E=embedding,
                columns="lib",
                target="lib",
            )
            error = pyEDM.ComputeError(prediction.Observations, prediction.Predictions)
            errors.append(error["MAE"])
        if np.any(np.isfinite(errors)):
            best_e = dimensions[int(np.nanargmin(errors))]
    return {"x": row.iloc[0], "y": row.iloc[1], "Ebest": best_e}


def process_tile(tile):
    input_dir = DATA_DIR + "GlobalLISMap1km/Results/4Neighborhood/"
    result_dir = DATA_DIR + "GlobalLISMap1km_Trend_py/Results/BestE/"
    os.makedirs(result_dir, exist_ok=True)
    total = pd.read_csv(f"{input_dir}Neigh_{tile}.csv")
    complete = total.notna().all(axis=1)
    for year in range(19):
        columns = [0, 1] + list(range(2 + year * 23, 2 + (year + 5) * 23))
        window = total.iloc[:, columns]
        output = window.iloc[:, :2].copy()
        output.columns = ["x", "y"]
        output["Ebest"] = np.nan
        if complete.any():
            records = [bestE_fun(row) for _, row in window.loc[complete].iterrows()]
            output.loc[complete, :] = pd.DataFrame(records).to_numpy()
        output.to_csv(f"{result_dir}BestE_{tile}_year{year}.csv", index=False)


def main():
    comm = MPI.COMM_WORLD
    residual_dir = DATA_DIR + "GlobalLISMap1km/Results/3stl/"
    tiles = sorted(
        int(name.split("_")[1].split(".")[0])
        for name in os.listdir(residual_dir)
        if name.startswith("Residual_") and name.endswith(".csv")
    )
    for tile in tiles[comm.Get_rank() :: comm.Get_size()]:
        process_tile(tile)


if __name__ == "__main__":
    main()
