import argparse
import logging

import cell2location
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def arguments():
    """Parse arguments."""
    parser = argparse.ArgumentParser(description="Cell2location arguments.")
    parser.add_argument("-a", "--alpha", type=int, help="Alpha", default=20)
    parser.add_argument(
        "-n", "--cell_number", type=int, help="Cell number per spot", default=30
    )
    parser.add_argument(
        "-c",
        "--celltypes",
        type=str,
        help="Cell type annotation type ('subtype' or 'celltype')",
        default="celltype",
    )
    args = parser.parse_args()

    return args


def prepare_data(adata_vis_path, inf_aver_path):
    """Prepare data for mapping."""
    logging.info("Reading and preparing data.")

    adata_vis = sc.read_h5ad(adata_vis_path)

    adata_vis.var["SYMBOL"] = adata_vis.var_names
    adata_vis.var.set_index("gene_ids", drop=True, inplace=True)

    inf_aver = pd.read_csv(
        inf_aver_path,
        sep="\t",
        index_col="gene_ids",
    )

    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    return adata_vis, inf_aver


def run_mapping(adata_vis, inf_aver, N_cells_per_location, detection_alpha):
    """Map cell types to visium data."""
    logging.info("Running spatial mapping")
    model_name = f"visium_model_alpha{alpha}_N{cell_number}_{celltype}"

    cell2location.models.Cell2location.setup_anndata(
        adata=adata_vis, batch_key="sample", layer="counts"
    )
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=cell_number,
        detection_alpha=alpha,
    )

    mod.train(max_epochs=30000, batch_size=None, train_size=1, accelerator="gpu")
    mod.save(
        f"chromothripsis/j462r/spatial_transcriptomics/cell2location/{model_name}",
        save_anndata=True,
    )

    mod.plot_history(1000)
    plt.savefig(
        f"chromothripsis/j462r/spatial_transcriptomics/cell2location/{model_name}/training_performance.png"
    )
    plt.close()

    return mod, adata_vis, model_name


def export_posterior(adata_vis, mod, model_name):
    """Export posterior distributions of cell abundances."""
    logging.info("Exporting posteriors")

    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={
            "num_samples": 1000,
            "batch_size": mod.adata.n_obs,
        },
    )

    mod.plot_QC()
    plt.savefig(
        f"chromothripsis/j462r/spatial_transcriptomics/cell2location/{model_name}/model_QC.png"
    )
    plt.close()

    adata_vis.obs[adata_vis.uns["mod"]["factor_names"]] = adata_vis.obsm[
        "q05_cell_abundance_w_sf"
    ]
    logging.info("Saving final adata")
    adata_vis.write_h5ad(
        f"chromothripsis/j462r/spatial_transcriptomics/cell2location/{model_name}/posteriors_adata.h5ad"
    )


if __name__ == "__main__":
    args = arguments()
    alpha = args.alpha
    celltype = args.celltypes
    cell_number = args.cell_number

    inf_aver_path = f"chromothripsis/j462r/spatial_transcriptomics/cell2location/inf_aver_{celltype}.tsv"
    adata_vis_path = "chromothripsis/j462r/spatial_transcriptomics/scripts/h5ad_objects/merged_samples.h5ad"

    adata_vis, inf_aver = prepare_data(
        adata_vis_path=adata_vis_path, inf_aver_path=inf_aver_path
    )

    mod, adata_vis, model_name = run_mapping(
        adata_vis=adata_vis,
        inf_aver=inf_aver,
        N_cells_per_location=cell_number,
        detection_alpha=alpha,
    )

    export_posterior(adata_vis=adata_vis, mod=mod, model_name=model_name)
