import argparse

import cell2location
import matplotlib.pyplot as plt
import scanpy as sc
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes


def arguments():
    """Parse arguments."""
    parser = argparse.ArgumentParser(description="Cell2location arguments.")
    parser.add_argument("-s", "--sample", type=str, help="Define the sample (dataset)")
    args = parser.parse_args()

    return args


def main(sample, ref_path, annotation, model_path):
    adata_ref = sc.read_h5ad(ref_path)
    adata_ref.var["SYMBOL"] = adata_ref.var.index
    adata_ref.var.set_index("gene_ids", drop=True, inplace=True)
    adata_ref.X = adata_ref.layers["counts"].copy()

    selected = filter_genes(
        adata_ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    adata_ref = adata_ref[:, selected].copy()

    cell2location.models.RegressionModel.setup_anndata(
        adata=adata_ref,
        # 10X reaction / sample / batch
        batch_key="sample",
        # cell type, covariate used for constructing signatures
        labels_key=annotation,
    )

    mod = RegressionModel(adata_ref)
    mod.train(max_epochs=500, accelerator="gpu")

    mod.save(
        model_path,
        save_anndata=True,
    )

    print(f"Working on model")
    mod.plot_history(100)
    plt.savefig(f"{model_path}/training_curve.png")
    plt.close()

    adata_ref = mod.export_posterior(
        adata_ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
    )
    mod.plot_QC()
    plt.savefig(f"{model_path}/model_QC.png")
    plt.close()

    # export estimated expression in each cluster
    if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in adata_ref.uns["mod"]["factor_names"]
            ]
        ].copy()
    else:
        inf_aver = adata_ref.var[
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in adata_ref.uns["mod"]["factor_names"]
            ]
        ].copy()
    inf_aver.columns = adata_ref.uns["mod"]["factor_names"]

    print(inf_aver.iloc[0:5, 0:5])

    inf_aver.to_csv(
        f"chromothripsis/j462r/spatial_transcriptomics/cell2location/inf_aver_{sample}_merged.tsv",
        sep="\t",
    )


if __name__ == "__main__":
    args = arguments()
    sample = args.sample

    if sample == "chen":
        ref_path = "chromothripsis/j462r/spatial_transcriptomics/scRNA_references/chen_nature/results/chen_annotated.h5ad"
        annotation = "broad_celltypes"
        model_path = "Chen_reference_model_celltype"
    elif sample == "gouin":
        ref_path = "chromothripsis/j462r/spatial_transcriptomics/scRNA_references/gouin_nature/merged.h5ad"
        annotation = "subtype"
        model_path = "Chen_reference_model_celltype"

    main(sample=sample, ref_path=ref_path, annotation=annotation, model_path=model_path)
