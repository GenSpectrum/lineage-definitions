import argparse
import os
from pathlib import Path
import pandas as pd
import requests


def download_file(url, save_path):
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()

        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        Path(save_path).write_bytes(response.content)

        print(f"File downloaded successfully and saved to '{save_path}'")

    except requests.exceptions.RequestException as e:
        msg = f"Failed to download path: {e}"
        print(msg)
        raise RuntimeError(msg) from e


def load_lineage_files(url: str, output_dir: str) -> dict[str, dict[str, list[str]]]:
    save_path = os.path.join(output_dir, "clades.tsv")
    download_file(url, save_path)
    df = pd.read_csv(save_path, sep="\t")
    unique_clades = df["clade"].dropna().unique()

    print(unique_clades)

    silo_lineages: dict = {}
    silo_lineages["unassigned"] = {"aliases": [], "parents": []}
    silo_lineages["None"] = {"aliases": [], "parents": []}
    for clade in unique_clades:
        genotype = clade.split("_")[0]
        if genotype not in silo_lineages:
            silo_lineages[genotype] = {
                "aliases": [],
                "parents": [],
            }
        if len(clade.split("_")) > 1:
            major_lineage = clade.split(".")[0]
            if major_lineage not in silo_lineages:
                silo_lineages[major_lineage] = {
                    "aliases": [],
                    "parents": [genotype],
                }
            if len(clade.split(".")) > 1:
                minor_lineage = clade
                if minor_lineage not in silo_lineages:
                    silo_lineages[minor_lineage] = {
                        "aliases": [],
                        "parents": [major_lineage],
                    }
    return silo_lineages


def parse_label(clade: str) -> tuple[str, str, str]:
    """Parse clade to remove unwanted characters."""
    genotype = clade.split("_")[0]
    major_lineage = ""
    minor_lineage = ""
    if len(clade.split("_")) > 1:
        major_lineage = clade.split(".")[0]
        if len(clade.split(".")) > 1:
            minor_lineage = clade
    return (genotype, major_lineage, minor_lineage)


def generate_hierarchy_file(
    lineage_file_url: str, output_dir: str, dataset_tag: str | None = None
) -> None:
    """Generate hierarchy file from github repo."""

    print("Loading lineage files...")
    lineages = load_lineage_files(lineage_file_url, output_dir)

    # Create output directory structure
    if dataset_tag:
        subtype_dir = os.path.join(output_dir, dataset_tag)
    else:
        subtype_dir = Path(output_dir)
    os.makedirs(subtype_dir, exist_ok=True)

    # Write hierarchy file
    output_file = os.path.join(subtype_dir, "lineages.yaml")

    with open(output_file, "w") as f:
        # Convert OrderedDict to regular dict for clean YAML output
        # but maintain the topological order by iterating through our sorted structure
        for name in sorted(lineages.keys(), key=parse_label):
            f.write(f"{name}:\n")
            f.write(f"  aliases: {lineages[name]['aliases']}\n")
            if lineages[name]["parents"]:
                f.write(f"  parents:\n")
                for parent in lineages[name]["parents"]:
                    f.write(f"  - {parent}\n")
            else:
                f.write(f"  parents: []\n")

        print(f"Generated hierarchy file: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Dengue hierarchy definition YAML files"
    )

    parser.add_argument(
        "--dataset-tag",
        help="Dataset tag to organize files in subfolders (e.g., nextclade dataset tag)",
    )

    args = parser.parse_args()

    url = "https://raw.githubusercontent.com/V-GEN-Lab/nextclade-datasets-workflow/main/denv/resources/"

    for subtype in ["denv1", "denv2", "denv3", "denv4"]:
        print(f"🚀 Starting {subtype} hierarchy generation...")
        generate_hierarchy_file(
            url + subtype + "/clades.tsv", subtype, args.dataset_tag
        )
        print(f"✅ {subtype} hierarchy generation completed successfully!")


if __name__ == "__main__":
    main()
