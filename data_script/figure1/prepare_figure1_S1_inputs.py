#!/usr/bin/env python3
import csv
import hashlib
from collections import Counter
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font


ORDER_LEVELS = [
    "Anguilliformes", "Aulopiformes", "Carangaria incertae sedis",
    "Carangiformes", "Chaetodontiformes", "Clupeiformes",
    "Eupercaria incertae sedis", "Gadiformes", "Gobiiformes",
    "Lutjaniformes", "Mugiliformes", "Myctophiformes", "Perciformes",
    "Pleuronectiformes", "Scombriformes", "Siluriformes", "Spariformes",
    "Tetraodontiformes", "Carcharhiniformes", "Myliobatiformes",
]


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    figure1 = Path(__file__).resolve().parent
    figure_s1 = figure1.parent / "figureS1"
    source = figure1 / "Table_S1_final_host_reads.xlsx"
    workbook = load_workbook(source, data_only=True)
    sheet = workbook["main"]
    header = [cell.value for cell in sheet[2]]
    source_rows = [dict(zip(header, row)) for row in sheet.iter_rows(min_row=3, values_only=True)
                   if row[1] is not None]
    library_ids = [str(row["lib_index"]).strip() for row in source_rows]
    if len(source_rows) != 1027 or len(set(library_ids)) != 1027:
        raise ValueError("Table S1 must contain 1027 unique library IDs.")

    order_index = {name: index + 1 for index, name in enumerate(ORDER_LEVELS)}
    libraries = []
    individuals = {}
    for row in source_rows:
        library_id = str(row["lib_index"]).strip()
        individual_id = library_id.rsplit("-", 1)[0]
        host_order = row["order / order-level lineage"].strip()
        if host_order not in order_index:
            raise ValueError(f"Add the observed host order to ORDER_LEVELS: {host_order}")
        individual = {
            "individual_id": individual_id,
            "species_abbreviation": row["species abbreviation"].strip(),
            "species": row["species"].strip(),
            "genus": row["genus"].strip(),
            "order_order_level_lineage": host_order,
            "class": row["class"].strip(),
        }
        if individual_id in individuals and individuals[individual_id] != individual:
            raise ValueError(f"Conflicting host annotations for {individual_id}")
        individuals[individual_id] = individual

        total = int(row["total_reads"])
        qc = int(row["after_qc_reads"])
        simple = int(row["remove_simple_region_reads"])
        non_rrna = int(row["Non-rRNA reads used for host mapping"])
        host = int(row["Concordantly host-mapped reads"])
        remaining = non_rrna - host
        if not (total >= qc >= simple >= non_rrna >= remaining >= 0 and host >= 0):
            raise ValueError(f"Invalid read-count progression for {library_id}")
        if total == 0 or non_rrna == 0:
            raise ValueError(f"Zero read-count denominator for {library_id}")
        supplied_remaining = row["Reads retained after host filtering"]
        if supplied_remaining is not None and int(supplied_remaining) != remaining:
            raise ValueError(f"Retained count disagrees with host counts for {library_id}")
        filtered = total - non_rrna
        qc_removed = total - simple
        rrna_filtering = simple - non_rrna
        libraries.append({
            "library_id": library_id,
            "individual_id": individual_id,
            "host_class": row["class"],
            "host_order": host_order,
            "host_family": row["family"],
            "host_genus": row["genus"],
            "host_species": row["species"],
            "species_abbreviation": row["species abbreviation"],
            "tissue_type": row["tissue_types"],
            "sample_longitude": row["sample_longitude"],
            "sample_latitude": row["sample_latitude"],
            "total_reads": total,
            "after_qc_reads": qc,
            "after_simple_filter_reads": simple,
            "non_rrna_input_reads": non_rrna,
            "concordant_host_reads": host,
            "remaining_reads": remaining,
            "qc_removed_reads": qc_removed,
            "rrna_filtering_reads": rrna_filtering,
            "qc_removed_pct_total": 100.0 * qc_removed / total,
            "rrna_filtering_pct_total": 100.0 * rrna_filtering / total,
            "qc_rrna_removed_reads": filtered,
            "qc_rrna_removed_pct_total": 100.0 * filtered / total,
            "host_mapped_pct_total": 100.0 * host / total,
            "remaining_pct_total": 100.0 * remaining / total,
            "host_mapped_pct_non_rrna": 100.0 * host / non_rrna,
            "remaining_pct_non_rrna": 100.0 * remaining / non_rrna,
            "order_position": order_index[host_order],
        })
    if len(individuals) != 134:
        raise ValueError(f"Expected 134 fish individuals, found {len(individuals)}")
    libraries.sort(key=lambda row: (row["order_position"], -row["remaining_reads"],
                                   row["library_id"]))
    for position, row in enumerate(libraries, 1):
        row["library_position"] = position

    species_counts = Counter((row["host_order"], row["host_species"],
                              row["species_abbreviation"]) for row in libraries)
    species_rows = []
    for position, key in enumerate(sorted(species_counts, key=lambda k: (order_index[k[0]], k[2])), 1):
        species_rows.append({"host_order": key[0], "species": key[1],
                             "species_abbreviation": key[2],
                             "library_count": species_counts[key],
                             "order_position": order_index[key[0]],
                             "species_position": position})
    if len(species_rows) != 51:
        raise ValueError(f"Expected 51 host species, found {len(species_rows)}")

    genus_counts = Counter((row["order_order_level_lineage"], row["genus"])
                           for row in individuals.values())
    order_totals = Counter(row["order_order_level_lineage"] for row in individuals.values())
    genera = sorted({row["genus"] for row in individuals.values()})
    composition = Workbook()
    main_sheet = composition.active
    main_sheet.title = "main"
    main_sheet.append(["order_order_level_lineage"] + genera)
    for order in sorted(order_totals):
        main_sheet.append([order] + [genus_counts[order, genus] / order_totals[order]
                                    for genus in genera])
    counts_sheet = composition.create_sheet("individual_counts")
    counts_sheet.append(["order_order_level_lineage", "host_genus", "individual_count",
                         "within_order_proportion"])
    for (order, genus), count in sorted(genus_counts.items()):
        counts_sheet.append([order, genus, count, count / order_totals[order]])
    mapping_sheet = composition.create_sheet("host_mapping")
    mapping_sheet.append(list(next(iter(individuals.values()))))
    for individual_id in sorted(individuals):
        mapping_sheet.append(list(individuals[individual_id].values()))
    sources = composition.create_sheet("sources")
    sources.append(["source", "sha256", "counting_unit"])
    sources.append(["../figure1/Table_S1_final_host_reads.xlsx",
                    hashlib.sha256(source.read_bytes()).hexdigest(), "one fish individual"])

    (figure1 / "figure1c").mkdir(exist_ok=True)
    figure_s1.mkdir(exist_ok=True)
    write_tsv(figure1 / "figure1c" / "figure1C_species_library_count.tsv", species_rows)
    write_tsv(figure_s1 / "FigureS1C_library_read_qc.tsv", libraries)
    for sheet in composition:
        for row in sheet:
            for cell in row:
                cell.font = Font(color="000000")
    composition.save(figure_s1 / "FigureS1B_host_genus_composition_input.xlsx")
    print(f"Prepared {len(libraries)} libraries, {len(species_rows)} species, "
          f"{len(individuals)} individuals and {len(order_totals)} host orders.")
    print("Figure S1C partitions total reads into QC-filtered, further-filtering/rRNA, host-mapped and remaining reads.")


if __name__ == "__main__":
    main()
