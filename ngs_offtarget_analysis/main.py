import os
import time
import csv
from pathlib import Path
from Bio.SeqUtils import GC
from scipy.stats import chisquare

from utils.io import read_info_file, read_ngs_output, get_genome
from utils.methods import get_target_site, add_cds_to_reads
from utils.constants import GENOME_BUFFER, AT_PCT_TEST_SIZE
from analysis_methods.at_enrichment import get_at_test_data, get_at_likelihood
from analysis_methods.alternate_targets import get_likely_target
from analysis_methods.sample_correlations import correlation_scores

from parameters import (
    info_file,
    samples,
    output_subdirectory,
    AT_enrichment_sizes,
    TARGET_WINDOW,
    INTEGRATION_SITE_DISTANCE,
    analysis_type,
    CORRELATIONS_BIN_SIZE,
)


### The code is here ###
def main():
    start_time = time.perf_counter()

    if analysis_type not in ["both", "correlation", "offtarget"]:
        raise Exception(
            "The 'analysis_type' parameters must be 'offtarget', 'correlation', or 'both'"
        )
    print("Starting the run...")
    outputs_directory = os.path.join(
        Path(__file__).parent.absolute(), "outputs", output_subdirectory
    )
    os.makedirs(outputs_directory, exist_ok=True)

    # Load E.Coli essential gene list
    essential_genes = []
    with open("Ecoli_essential_genes.csv", "r", encoding="utf-8-sig") as infile:
        reader = csv.DictReader(infile)
        essential_genes = [r["Gene"] for r in reader if r["Essential"] == "TRUE"]

    all_samples = read_info_file(info_file)
    samples_to_run = samples
    if not samples_to_run:
        samples_to_run = list(all_samples.keys())

    if analysis_type in ["both", "offtarget"]:
        results = {}

        for sample in samples_to_run:
            os.makedirs(os.path.join(outputs_directory, sample), exist_ok=True)
            print(f"Analyzing AT richness and off-target binding for {sample}...")
            sample_results = {"Sample": sample}

            sample_info = all_samples[sample]
            reads = read_ngs_output(sample)
            genome = get_genome(sample_info)

            target_site = get_target_site(
                sample_info, genome, INTEGRATION_SITE_DISTANCE
            )

            offtarget_reads = [
                r
                for r in reads
                if abs(int(r["position"]) - target_site) >= TARGET_WINDOW
            ]

            # buffer the genome for sequences at start/end
            buffered_genome = (
                genome[-GENOME_BUFFER:].seq + genome.seq + genome[:GENOME_BUFFER].seq
            ).upper()

            # First do AT enrichment tests
            for at_window_size in AT_enrichment_sizes:
                # The name of the key to save this data as
                at_pct_key = f"at_pct_{at_window_size}bp"
                # Get test data
                offtarget_reads_ct = sum([int(r["reads"]) for r in offtarget_reads])
                at_sampling = get_at_test_data(genome, at_window_size)
                with open(
                    os.path.join(
                        outputs_directory,
                        sample,
                        f"{at_window_size}bp_random_sampling.csv",
                    ),
                    "w",
                ) as outfile:
                    writer = csv.DictWriter(
                        outfile, fieldnames=["at_pct", "count", "scaled_to_sample_ct"]
                    )
                    writer.writeheader()
                    for (at_pct, read_ct) in at_sampling.items():
                        scaled_ct = round(
                            read_ct * (offtarget_reads_ct / AT_PCT_TEST_SIZE)
                        )
                        writer.writerow(
                            {
                                "at_pct": at_pct,
                                "count": read_ct,
                                "scaled_to_sample_ct": scaled_ct,
                            }
                        )

                # Calculate AT pct and relative odds at each position
                for row in reads:
                    row[at_pct_key] = None
                    row[f"{at_pct_key}_odds"] = None
                for row in offtarget_reads:
                    position = int(row["position"])
                    window_up = buffered_genome[
                        position
                        + GENOME_BUFFER
                        - at_window_size : position
                        + GENOME_BUFFER
                    ]
                    window_down = buffered_genome[
                        position
                        + GENOME_BUFFER : position
                        + GENOME_BUFFER
                        + at_window_size
                    ]
                    row[at_pct_key] = round(
                        100 - min([GC(window_up), GC(window_down)]), 2
                    )
                    row[f"{at_pct_key}_odds"] = get_at_likelihood(
                        row[at_pct_key], at_sampling
                    )

                # Get chisquare for the actual vs random AT-distribution counts
                possible_at_pcts = sorted(
                    set(
                        [round(row[at_pct_key], 2) for row in offtarget_reads]
                        + list(at_sampling.keys())
                    )
                )
                expected_dist = []
                actual_dist = []
                for possible_pct in possible_at_pcts:
                    read_counts_at_pct = sum(
                        [
                            int(r["reads"])
                            for r in offtarget_reads
                            if r[at_pct_key] == possible_pct
                        ]
                    )
                    expected_dist += [
                        round(
                            at_sampling.get(possible_pct, 0)
                            * (offtarget_reads_ct / AT_PCT_TEST_SIZE),
                            8,
                        )
                    ]
                    actual_dist += [read_counts_at_pct]
                chisquare_probability = chisquare(actual_dist, expected_dist)
                sample_results["chisquare_statistic_at_pct_distribution"] = round(
                    chisquare_probability.statistic, 6
                )
                sample_results["chisquare_pval_at_pct_distribution"] = round(
                    chisquare_probability.pvalue, 6
                )

            # Evaluate off-target Cascade binding sites
            for row in reads:
                position = int(row["position"])
                if abs(position - target_site) < TARGET_WINDOW:
                    row["offtarget_seq"] = ""
                    row["best_offtarget_seed_match"] = "In target window"
                    row["best_offtarget_total_match"] = "In target window"
                else:
                    tgt = get_likely_target(
                        position,
                        buffered_genome,
                        sample_info["Spacer"],
                        INTEGRATION_SITE_DISTANCE,
                    )
                    row["offtarget_seq"] = tgt["seq"]
                    row["best_offtarget_seed_match"] = tgt["seed_match"]
                    row["best_offtarget_total_match"] = tgt["total_match"]

            reads = add_cds_to_reads(reads, genome, essential_genes)

            # Write the position-level analyses to an output file
            sample_output_path = Path(
                os.path.join(outputs_directory, sample, f"offtargets_data.csv")
            )
            os.makedirs(sample_output_path.parent, exist_ok=True)
            with open(sample_output_path, "w") as sample_out:
                fieldnames = reads[0].keys()
                writer = csv.DictWriter(sample_out, fieldnames=fieldnames)
                writer.writeheader()
                for row in reads:
                    writer.writerow(row)

            # Calculate the sample-level AT enrichment:
            # A weighted-average of the AT-richness odds
            total_reads = sum([int(row["reads"]) for row in reads])
            sample_results["total_reads"] = total_reads
            # Only use off-target and non-T7 RNAP reads
            at_enrichment_reads = 0

            # Count the reads in essential genes
            essential_gene_reads = 0
            for row in reads:
                if "In target window" in str(row["best_offtarget_seed_match"]):
                    row["skip_AT_odds"] = True
                    continue
                if "CDS" in row and "T7 DNA-directed RNA polymerase" in row["CDS"]:
                    row["skip_AT_odds"] = True
                    continue
                row["skip_AT_odds"] = False
                at_enrichment_reads += int(row["reads"])
                if "Essential" in row and row["Essential"]:
                    essential_gene_reads += int(row["reads"])
            sample_results["reads_ct_for_at_enrichment"] = at_enrichment_reads
            sample_results["essential_offtarget_reads"] = essential_gene_reads
            sample_results["essential_offtarget_reads_pct"] = (
                round(essential_gene_reads / at_enrichment_reads, 6) * 100
            )

            for at_window_size in AT_enrichment_sizes:
                overall_odds = 0
                for row in reads:
                    if not row["skip_AT_odds"]:
                        overall_odds += row[f"at_pct_{at_window_size}bp_odds"] * int(
                            row["reads"]
                        )
                overall_odds = round(overall_odds / at_enrichment_reads, 4)
                sample_results[
                    f"at_weighted_avg_enrichment_{at_window_size}bp"
                ] = overall_odds

            # Calculate chisquare probably for insertions into essential genes
            T7RNAP_length = None
            essential_gene_length = 0
            for feature in genome.features:
                if feature.type != "CDS":
                    continue
                if "T7 DNA-directed RNA polymerase" in "".join(
                    feature.qualifiers["label"]
                ):
                    T7RNAP_length = feature.location.end - feature.location.start
                if feature.qualifiers["label"][0].split()[1] == "CDS":
                    gene_id = feature.qualifiers["label"][0].split()[0]
                    if gene_id and gene_id in essential_genes:
                        essential_gene_length += (
                            feature.location.end - feature.location.start
                        )

            if T7RNAP_length:
                relevant_length = len(genome.seq) - 2 * TARGET_WINDOW - T7RNAP_length
            else:
                print("No T7RNAP gene found in the provided genbank file")
                relevant_length = len(genome.seq) - 2 * TARGET_WINDOW
            essential_gene_relevant_pct = essential_gene_length / relevant_length
            expected_dist = [
                essential_gene_relevant_pct * at_enrichment_reads,
                (1 - essential_gene_relevant_pct) * at_enrichment_reads,
            ]
            actual_dist = [
                essential_gene_reads,
                at_enrichment_reads - essential_gene_reads,
            ]
            chisquare_probability = chisquare(actual_dist, expected_dist)
            sample_results["chisquare_statistic_essentialInsertions"] = round(
                chisquare_probability.statistic, 6
            )
            sample_results["chisquare_pval_essentialInsertions"] = round(
                chisquare_probability.pvalue, 6
            )

            # Persist the sample-level results
            results[sample] = sample_results
            print(
                f"Finished with sample {sample}: {time.perf_counter()-start_time} seconds"
            )

        # Output the sample-level results
        with open(
            os.path.join(outputs_directory, "summary_results.csv"), "w"
        ) as outfile:
            fieldnames = results[samples_to_run[0]].keys()
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            for (sample, summary_results) in results.items():
                writer.writerow(summary_results)

    # Do correlation analysis
    if analysis_type in ["correlation", "both"]:
        print("Starting the correlation analysis...")
        all_samples_reads = {}
        # Load all the samples
        genome = None
        for sample in samples_to_run:
            sample_info = all_samples[sample]
            reads = read_ngs_output(sample)
            genome = get_genome(sample_info)

            reads = add_cds_to_reads(reads, genome, essential_genes)

            # Remove reads that are on-target or in T7RNAP
            target_site = get_target_site(
                sample_info, genome, INTEGRATION_SITE_DISTANCE
            )
            for (idx, row) in enumerate(reads):
                position = int(row["position"])
                if abs(position - target_site) < TARGET_WINDOW or (
                    "CDS" in row and "T7 DNA-directed RNA polymerase" in row["CDS"]
                ):
                    reads[idx] = None
            reads = [r for r in reads if r is not None]

            sample_reads = {int(r["position"]): int(r["reads"]) for r in reads}
            all_samples_reads[sample] = sample_reads

        # Run correlation analysis
        correlation_scores(
            all_samples_reads,
            binsize=CORRELATIONS_BIN_SIZE,
            genome=genome,
            outputs_directory=outputs_directory,
        )
    print(f"Finished in {round(time.perf_counter()-start_time,2)} seconds")


main()
